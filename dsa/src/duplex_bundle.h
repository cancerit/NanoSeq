#ifndef DUPLEX_BUNDLE_H_
#define DUPLEX_BUNDLE_H_

#include "constants.h"
#include "probs.h"
#include "bundle.h"
#include "base.h"

static inline void read_info_duplex_init(read_info_t *r, const bam1_t *read) {
	read_info_init(r, read);
	r->is_5p_clipped = get_is_5p_clipped(read, r->strand);
}

typedef struct duplex_bundle_t {
	bundle_t bundle;
    int64_t clip[RTYPE_COUNT] = {};

    uint64_t duplex_depth[STRAND_COUNT][READ_TYPE_COUNT] = {};
    double duplex_consensus_quality_accum[RTYPE_COUNT][ALPH_LEN] = {};
} duplex_bundle_t;

static inline int high_duplex_depth(const duplex_bundle_t *b, const int strand, const uint64_t min_dplx_depth) {
    return static_cast<int>(
        (b->duplex_depth[strand][READ_TYPE_INDEX_READ_1] >= min_dplx_depth) &&
        (b->duplex_depth[strand][READ_TYPE_INDEX_READ_2] >= min_dplx_depth));
}

static inline int duplex_base_get_bundle_type(const duplex_bundle_t *b, const uint64_t min_dplx_depth) {
    // BEWARE: only use on duplex bundles (not bulk)!
    return
        (high_duplex_depth(b, STRAND_INDEX_REVERSE, min_dplx_depth) << 1) |
        (high_duplex_depth(b, STRAND_INDEX_FORWARD, min_dplx_depth) << 0);
}

static inline void duplex_bundle_update(duplex_bundle_t *bundle, const read_info_t *r, const probs_t *probs, const base_t *base) {
    // Duplex-specific stats
    // NOTE: zero as default should replicate missing key (strand) in the original implementation
    int r_type = 0;
    if (r->strand != STRAND_INDEX_IGNORE) {
        r_type = RTYPES[r->strand][r->read_index];
        bundle->duplex_depth[r->strand][r->read_index]++;
    }
    bundle->clip[r_type] += r->is_5p_clipped;
    bundle_update(&bundle->bundle, r, r_type);

    bundle->bundle.allele_counts[r_type][base->base]++;
    probs_add_p_error(
    	probs, base->qual, base->base,
     	bundle->duplex_consensus_quality_accum[r_type]);
}

static inline void duplex_bundle_finalise(duplex_bundle_t *bundle, pos_final_stats_t *s) {
	const bundle_t *b = &bundle->bundle;
	const bool has_a = b->read_counts[RTYPE_A] != 0;
    const bool has_b = b->read_counts[RTYPE_B] != 0;
	pos_stats_set_common(s, b);

	// NM
    if (has_a && has_b) {
        s->nm = std::max(
            bundle_mean(b, RTYPE_A, b->nm),
            bundle_mean(b, RTYPE_B, b->nm));
    } else if (has_a) {
        s->nm = bundle_mean(b, RTYPE_A, b->nm);
    } else if (has_b) {
        s->nm = bundle_mean(b, RTYPE_B, b->nm);
    } else {
        s->nm = 0.0f;
    }

    // 5'-clipping
    {
	    const uint64_t total = b->read_counts[RTYPE_A] + b->read_counts[RTYPE_B];
	    if (total != 0) {
	        s->clip = static_cast<double>(bundle->clip[RTYPE_A] + bundle->clip[RTYPE_B]) / static_cast<double>(total);
	    } else {
	        s->clip = 0.0;
	    }
    }

    // Duplex consensus base quality scores (overrides accumulator!)
    {
	    uint64_t total;
	    for (int i = 0; i < RTYPE_COUNT; ++i) {
	        total = 0;
	        // NOTE: skipping the first allele, which is a placeholder for discarded alleles
	        for (int j = 1; j < ALLELE_COUNT; ++j) {
	            total += bundle->bundle.allele_counts[i][j];
	        }
	        if (total != 0) {
	            finalise_consensus_quality_scores(bundle->duplex_consensus_quality_accum[i]);
	        }
	    }
    }
}

#endif
