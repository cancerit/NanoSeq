#ifndef BUNDLE_H_
#define BUNDLE_H_

#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "constants.h"
#include "duplex_tag_info.h"
#include "base.h"
#include "probs.h"
#include "read_info.h"

typedef struct bundle_open_t {
    int64_t asxs_accum[RTYPE_COUNT] = {0, 0};
    int64_t clip_accum[RTYPE_COUNT] = {0, 0};
    int64_t nmms_accum[RTYPE_COUNT] = {0, 0};

    int64_t rtype_ppair_counts[RTYPE_COUNT] = {0, 0};
    uint64_t rtype_read_counts[RTYPE_COUNT] = {0, 0};  // then divide ppair to get the averages
} bundle_open_t;

static inline bool bundle_open_is_empty(const bundle_open_t *b) {
    return b->rtype_read_counts[RTYPE_A] == 0 && b->rtype_read_counts[RTYPE_B] == 0;
}

static inline void bundle_open_bulk_update(bundle_open_t *bundle, const bam1_t *read, const int32_t strand) {
    // ASSUMPTION: strand has been sanitised already
    bundle->asxs_accum[strand] += get_as_minus_xs(read);
    bundle->nmms_accum[strand] += get_nm(read);
    bundle->rtype_ppair_counts[strand] += read_is_in_proper_pair(read);
    bundle->rtype_read_counts[strand]++;
}

static inline int bundle_open_duplex_update(bundle_open_t *bundle, const bam1_t *read) {
    const int strand = get_strand_index(read);
    const int read_index = get_read_type_index(read);

    // Duplex-specific stats
    // NOTE: zero as default should replicate missing key (strand) in the original implementation
    int r_type = 0;
    if (strand != STRAND_INDEX_IGNORE) {
        r_type = RTYPES[strand][read_index];
        // bundle->dplx_depth[strand][read_index]++;
    }
    bundle->clip_accum[r_type] += get_is_5p_clipped(read, strand);

    // Stats shared with bulk bundles
    // ASSUMPTION: strand and read type have been sanitised already
    bundle->asxs_accum[r_type] += get_as_minus_xs(read);
    bundle->nmms_accum[r_type] += get_nm(read);
    bundle->rtype_ppair_counts[r_type] += read_is_in_proper_pair(read);
    bundle->rtype_read_counts[r_type]++;

    return r_type;
}

typedef struct bundle_closed_t {
    // duplex_tag_info duplex_tag_info = {};
	// double consensus_qualities[RTYPE_COUNT][ALPH_LEN] = {};
    double asxs = 0.0;
    double nm = 0.0;
    double clip = 0.0;
    double proper_pairs = 0.0;
} bundle_closed_t;

void bundle_closed_bulk_init(bundle_closed_t *s, const bundle_open_t *b);
void bundle_closed_duplex_init(bundle_closed_t *s, const bundle_open_t *b);

/* PER-BASE DUPLEX STATS */

typedef struct duplex_base_t {
    // Duplex-only
    double duplex_consensus_quality_accum[RTYPE_COUNT][ALPH_LEN];
    uint64_t duplex_depth[STRAND_COUNT][READ_TYPE_COUNT] = {{0, 0}, {0, 0}};
    // Duplex and bulk
    uint64_t counts[BUNDLE_TYPES_COUNT][RTYPE_COUNT][ALLELE_COUNT];
} duplex_base_t;

static inline int high_duplex_depth(const duplex_base_t *b, const int strand, const uint64_t min_dplx_depth) {
    return static_cast<int>(
        (b->duplex_depth[strand][READ_TYPE_INDEX_READ_1] >= min_dplx_depth) &&
        (b->duplex_depth[strand][READ_TYPE_INDEX_READ_2] >= min_dplx_depth));
}

static inline int duplex_base_get_bundle_type(const duplex_base_t *b, const uint64_t min_dplx_depth) {
    // BEWARE: only use on duplex bundles (not bulk)!
    return
        (high_duplex_depth(b, STRAND_INDEX_REVERSE, min_dplx_depth) << 1) |
        (high_duplex_depth(b, STRAND_INDEX_FORWARD, min_dplx_depth) << 0);
}

// Calculate the duplex consensus quality scores (overrides accumulator)
static inline void duplex_base_finalise(duplex_base_t *base) {
    for (int i = 0; i < RTYPE_COUNT; ++i) {
        if (base->counts[DUPLEX_INDEX][i] != 0) {
            finalise_consensus_quality_scores(base->duplex_consensus_quality_accum[i]);
        }
    }
}

/* BUNDLE PAIR */

typedef struct bundle_closed_pair_t {
    duplex_tag_info duplex_tag_info;
    bundle_closed_t bundles[BUNDLE_TYPES_COUNT];
} bundle_closed_pair_t;

std::string bundle_closed_pair_to_dsa_row(const bundle_closed_pair_t *bp, const std::string pos_prefix, const duplex_base_t *base, const uint8_t bundle_type);

#endif
