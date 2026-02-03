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

/* DUPLEX STATS */

typedef struct bundle_open_t {
    duplex_tag_info duplex_tag_info = {};

    int64_t asxs_accum[RTYPE_COUNT] = {0, 0};
    int64_t clip_accum[RTYPE_COUNT] = {0, 0};
    int64_t nmms_accum[RTYPE_COUNT] = {0, 0};

    int64_t rtype_ppair_counts[RTYPE_COUNT] = {0, 0};
    uint64_t rtype_read_counts[RTYPE_COUNT] = {0, 0};  // then divide ppair to get the averages
} bundle_open_t;

static inline bool bundle_open_is_empty(const bundle_open_t *b) {
    return b->rtype_read_counts[RTYPE_A] == 0 && b->rtype_read_counts[RTYPE_B] == 0;
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
    duplex_tag_info duplex_tag_info = {};
    double asxs = 0.0;
    double nm = 0.0;
    double clip = 0.0;
    double proper_pairs = 0.0;
} bundle_closed_t;

// void bundle_closed_bulk_init(bundle_closed_t *s, const bundle_open_t *b);
void bundle_closed_duplex_init(bundle_closed_t *s, const bundle_open_t *b);

typedef struct duplex_base_t {
    uint64_t counts[RTYPE_COUNT][ALLELE_COUNT] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
    uint64_t duplex_depth[STRAND_COUNT][READ_TYPE_COUNT] = {{0, 0}, {0, 0}};
    double duplex_consensus_quality_accum[RTYPE_COUNT][ALPH_LEN] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
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
    uint64_t total;
    for (int i = 0; i < RTYPE_COUNT; ++i) {
        total = 0;
        // NOTE: skipping the first allele, which is a placeholder for discarded alleles
        for (int j = 1; j < ALLELE_COUNT; ++j) {
            total += base->counts[i][j];
        }
        if (total != 0) {
            finalise_consensus_quality_scores(base->duplex_consensus_quality_accum[i]);
        }
    }
}

/* BULK STATS */

typedef struct bulk_base_open_t {
    uint64_t counts[RTYPE_COUNT][ALLELE_COUNT] = {{0, 0}, {0, 0}};

    int64_t asxs_accum[RTYPE_COUNT] = {0, 0};
    int64_t nm_accum[RTYPE_COUNT] = {0, 0};
    int64_t ppair_accum[RTYPE_COUNT] = {0, 0};

    uint64_t read_counts[RTYPE_COUNT] = {0, 0};
} bulk_base_open_t;

typedef struct bulk_base_closed_t {
    uint64_t counts[RTYPE_COUNT][ALLELE_COUNT] = {{0, 0}, {0, 0}};
    double asxs = 0.0;
    double nm = 0.0;
    double clip = 0.0;
    double proper_pairs = 0.0;
} bulk_base_closed_t;

void bulk_base_closed_init(bulk_base_closed_t *c, const bulk_base_open_t *o);

typedef struct bulk_read_info_t {
    int64_t asxs = 0;
    int64_t nm = 0;
    int64_t proper_pair = 0;
} bulk_read_info_t;

static inline void bulk_read_info_init(bulk_read_info_t *r, const bam1_t *read) {
    r->asxs = get_as_minus_xs(read);
    r->nm = get_nm(read);
    r->proper_pair = read_is_in_proper_pair(read);
}

static inline void bulk_base_update(bulk_base_open_t *b, const bulk_read_info_t *read_info, const base_t *bi, const int32_t strand) {
    // ASSUMPTION: strand has been sanitised already
    b->asxs_accum[strand] += read_info->asxs;
    b->nm_accum[strand] += read_info->nm;
    b->ppair_accum[strand] += read_info->proper_pair;
    b->read_counts[strand]++;
    b->counts[strand][bi->base]++;
}

std::string bulk_base_get_dsa_chunk(bulk_base_closed_t *bulk, const std::string pos_prefix);

/* GENOMIC POSITION STATS */

typedef struct pos_stats_t {
    bulk_base_open_t bulk_base = {};
    std::map<uint64_t, duplex_base_t> duplex_bases = {};
} pos_stats_t;

void dsa_push_row(
    std::stringstream &s,
    const bundle_closed_t *duplex_bundle,
    const duplex_base_t *duplex_base,
    const bulk_base_closed_t *bulk_base,
    const std::string pos_prefix,
    const uint8_t bundle_type
);

#endif
