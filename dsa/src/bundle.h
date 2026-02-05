#ifndef BUNDLE_H_
#define BUNDLE_H_

#include <stdint.h>
#include <htslib/sam.h>
#include "constants.h"
#include "read_info.h"
#include "pos_stats.h"

/// Stats common to duplex and bulk bundles
typedef struct bundle_t {
    uint64_t read_counts[RTYPE_COUNT] = {};
    int64_t asxs[RTYPE_COUNT] = {};
    int64_t nm[RTYPE_COUNT] = {};
    int64_t proper_pairs[RTYPE_COUNT] = {};

    uint64_t allele_counts[RTYPE_COUNT][ALLELE_COUNT] = {};
} bundle_t;

static inline void bundle_update(bundle_t *bundle, const read_info_t *r, const int32_t i) {
	// Stats shared with bulk bundles
    // ASSUMPTION: strand and read type have been sanitised already
    bundle->read_counts[i]++;
    bundle->asxs[i] += r->asxs;
    bundle->nm[i] += r->nm;
    bundle->proper_pairs[i] += r->proper_pair;
}

static inline double bundle_mean(const bundle_t *b, const uint64_t r_type, const int64_t x[2]) {
    return static_cast<double>(x[r_type]) / static_cast<double>(b->read_counts[r_type]);
}

static inline void pos_stats_set_common(pos_final_stats_t *s, const bundle_t *b) {
    const bool has_a = b->read_counts[RTYPE_A] != 0;
    const bool has_b = b->read_counts[RTYPE_B] != 0;

    // Mean AS - XS and proper pair counts
    if (has_a && has_b) {

        s->asxs = std::min(
            bundle_mean(b, RTYPE_A, b->asxs),
            bundle_mean(b, RTYPE_B, b->asxs));

        s->proper_pairs = std::min(
            bundle_mean(b, RTYPE_A, b->proper_pairs),
            bundle_mean(b, RTYPE_B, b->proper_pairs));

    } else if (has_a) {
        s->asxs = bundle_mean(b, RTYPE_A, b->asxs);
        s->proper_pairs = bundle_mean(b, RTYPE_A, b->proper_pairs);
    } else if (has_b) {
        s->asxs = bundle_mean(b, RTYPE_B, b->asxs);
        s->proper_pairs = bundle_mean(b, RTYPE_B, b->proper_pairs);
    } else {
        s->asxs = 0.0;
        s->proper_pairs = 0.0;
    }
}

#endif
