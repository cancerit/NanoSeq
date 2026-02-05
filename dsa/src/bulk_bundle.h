#ifndef BULK_BUNDLE_H_
#define BULK_BUNDLE_H_

#include "base.h"
#include "bundle.h"

typedef bundle_t bulk_bundle_t;

static inline void bulk_bundle_update(bulk_bundle_t *bundle, const read_info_t *r, const base_t *base) {
	bundle_update(bundle, r, r->strand);
	bundle->allele_counts[r->strand][base->base]++;
}

static inline void bulk_bundle_finalise(bulk_bundle_t *b, pos_final_stats_t *s) {
	const bool has_a = b->read_counts[RTYPE_A] != 0;
    const bool has_b = b->read_counts[RTYPE_B] != 0;
	pos_stats_set_common(s, b);

	// NM
    if (has_a && has_b) {
        // TODO: should this be max(a, b) instead (same as duplex NM)?
        s->nm = bundle_mean(b, RTYPE_A, b->nm);
    } else if (has_a) {
        s->nm = bundle_mean(b, RTYPE_A, b->nm);
    } else if (has_b) {
        s->nm = bundle_mean(b, RTYPE_B, b->nm);
    } else {
        s->nm = 0.0f;
    }
}

#endif
