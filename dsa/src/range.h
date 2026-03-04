#ifndef RANGE_H_
#define RANGE_H_

#include <algorithm>
#include <cassert>
#include <stdint.h>

// 0-indexed end-exclusive aka half-open
// (as BED spec)
typedef struct {
    int32_t start;  // signed ints per htslib
    int32_t end;
} range_t ;

inline void range_is_regular(const range_t *r) {
    assert(r->start >= 0 && r->end > r->start);
}

inline int32_t range_length(const range_t *r) {
    return r->end - r->start;
}

inline void range_regularise(range_t *r) {
    if (r->end < r->start) {
        const auto t = r->start;
        r->start = r->end;
        r->end = t;
    }
}

inline void range_clamp(range_t *r, const range_t *t) {
    // ASSUMPTION: t is regularised
    range_regularise(r);
    r->start = std::max(r->start, t->start);
    r->end = std::min(r->end, t->end);
}

/// Grow range by two units either side to accommodate for triplet retrieval
inline range_t range_triplet_grow(const range_t *r) {
    return {
        r->start <= 2 ? 0 : (r->start - 2),
        r->end + 2
    };
}

inline bool range_contains(const range_t *r, const int32_t pos) {
    return pos >= r->start && pos < r->end;  // for half-open range
}

typedef struct {
    range_t grange;
    int32_t tid;
} genomic_region_t;



#endif
