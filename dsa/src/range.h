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

inline bool range_is_valid(const range_t *r) {
    return r->start >= 0 && r->end > r->start;
}

inline int32_t range_length(const range_t *r) {
    assert (range_is_valid(r));
    return r->end - r->start;
}

inline void range_clamp(range_t *r, const range_t *t) {
    assert (range_is_valid(r));
    assert (range_is_valid(t));
    r->start = std::max(r->start, t->start);
    r->end = std::min(r->end, t->end);
}

/// Grow range by two units either side to accommodate for triplet retrieval
[[nodiscard]] inline range_t range_triplet_grow(const range_t *r) {
    assert (range_is_valid(r));
    return {
        r->start <= 2 ? 0 : (r->start - 2),
        r->end + 2
    };
}

inline bool range_contains(const range_t *r, const int32_t pos) {
    assert (range_is_valid(r));
    return pos >= r->start && pos < r->end;  // for half-open range
}

typedef struct {
    range_t grange;
    int32_t tid;
} genomic_region_t;


#endif
