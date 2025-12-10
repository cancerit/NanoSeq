#ifndef RANGE_H_
#define RANGE_H_

#include <stdint.h>

typedef struct {
    int32_t start;
    int32_t end;
} range_t;

static inline void range_validate(const range_t *r) {
    assert(r->end > r->start && r->end > 0);
}

static inline int32_t range_length(const range_t *r) {
    return r->end - r->start;
}

static inline void range_regularise(range_t *r) {
    if (r->end < r->start) {
        const int32_t t = r->start;
        r->start = r->end;
        r->end = t;
    }
}

static inline void range_clamp(range_t *r, const range_t *t) {
    // ASSUMPTION: t is regularised
    range_regularise(r);
    if (r->start < t->start) {
        r->start = t->start;
    }
    if (r->end > t->end) {
        r->end = t->end;
    }
}

static inline const range_t range_grow(const range_t *r) {
    return {
        r->start <= 1 ? 0 : (r->start - 1),
        r->end + 1
    };
}

typedef struct {
    int32_t start;
    int32_t end;
    int32_t tid;
    int32_t __padding;
} range_tid_t;

#endif
