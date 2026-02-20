#ifndef PILEUP_CUSTOM_H
#define PILEUP_CUSTOM_H

#include <htslib/sam.h>
#include "constants.h"
#include "base.h"

typedef struct base_array_t {
    int32_t start;  // Alignment start position, used to find anchors for indels on the first base
    uint64_t count;
    uint64_t capacity;
    base_t *bases;
} base_array_t;

static inline int base_info_array_reset(base_array_t *a, const uint64_t capacity, const int32_t start) {
    a->count = 0;
    a->start = start;
    if (capacity > a->capacity) {
        a->bases = (base_t*)realloc(a->bases, capacity * sizeof(base_t));
        a->capacity = capacity;
    }
    return a->bases == NULL;
}

static inline int base_info_array_init(base_array_t *a, const uint64_t capacity) {
    return base_info_array_reset(a, capacity, -1);
}

static inline base_t *base_info_array_get_next(const base_array_t *a) {
    return &a->bases[a->count];
}

static inline base_t *base_info_get_last(const base_array_t *a) {
    return a->count != 0 ? &a->bases[a->count - 1] : NULL;
}

int base_array_update(base_array_t *ba, bam1_t *read);

#endif
