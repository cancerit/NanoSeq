#ifndef PILEUP_CUSTOM_H
#define PILEUP_CUSTOM_H

#include <htslib/sam.h>
#include "constants.h"
#include "base.h"

typedef struct base_array_t {
    uint64_t count;
    uint64_t capacity;
    base_t *bases;
} base_array_t;

static inline int base_info_array_reset(base_array_t *a, const uint64_t capacity) {
    a->count = 0;
    if (capacity > a->capacity) {
        a->bases = (base_t*)realloc(a->bases, capacity * sizeof(base_t));
        a->capacity = capacity;
    }
    return a->bases == NULL;
}

int base_array_update(base_array_t *ba, bam1_t *read, const uint8_t min_qual);

#endif
