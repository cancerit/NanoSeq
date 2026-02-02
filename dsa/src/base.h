#ifndef BASE_H_
#define BASE_H_

#include <stdint.h>

typedef struct base_t {
    int16_t read_pos;
    int32_t aln_pos;
    uint8_t base;
    uint8_t qual;
} base_t;

#endif
