#ifndef RANGE_H_
#define RANGE_H_

#include <stdint.h>

typedef struct {
  int32_t start;
  int32_t end;
} range_t;

static inline int32_t range_length(const range_t *r) {
  return r->end - r->start;
}

static inline const range_t range_grow(const range_t *r) {
  return {
    r->start <= 1 ? 0 : (r->start - 1),
    r->end + 1
  };
}

typedef struct {
  int32_t tid;
  int32_t start;
  int32_t end;
} range_tid_t;

#endif
