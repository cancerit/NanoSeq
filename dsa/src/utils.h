#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <numeric>
#include <cmath>
#include "constants.h"

#define flag_is_set(a,b) (((a) & (b)) != 0)

static inline void upper(std::string &str) {
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

static inline void get_region(const char *contig, const int start, const int end, char region[MAX_REGION_STR_LENGTH]) {
  snprintf(region, MAX_REGION_STR_LENGTH, "%s:%d-%d", contig, start, end);
}

static inline int custom_round(const double d) {
  const int n = std::round(d);
  // TODO: verify how this can occur
  return (n == -0) ? 0 : n;
}

#endif
