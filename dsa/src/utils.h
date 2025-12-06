#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <numeric>
#include <cmath>
#include "constants.h"

#define flag_is_set(a,b) (((a) & (b)) != 0)

static inline void get_region(const char *contig, const int start, const int end, char region[MAX_REGION_STR_LENGTH]) {
  snprintf(region, MAX_REGION_STR_LENGTH, "%s:%d-%d", contig, start, end);
}

static inline int custom_round(const double d) {
  const int n = std::round(d);
  // TODO: verify how this can occur
  return (n == -0) ? 0 : n;
}

static inline int64_t vector_sum(std::vector<int> v) {
  return std::reduce(v.begin(), v.end(), int64_t{0});
}

static inline float vector_mean(std::vector<int> v) {
  const size_t n = v.size();
  if (n == 0) {
    // TODO: reconsider
    return 0.0f;
  }
  return static_cast<float>(
    static_cast<double>(vector_sum(v)) /
    static_cast<double>(n));
}

static float vector_pair_mean(std::vector<int> v1, std::vector<int> v2) {
  const size_t n = v1.size() + v2.size();
  if (n == 0) {
    // TODO: reconsider
    return 0.0f;
  }
  const int64_t total = vector_sum(v1) + vector_sum(v2);
  return static_cast<float>(
    static_cast<double>(total) /
    static_cast<double>(n));
}

#endif
