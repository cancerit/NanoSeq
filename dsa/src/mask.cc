#include <algorithm>  // fill
#include <assert.h>
#include "mask.h"

void Mask::Update(const int32_t start, const int32_t end, const uint8_t flag) {
  const int32_t a = start - this->start;
  if (end == start) {
    this->mask[a] |= flag;
  } else {
    for (int i = a; i <= (end - start); ++i) {
      this->mask[i] |= flag;
    }
  }
}

void Mask::Reset(const int32_t start, const int32_t end) {
  const int32_t m = end - start;

  const int32_t n = this->end - this->start;
  assert(this->mask.size() == n);

  // Expand the mask (if necessary) and reset it to zero
  if (m > n) {
    this->mask.resize(m);
  }
  std::fill(this->mask.begin(), this->mask.end(), 0);
}
