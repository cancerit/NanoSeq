#include <algorithm>  // fill
#include <assert.h>
#include "mask.h"

void Mask::Update(const range_t range, const uint8_t flag) {
  const int32_t a = range.start - this->range.start;
  if (range.end == range.start) {
    this->mask[a] |= flag;
  } else {
    for (int i = a; i <= range_length(&range); ++i) {
      this->mask[i] |= flag;
    }
  }
}

void Mask::Reset(const range_t range) {
  const int32_t m = range_length(&range);
  const int32_t n = range_length(&this->range);
  assert(this->mask.size() == n);

  // Expand the mask (if necessary) and reset it to zero
  if (m > n) {
    this->mask.resize(m);
  }
  std::fill(this->mask.begin(), this->mask.end(), 0);
}
