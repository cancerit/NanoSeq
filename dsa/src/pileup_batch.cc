#include <algorithm>  // fill
#include "pileup_batch.h"

void PileupBatch::Update(const int start, const int end) {
  const int32_t m = end - start;

  const int32_t n = this->end - this->start;
  assert(this->mask.size() == n);

  // Expand the mask (if necessary) and reset it to zero
  if (m > n) {
    this->mask.resize(m);
  }
  std::fill(this->mask.begin(), this->mask.end(), 0);
}
