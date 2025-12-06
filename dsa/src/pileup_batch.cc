#include <algorithm>  // fill
#include "pileup_batch.h"
#include "utils.h"

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

void PileupBatch::EvalPos(const int pos) {
  const uint8_t mask_flag = this->mask[pos - this->start];
  const int is_snp = flag_is_set(mask_flag, MASK_INDEX_SNP);
  const int is_masked = flag_is_set(mask_flag, MASK_INDEX_NOISE);

  // ...
  // Have separate structures (and maybe pre-cached TSV slices) for the three sets of columns: duplex info, duplex stats, bulk stats
}
