#include <algorithm>  // fill
#include "pileup_batch.h"
#include "utils.h"

void PileupBatch::Update(const int start, const int end) {
  this->mask.Reset({start, end});
}

void PileupBatch::EvalPos(const int pos) {
  const uint8_t mask_flag = this->mask.GetFlag(pos);
  const int is_snp = flag_is_set(mask_flag, MASK_INDEX_SNP);
  const int is_masked = flag_is_set(mask_flag, MASK_INDEX_NOISE);

  // ...
  // Have separate structures (and maybe pre-cached TSV slices) for the three sets of columns: duplex info, duplex stats, bulk stats
}
