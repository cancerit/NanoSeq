#include <algorithm>  // fill
#include "pileup_batch.h"
#include "utils.h"
#include <iostream>
#include <assert.h>

void PileupBatch::Update(const char *contig, const range_t range, MaskLoader mls[2], Ref *ref) {
    assert(range_length(&range) > 0);
    this->contig = contig;
    this->range = range;

    // Load masks
    this->mask.Reset(this->range);
    for (int i = 0; i < 2; ++i) {
        mls[i].LoadMask(this->contig, this->range.start, this->range.end, this->mask);
    }
    std::cerr << std::format("Masked positions: {}\n", mask.CountBytesSet());

    // Load reference sequence
    const range_t ref_range = range_grow(&range);
    std::cerr << std::format("SLICE: {}:{}-{}\n", contig, range.start, range.end);
    std::cerr << std::format("REF: {}:{}-{}\n", contig, ref_range.start, ref_range.end);
    ref->Fetch(contig, ref_range);
    std::cerr << "[" << ref->ToString() << "]" << std::endl;
    std::cerr << "<" << ref->GetTripletAround(range.start) << ">" << std::endl;
    std::cerr << "<" << ref->GetTripletAround(range.start + 1) << ">" << std::endl;
    std::cerr << "<" << ref->GetTripletAround(range.end) << ">" << std::endl;
}

void PileupBatch::EvalPos(const int pos) {
    const uint8_t mask_flag = this->mask.GetFlag(pos);
    const int is_snp = flag_is_set(mask_flag, MASK_INDEX_SNP);
    const int is_masked = flag_is_set(mask_flag, MASK_INDEX_NOISE);

  // ...
  // Have separate structures (and maybe pre-cached TSV slices) for the three sets of columns: duplex info, duplex stats, bulk stats
}
