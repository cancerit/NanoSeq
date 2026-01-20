#ifndef PILEUP_BATCH_H_
#define PILEUP_BATCH_H_

#include "constants.h"
#include "mask.h"
#include "mask_loader.h"
#include "range.h"
#include "ref.h"
#include "writeout.h"
#include "htslib/sam.h"
#include "static_string_builder.hpp"

class PileupBatch {
private:
    const char *contig;
    range_t range;
    StaticStringBuilder<MAX_DSA_LINE_LENGHT> ssb;

    const std::string PositionString(const char *contig, const int pos, Ref *ref, const uint8_t mask_values[MASK_COUNT]);

public:
    Mask mask;
    PileupBatch() : contig(nullptr), range({0, 0}) {};
    void Update(const char *contig, const range_t range, MaskLoader mls[2], Ref *ref);
    void Pileup(bam_mplp_t mplp, Ref *ref, WriteOut *out);
};

#endif
