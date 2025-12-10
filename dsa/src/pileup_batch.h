#ifndef PILEUP_BATCH_H_
#define PILEUP_BATCH_H_

#include "constants.h"
#include "mask.h"
#include "mask_loader.h"
#include "range.h"
#include "ref.h"

class PileupBatch {
private:
    const char *contig;
    range_t range;

public:
    Mask mask;
    PileupBatch() : contig(nullptr), range({0, 0}) {};
    void Update(const char *contig, const range_t range, MaskLoader mls[2], Ref *ref);
    void EvalPos(const int pos);
    // void MultiplePileup();
};

#endif
