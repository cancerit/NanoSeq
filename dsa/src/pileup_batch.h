#ifndef PILEUP_BATCH_H_
#define PILEUP_BATCH_H_

#include "constants.h"
#include "mask.h"
#include "mask_loader.h"
#include "range.h"
#include "htslib/sam.h"
#include "pileup_state.h"

class PileupBatch {
private:
    const char *contig = NULL;
    int32_t tid = -1;
    range_t range = {0, 0};

    const std::string PositionString(const char *contig, const int pos, Ref *ref, const uint8_t mask_values[MASK_COUNT]);

public:
    Mask mask;
    PileupBatch() : contig(nullptr), tid(-1), range({0, 0}) {};
    void PileupBulk(const Options *opts, pileup_state_t *state);
    void PileupDumpPosition(const Options *opts, pileup_state_t *state, const int32_t pos);
    void Update(const char *contig, const range_tid_t range, MaskLoader mls[2], Ref *ref);
    void Pileup(pileup_state_t *state);
};

#endif
