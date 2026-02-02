#ifndef PILEUP_BATCH_H_
#define PILEUP_BATCH_H_

#include "constants.h"
#include "mask.h"
#include "mask_loader.h"
#include "range.h"
#include "ref.h"
// #include "writeout.h"
#include "htslib/sam.h"
#include "static_string_builder.hpp"
#include "aux.h"
#include "compressor.h"
#include "options.h"

class PileupBatch {
private:
    const char *contig;
    int32_t tid;
    range_t range;
    StaticStringBuilder<MAX_DSA_LINE_LENGHT> ssb;

    const std::string PositionString(const char *contig, const int pos, Ref *ref, const uint8_t mask_values[MASK_COUNT]);

public:
    Mask mask;
    PileupBatch() : contig(nullptr), tid(-1), range({0, 0}) {};
    void Update(const char *contig, const range_tid_t range, MaskLoader mls[2], Ref *ref);
    void Pileup(aux_t **data, Ref *ref, const Options *opts, GzipCompressor *compressor);
    // void Pileup(aux_t **data, Ref *ref, WriteOut *out);
    // void Pileup(bam_mplp_t mplp, Ref *ref, WriteOut *out);
};

#endif
