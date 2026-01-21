#include <algorithm>  // fill
#include "pileup_batch.h"
#include "utils.h"
#include <iostream>
#include <assert.h>
#include "read_bundler.h"

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

    // DEBUG ONLY!
    /*
    std::cerr << "[" << ref->ToString() << "]" << std::endl;
    std::cerr << "<" << ref->GetTripletAround(range.start) << ">" << std::endl;
    std::cerr << "<" << ref->GetTripletAround(range.start + 1) << ">" << std::endl;
    std::cerr << "<" << ref->GetTripletAround(range.end) << ">" << std::endl;
    */
}

const std::string PileupBatch::PositionString(const char *contig, const int pos, Ref *ref, const uint8_t mask_values[MASK_COUNT]) {
    const std::string_view ctx = ref->GetTripletAround(pos);

    std::stringstream ss;
    ss << contig;
    ss << "\t";
    ss << pos;
    ss << "\t";
    ss << pos + 1;
    ss << "\t";
    ss << ctx;
    ss << "\t";
    ss << static_cast<int>(mask_values[MASK_INDEX_SNP]);
    ss << "\t";
    ss << static_cast<int>(mask_values[MASK_INDEX_NOISE]);

    return ss.str();
}

void PileupBatch::Pileup(bam_mplp_t mplp, Ref *ref, WriteOut *out) {
    int pos, tid;
    const Options *opts = out->opts;
    int n_plp[BUNDLE_TYPES_COUNT];
    std::vector<const bam_pileup1_t *> plps[BUNDLE_TYPES_COUNT];
    const bam_pileup1_t *plp[BUNDLE_TYPES_COUNT];
    uint8_t mask_flag = 0;
    uint8_t mask_values[MASK_COUNT] = {0, 0};
    std::string posn;
    ReadBundler rb;

    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {
        plps[BUNDLE_TYPE_BULK].clear();
        plps[BUNDLE_TYPE_DUPLEX].clear();

        // TODO: verify end inclusiveness convention!
        //  Originally: ((pos >= opts->beg) && (pos <= opts->end))
        if (pos < range.start) {
            continue;
        } else if (pos >= range.end) {
            break;
        }

        // Pileup
        for (int i = 0; i < BUNDLE_TYPES_COUNT; i++) {
            for (int j = 0; j < n_plp[i]; ++j) {
                plps[i].push_back(plp[i] + j);
            }
        }

        // Bundle reads
        bundles dplx = rb.DplxBundles(pos, opts->offset, opts->min_dplx_depth, plps[BUNDLE_TYPE_DUPLEX]);
        if (dplx.size() == 0) {
            continue;
        }
        bundle bulk = rb.BulkBundle(plps[BUNDLE_TYPE_BULK], opts->min_base_quality);

        // Generate DSA table row prefix
        mask_flag = this->mask.GetFlag(pos);
        mask_values[MASK_INDEX_SNP] = flag_is_set(mask_flag, MASK_FLAG_SNP);
        mask_values[MASK_INDEX_NOISE] = flag_is_set(mask_flag, MASK_FLAG_NOISE);

        // TODO: avoid string reallocation!
        posn = PositionString(contig, pos, ref, mask_values);

        // Push DSA table rows to compressor
        out->WriteRows(bulk, dplx, posn);
    }
}
