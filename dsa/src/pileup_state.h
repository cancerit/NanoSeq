#ifndef PILEUP_STATE_H_
#define PILEUP_STATE_H_

#include "aux.h"
#include "bulk_bundle.h"
#include "compressor.h"
#include "duplex_bundle.h"
#include "duplex_tag_info.h"
#include "options.h"
#include "pileup_custom.h"
#include "range.h"
#include "ref.h"
#include <map>

typedef struct pos_stats_t {
    bulk_bundle_t bulk = {};
    std::map<uint64_t, duplex_bundle_t> duplexes = {};
} pos_stats_t;

typedef struct pileup_state_t {
    aux_t *bulk_aux;
    aux_t *duplex_aux;
    Ref *ref;
    Options *opts;
    GzipCompressor *compressor;

    std::map<int32_t, pos_stats_t> pos_bundles = {};
    // Genomic position -> bundle indices
    std::vector<duplex_tag_info_t> bundle_id_decoder = {};

    probs_t probs = {};  // precomputed quality score stats
    base_array_t base_buffer = {};  // buffer for usage by the CIGAR stepper

    bam1_t *read = NULL;
    uint64_t dsa_row_count = 0;
    int32_t read_start = -1;

    range_tid_t range = {};

    FILE *debug_pos_duplexes_f = NULL;
} pileup_state_t;

static void pileup_state_init(pileup_state_t *state) {
    state->read = bam_init1();

    if (base_info_array_reset(&state->base_buffer, 256)) {
        throw std::runtime_error("Failed to allocate base info array!");
    }

    probs_init(&state->probs);
}

#endif
