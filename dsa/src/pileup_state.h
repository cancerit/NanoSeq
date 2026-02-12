#ifndef PILEUP_STATE_H_
#define PILEUP_STATE_H_

#include "aux.h"
#include "bulk_bundle.h"
#include "compressor.h"
#include "duplex_bundle.h"
#include "duplex_tag_info.h"
#include "options.h"
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

    uint64_t dsa_row_count = 0;

    FILE *debug_pos_duplexes_f = NULL;
} pileup_state_t;

#endif
