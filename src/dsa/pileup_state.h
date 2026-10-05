/*########## LICENCE ##########
# Copyright (c) 2022, 2025, 2026 Genome Research Ltd
#
# Authors: Luca Barbon <lb29@sanger.ac.uk>, Alex Byrne <ab63@sanger.ac.uk>
#
# This file is part of NanoSeq.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
##########################*/

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
#include <sstream>

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

    probs_t probs = {};  // precomputed quality score stats
    base_array_t base_buffer = {};  // buffer for usage by the CIGAR stepper

    bam1_t *read = NULL;

    // To clear between batches (pileup_state_reset)
    genomic_region_t range = {};
    uint64_t dsa_row_count = 0;
    int32_t read_start = -1;
    std::map<int32_t, pos_stats_t> pos_bundles = {};
    // Genomic position -> bundle indices
    std::vector<duplex_tag_info_t> bundle_id_decoder = {};
    std::stringstream dsa_uncompressed_stream = {};

    FILE *debug_pos_duplexes_f = NULL;
} pileup_state_t;

inline void pileup_state_init(pileup_state_t *state) {
    state->read = bam_init1();

    if (base_info_array_init(&state->base_buffer, 256)) {
        throw std::runtime_error("Failed to allocate base info array!");
    }

    probs_init(&state->probs);
}

inline void pileup_state_reset_dsa_stream(pileup_state_t *state) {
    state->dsa_uncompressed_stream.str("");
    state->dsa_uncompressed_stream.clear();
}

inline void pileup_state_reset(pileup_state_t *state, const genomic_region_t *r) {
    const auto gr = r->grange;
    state->range.grange.start = gr.start;
    state->range.grange.end = gr.end;
    state->range.tid = r->tid;
    state->read_start = -1;
    state->dsa_row_count = 0;
    state->pos_bundles.clear();
    state->bundle_id_decoder.clear();
    pileup_state_reset_dsa_stream(state);
}

#endif
