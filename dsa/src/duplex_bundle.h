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

#ifndef DUPLEX_BUNDLE_H_
#define DUPLEX_BUNDLE_H_

#include "constants.h"
#include "probs.h"
#include "bundle.h"
#include "base.h"

inline void read_info_duplex_init(read_info_t *r, const bam1_t *read) {
    read_info_init(r, read);
    r->is_5p_clipped = get_is_5p_clipped(read, r->strand);
}

typedef struct duplex_bundle_t {
    bundle_t bundle;
    int64_t clip[RTYPE_COUNT] = {};

    uint64_t duplex_depth[STRAND_COUNT][READ_TYPE_COUNT] = {};
    double duplex_consensus_quality_accum[RTYPE_COUNT][ALPH_LEN] = {};
} duplex_bundle_t;

inline int high_duplex_depth(const duplex_bundle_t *b, const int strand, const uint64_t min_dplx_depth) {
    return static_cast<int>(
        (b->duplex_depth[strand][READ_TYPE_INDEX_READ_1] >= min_dplx_depth) &&
        (b->duplex_depth[strand][READ_TYPE_INDEX_READ_2] >= min_dplx_depth));
}

inline int duplex_base_get_bundle_type(const duplex_bundle_t *b, const uint64_t min_dplx_depth) {
    // BEWARE: only use on duplex bundles (not bulk)!
    return
        (high_duplex_depth(b, STRAND_INDEX_REVERSE, min_dplx_depth) << 1) |
        (high_duplex_depth(b, STRAND_INDEX_FORWARD, min_dplx_depth) << 0);
}

inline void duplex_bundle_update(duplex_bundle_t *bundle, const read_info_t *r, const probs_t *probs, const base_t *base) {
    // Duplex-specific stats
    // NOTE: zero as default should replicate missing key (strand) in the original implementation
    int r_type = 0;
    if (r->strand != STRAND_INDEX_IGNORE) {
        r_type = RTYPES[r->strand][r->read_index];
        bundle->duplex_depth[r->strand][r->read_index]++;
    }
    bundle->clip[r_type] += r->is_5p_clipped;
    bundle_update(&bundle->bundle, r, r_type);

    if (base->base != ALLELE_INVALID) {
        bundle->bundle.allele_counts[r_type][base->base]++;
        probs_add_p_error(
            probs, base->qual, base->base,
            bundle->duplex_consensus_quality_accum[r_type]);
    }
}

inline void duplex_bundle_finalise(duplex_bundle_t *bundle, pos_final_stats_t *s) {
    const bundle_t *b = &bundle->bundle;
    pos_stats_set_common(s, b);

    // NM
    s->nm = std::max(
        bundle_mean(b, RTYPE_A, b->nm),
        bundle_mean(b, RTYPE_B, b->nm));

    // 5'-clipping
    {
        const uint64_t total = b->read_counts[RTYPE_A] + b->read_counts[RTYPE_B];
        if (total != 0) {
            s->clip = custom_round(static_cast<float>(bundle->clip[RTYPE_A] + bundle->clip[RTYPE_B]) / static_cast<float>(total));
        } else {
            s->clip = 0;
        }
    }

    // Duplex consensus base quality scores (overrides accumulator!)
    for (int i = 0; i < RTYPE_COUNT; ++i) {
        finalise_consensus_quality_scores(bundle->duplex_consensus_quality_accum[i]);
    }
}

#endif
