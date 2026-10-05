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

#ifndef BUNDLE_H_
#define BUNDLE_H_

#include <stdint.h>
#include <htslib/sam.h>
#include "constants.h"
#include "read_info.h"
#include "pos_stats.h"
#include "utils.h"

/// Stats common to duplex and bulk bundles
typedef struct bundle_t {
    uint64_t read_counts[RTYPE_COUNT] = {};
    int64_t asxs[RTYPE_COUNT] = {};
    int64_t nm[RTYPE_COUNT] = {};
    int64_t proper_pairs[RTYPE_COUNT] = {};

    uint64_t allele_counts[RTYPE_COUNT][ALLELE_COUNT] = {};
} bundle_t;

inline void bundle_update(bundle_t *bundle, const read_info_t *r, const int32_t i) {
	// Stats shared with bulk bundles
    // ASSUMPTION: strand and read type have been sanitised already
    bundle->read_counts[i]++;
    bundle->asxs[i] += r->asxs;
    bundle->nm[i] += r->nm;
    bundle->proper_pairs[i] += r->proper_pair;
}

inline float bundle_mean(const bundle_t *b, const uint64_t r_type, const int64_t x[2]) {
    return static_cast<float>(x[r_type]) / static_cast<float>(b->read_counts[r_type]);
}

inline void pos_stats_set_common(pos_final_stats_t *s, const bundle_t *b) {
    const bool has_a = b->read_counts[RTYPE_A] != 0;
    const bool has_b = b->read_counts[RTYPE_B] != 0;

    // Mean AS - XS and proper pair counts
    if (has_a && has_b) {

        s->asxs = custom_round(std::min(
            bundle_mean(b, RTYPE_A, b->asxs),
            bundle_mean(b, RTYPE_B, b->asxs)));

        s->proper_pairs = custom_round(std::min(
            bundle_mean(b, RTYPE_A, b->proper_pairs),
            bundle_mean(b, RTYPE_B, b->proper_pairs)));

    } else if (has_a) {
        s->asxs = custom_round(bundle_mean(b, RTYPE_A, b->asxs));
        s->proper_pairs = custom_round(bundle_mean(b, RTYPE_A, b->proper_pairs));
    } else if (has_b) {
        s->asxs = custom_round(bundle_mean(b, RTYPE_B, b->asxs));
        s->proper_pairs = custom_round(bundle_mean(b, RTYPE_B, b->proper_pairs));
    } else {
        s->asxs = 0;
        s->proper_pairs = 0;
    }
}

#endif
