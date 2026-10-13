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

#ifndef BULK_BUNDLE_H_
#define BULK_BUNDLE_H_

#include "base.h"
#include "bundle.h"

typedef bundle_t bulk_bundle_t;

inline void bulk_bundle_update(bulk_bundle_t *bundle, const read_info_t *r, const base_t *base, const uint8_t min_qual) {
    if (base->base == ALLELE_DEL || base->qual >= min_qual) {
        bundle_update(bundle, r, r->strand);
        if (base->base != ALLELE_INVALID) {
            bundle->allele_counts[r->strand][base->base]++;
        }
    }
}

inline void bulk_bundle_finalise(bulk_bundle_t *b, pos_final_stats_t *s) {
    const bool has_a = b->read_counts[RTYPE_A] != 0;
    const bool has_b = b->read_counts[RTYPE_B] != 0;
    pos_stats_set_common(s, b);

    // NM
    if (has_a && has_b) {
        // TODO: should this be max(a, b) instead (same as duplex NM)?
        s->nm = bundle_mean(b, RTYPE_A, b->nm);
    } else if (has_a) {
        s->nm = bundle_mean(b, RTYPE_A, b->nm);
    } else if (has_b) {
        s->nm = bundle_mean(b, RTYPE_B, b->nm);
    } else {
        s->nm = 0.0f;
    }
}

#endif
