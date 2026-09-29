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

#ifndef RANGE_H_
#define RANGE_H_

#include <algorithm>
#include <cassert>
#include <stdint.h>

// 0-indexed end-exclusive aka half-open
// (as BED spec)
typedef struct {
    int32_t start;  // signed ints per htslib
    int32_t end;
} range_t ;

inline bool range_is_valid(const range_t *r) {
    return r->start >= 0 && r->end > r->start;
}

inline int32_t range_length(const range_t *r) {
    assert (range_is_valid(r));
    return r->end - r->start;
}

inline void range_clamp(range_t *r, const range_t *t) {
    assert (range_is_valid(r));
    assert (range_is_valid(t));
    r->start = std::max(r->start, t->start);
    r->end = std::min(r->end, t->end);
}

/// Grow range by two units either side to accommodate for triplet retrieval
[[nodiscard]] inline range_t range_triplet_grow(const range_t *r) {
    assert (range_is_valid(r));
    return {
        r->start <= 2 ? 0 : (r->start - 2),
        r->end + 2
    };
}

inline bool range_contains(const range_t *r, const int32_t pos) {
    assert (range_is_valid(r));
    return pos >= r->start && pos < r->end;  // for half-open range
}

typedef struct {
    range_t grange;
    int32_t tid;
} genomic_region_t;


#endif
