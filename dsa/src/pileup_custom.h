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

#ifndef PILEUP_CUSTOM_H
#define PILEUP_CUSTOM_H

#include <htslib/sam.h>
#include "base.h"

typedef struct base_array_t {
    int32_t start;  // Alignment start position, used to find anchors for indels on the first base
    uint64_t count;
    uint64_t capacity;
    base_t *bases;
} base_array_t;

inline int base_info_array_reset(base_array_t *a, const uint64_t capacity, const int32_t start) {
    a->count = 0;
    a->start = start;
    if (capacity > a->capacity) {
        a->bases = static_cast<base_t*>(realloc(a->bases, capacity * sizeof(base_t)));
        a->capacity = capacity;
    }
    return a->bases == NULL;
}

inline int base_info_array_init(base_array_t *a, const uint64_t capacity) {
    return base_info_array_reset(a, capacity, -1);
}

inline base_t *base_info_array_get_next(const base_array_t *a) {
    return &a->bases[a->count];
}

inline base_t *base_info_get_last(const base_array_t *a) {
    return a->count != 0 ? &a->bases[a->count - 1] : NULL;
}

int base_array_update(base_array_t *ba, bam1_t *read);

#endif
