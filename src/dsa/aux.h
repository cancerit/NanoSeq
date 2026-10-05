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

#ifndef AUX_H_
#define AUX_H_

#include "range.h"
#include <htslib/sam.h>

typedef struct {

    // Set at startup
    htsFile *fp;
    sam_hdr_t *head;
    hts_idx_t *idx;
    int min_mapQ;

    // Set per genomic range
    genomic_region_t range;
    hts_itr_t *iter;
    uint64_t iterations;

} aux_t;

inline void aux_init(aux_t *a, const char *fp) {
    a->fp = NULL;
    a->head = NULL;
    a->idx = NULL;
    a->iter = NULL;
    a->iterations = 0;

    a->fp = hts_open(fp, "r");
    assert(a->fp);
    a->head = sam_hdr_read(a->fp);
    assert(a->head);
    a->idx = sam_index_load(a->fp, fp);
    assert(a->idx);
}

inline void aux_duplex_init(aux_t *a, const char *fp, const int min_map_q) {
    aux_init(a, fp);
    a->min_mapQ = min_map_q;
}

inline void aux_bulk_init(aux_t *a, const char *fp) {
    aux_init(a, fp);
    a->min_mapQ = 0;
}

inline void aux_reset(aux_t *a) {
    // Iterator
    if (a->iter != NULL) {
        sam_itr_destroy(a->iter);
        a->iter = NULL;
    }

    // Stats
    a->iterations = 0;

    // Range
    a->range.tid = -1;
    a->range.grange.start = 0;
    a->range.grange.end = 0;
}

inline int aux_set_iterator(aux_t *a, const genomic_region_t r) {
    aux_reset(a);
    a->range = r;
    a->iter = sam_itr_queryi(a->idx, r.tid, r.grange.start, r.grange.end);
    if (a->iter == NULL) {
        fprintf(stderr, "Failed to initialise iterator!\n");
        return 1;
    }
    return 0;
}

inline int aux_iter(aux_t *a, bam1_t *b) {
    a->iterations++;
    return sam_itr_next(a->fp, a->iter, b);
}

#endif  // AUX_H_
