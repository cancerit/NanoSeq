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

#include <assert.h>

#include "pileup_custom.h"
#include "constants.h"

#define count_bits(x) __builtin_popcount(x)

// Lookup tables for CIGAR operations (indexed by operation code)
// cigar_consumes_query[op] = 1 if operation consumes query sequence
static const uint8_t cigar_consumes_query[16] = {
    [BAM_CMATCH] = 1,
    [BAM_CINS] = 1,
    [BAM_CSOFT_CLIP] = 1,
    [BAM_CEQUAL] = 1,
    [BAM_CDIFF] = 1
};
// cigar_consumes_ref[op] = 1 if operation consumes reference sequence
static const uint8_t cigar_consumes_ref[16] = {
    [BAM_CMATCH] = 1,
    [BAM_CDEL] = 1,
    [BAM_CREF_SKIP] = 1,
    [BAM_CEQUAL] = 1,
    [BAM_CDIFF] = 1
};
// cigar_maps_to_ref[op] = 1 if operation creates actual alignment to reference
static const uint8_t cigar_maps_to_ref[16] = {
    [BAM_CMATCH] = 1,
    [BAM_CEQUAL] = 1,
    [BAM_CDIFF] = 1
};

static inline void base_mark_as_indel(base_t *b) {
    b->base = ALLELE_DEL;
    b->qual = 0;  // NOTE: this is probably unnecessary
}

static inline void patch_anchor_base(base_array_t *a, const int32_t pos) {
    /*
    Mark the anchor base of an indel as a deletion, matching the original behaviour:

    if ((p->is_del) || (p->indel != 0)) {
        return std::make_pair('i', -1);
    }
    */

    // NOTE: if the anchor base is non-canonical, it will not be present;
    //  therefore, the alignment position needs to be tested before patching.

    // Do not create any base before the start (to verify)
    if (pos == a->start) {
        return;
    }

    base_t *b = base_info_get_last(a);
    if (b == NULL || b->aln_pos != pos - 1) {

        // Create anchor
        // NOTE: verify under which conditions this is possible now that non-canonical and low-quality bases are still pushed to the array of bases
        b = base_info_array_get_next(a);
        b->aln_pos = pos - 1;
        base_mark_as_indel(b);
        a->count++;

    } else {

        // Patch anchor (may override an indel allele, but won't change the result)
        base_mark_as_indel(b);

    }
}

int base_array_update(base_array_t *ba, bam1_t *read) {
    const bam1_core_t *c = &read->core;
    ba->start = static_cast<int32_t>(c->pos);

    // Expand output array if necessary
    const int32_t read_length = c->l_qseq;
    // NOTE: verify soft-clipping implications on read start...
    if (base_info_array_reset(ba, static_cast<uint64_t>(read_length), c->pos)) {
        fprintf(stderr, "Failed to allocate base info array!\n");
        return 1;
    }

    uint32_t *cigar = bam_get_cigar(read);

    int32_t query_pos = 0;
    int32_t ref_offset = 0;
    const uint8_t *seq =  bam_get_seq(read);
    const uint8_t *qual = bam_get_qual(read);

    uint8_t op;
    int32_t len;
    uint8_t bam_nt;

    int32_t query_pos_end;
    int32_t ref_offset_end;

    base_t *bi = NULL;

    for (uint32_t i = 0; i < c->n_cigar; ++i) {
        op =  cigar[i] & BAM_CIGAR_MASK;
        len = cigar[i] >> BAM_CIGAR_SHIFT;
        if (cigar_maps_to_ref[op]) {
            query_pos_end = query_pos + len;
            for (; query_pos < query_pos_end; ++query_pos) {
                bam_nt = bam_seqi(seq, query_pos);

                // A. Push canonical or invalid base
                bi = base_info_array_get_next(ba);
                bi->aln_pos = c->pos + ref_offset;
                bi->qual = qual[query_pos];
                bi->base = count_bits(bam_nt) == 1 ?
                    canonical_nt16_minus_one_to_allele[bam_nt - 1] :
                    ALLELE_INVALID;
                ba->count++;

                ref_offset++;
            }
        } else if (op == BAM_CDEL) {
            // ANCHOR BASE PATCH (ONLY TO MATCH CURRENT OUTPUT!)
            patch_anchor_base(ba, c->pos + ref_offset);

            // Consider that in samtools pileup skips are marked as deletions (bug)!

            // NOTE: matching the bam_plp_* behaviour by ignoring non-reference bases
            ref_offset_end = ref_offset + len;
            for (; ref_offset < ref_offset_end; ++ref_offset) {

                // B. Push deletion
                bi = base_info_array_get_next(ba);
                base_mark_as_indel(bi);
                bi->aln_pos = c->pos + ref_offset;
                ba->count++;

            }
        } else if (op == BAM_CINS) {
            // Skip ahead (this avoids processing any actual inserted base, which would have no reference position)
            query_pos += len;

            // ANCHOR BASE PATCH (ONLY TO MATCH CURRENT OUTPUT!)
            // TODO: verify whether this case should be skipping any preceding deletion up to the anchor...?
            patch_anchor_base(ba, c->pos + ref_offset);

        } else {
            // TODO: consider behaviour with soft-clipped bases!
            query_pos  += static_cast<int32_t>(cigar_consumes_query[op]) * len;
            ref_offset += static_cast<int32_t>(cigar_consumes_ref[op])   * len;
        }
    }

    return 0;
}
