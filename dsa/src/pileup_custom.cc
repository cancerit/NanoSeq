#include "pileup_custom.h"
#include <assert.h>

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

static inline void base_init_indel(base_t *b) {
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

    base_t *b = base_info_get_last(a);
    if (
        b != NULL &&              // not the first base
        b->aln_pos == pos - 1 &&  // immediately preceding position
        b->base != ALLELE_DEL     // neither an indel nor a [patched] anchor
    ) {
        // NOTE: the quality score is not being overridden
        // ASSUMPTION: the quality score of an indel is never evaluated
        base_init_indel(b);
    }
}

/*
static inline void cigar_print(const bam1_t *read) {
    uint32_t *cigar = bam_get_cigar(read);
    for (uint32_t i = 0; i < read->core.n_cigar; ++i) {
        fprintf(stderr, "%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));
    }
    printf("\n");
}
*/

int base_array_update(base_array_t *ba, bam1_t *read, const uint8_t min_qual) {
    const bam1_core_t *c = &read->core;
    ba->start = static_cast<int32_t>(c->pos);

    // Expand output array if necessary
    const int32_t read_length = c->l_qseq;
    if (base_info_array_reset(ba, (uint64_t)read_length)) {
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

                // Verify it is a canonical base
                if (count_bits(bam_nt) == 1 && qual[query_pos] >= min_qual) {

                    // A. Push canonical base
                    bi = base_info_array_get_next(ba);
                    // bi->read_pos = (int16_t)query_pos;
                    bi->aln_pos = c->pos + ref_offset;
                    bi->base = canonical_nt16_minus_one_to_allele[bam_nt - 1];
                    bi->qual = qual[query_pos];
                    ba->count++;

                }

                ref_offset++;
            }
        } else if (op == BAM_CDEL) {
            /*
            cigar_print(read);
            fprintf(stderr, "\t%d\tDELETION at %d/%d\n", ba->start, (int32_t)read->core.pos + ref_offset, query_pos);
            */
            // ANCHOR BASE PATCH (ONLY TO MATCH CURRENT OUTPUT!)
            patch_anchor_base(ba, c->pos + ref_offset);
            // ref_offset += len;

            // Consider that in samtools pileup skips are marked as deletions (bug)!

            // NOTE: matching the bam_plp_* behaviour by ignoring non-reference bases
            ref_offset_end = ref_offset + len;
            for (; ref_offset < ref_offset_end; ++ref_offset) {

                // B. Push deletion
                bi = base_info_array_get_next(ba);
                // bi->read_pos = (int16_t)query_pos;
                base_init_indel(bi);
                bi->aln_pos = c->pos + ref_offset;
                ba->count++;

            }
        } else if (op == BAM_CINS) {
       	    /*
            cigar_print(read);
            fprintf(stderr, "\t%d\tINSERTION at %d/%d\n", ba->start, (int32_t)read->core.pos + ref_offset, query_pos);
            */
            // Skip ahead (this avoids processing any actual inserted base, which would have no reference position)
            query_pos  += len;

            // ANCHOR BASE PATCH (ONLY TO MATCH CURRENT OUTPUT!)
            // TODO: verify whether this case should be skipping any preceding deletion up to the anchor...?
            patch_anchor_base(ba, c->pos + ref_offset);

        } else {
            // TODO: consider behaviour with soft-clipped bases!
            query_pos  += (int32_t)cigar_consumes_query[op] * len;
            ref_offset += (int32_t)cigar_consumes_ref[op]   * len;
        }
    }

    return 0;
}
