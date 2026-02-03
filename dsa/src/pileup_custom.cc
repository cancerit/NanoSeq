#include "pileup_custom.h"

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

static inline void patch_anchor_base(const base_array_t *a) {
    /*
    Mark the anchor base of an indel as a deletion, matching the original behaviour:

    if ((p->is_del) || (p->indel != 0)) {
        return std::make_pair('i', -1);
    }
    */
    base_t *b = base_info_get_last(a);
    if (b != NULL && b->base != ALLELE_DEL) {
        b->base = ALLELE_DEL;
        b->qual = 0;
    }
}

int base_array_update(base_array_t *ba, bam1_t *read, const uint8_t min_qual) {
    const bam1_core_t *c = &read->core;

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
                    bi->aln_pos = read->core.pos + ref_offset;
                    bi->base = nt16_allele[bam_nt];
                    bi->qual = qual[query_pos];
                    ba->count++;

                }

                ref_offset++;
            }
        } else if (op == BAM_CDEL) {

            // ANCHOR BASE PATCH (ONLY TO MATCH CURRENT OUTPUT!)
            patch_anchor_base(ba);

            // Consider that in samtools pileup skips are marked as deletions (bug)!
            ref_offset_end = ref_offset + len;
            for (; ref_offset < ref_offset_end; ++ref_offset) {

                // B. Push deletion
                bi = base_info_array_get_next(ba);
                // bi->read_pos = (int16_t)query_pos;
                bi->aln_pos = read->core.pos + ref_offset;
                bi->base = ALLELE_DEL;
                bi->qual = 0;
                ba->count++;

            }
        } else if (op == BAM_CINS) {

            // ANCHOR BASE PATCH (ONLY TO MATCH CURRENT OUTPUT!)
            patch_anchor_base(ba);

        } else {
            // TODO: consider behaviour with soft-clipped bases!
            query_pos  += (int32_t)cigar_consumes_query[op] * len;
            ref_offset += (int32_t)cigar_consumes_ref[op]   * len;
        }
    }

    return 0;
}
