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
