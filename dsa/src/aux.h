#ifndef AUX_H_
#define AUX_H_

typedef struct {

    // Set at startup
    htsFile *fp;
    sam_hdr_t *head;
    hts_idx_t *idx;
    int min_mapQ;
    int duplex;

    // Set per genomic range
    range_tid_t range;
    hts_itr_t *iter;
    uint64_t iterations;

} aux_t;

static void aux_init(aux_t *a, const char *fp) {
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

static void aux_duplex_init(aux_t *a, const char *fp, const int min_map_q) {
    aux_init(a, fp);
    a->min_mapQ = min_map_q;
    a->duplex = 1;
}

static void aux_bulk_init(aux_t *a, const char *fp) {
    aux_init(a, fp);
    a->min_mapQ = 0;
    a->duplex = 0;
}

static void aux_reset(aux_t *a) {
    // Iterator
    if (a->iter != NULL) {
        sam_itr_destroy(a->iter);
        a->iter = NULL;
    }

    // Stats
    a->iterations = 0;

    // Range
    a->range.tid = -1;
    a->range.start = 0;
    a->range.end = 0;
}

static void aux_set_iterator(aux_t *a, const range_tid_t range) {
    aux_reset(a);
    a->range = range;
    a->iter = sam_itr_queryi(a->idx, range.tid, range.start, range.end + 1);
    if (a->iter == NULL) {
        fprintf(stderr, "Failed to initialise iterator!\n");
    }
}

static int aux_iter(aux_t *a, bam1_t *b) {
    a->iterations++;
    return sam_itr_next(a->fp, a->iter, b);
}

#endif  // AUX_H_
