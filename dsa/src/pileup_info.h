#ifndef PILEUP_INFO_H_
#define PILEUP_INFO_H_

#include <htslib/sam.h>
#include "aux.h"

typedef struct pileup_t {
    bam_plp_t plp;

    int32_t count;
    bam_pileup1_t *plps;

} pileup_t;

static void pileup_init(pileup_t *p,  bam_plp_auto_f f, aux_t *aux, const int32_t max_count) {
    p->count = 0;
    p->plps = NULL;

    p->plp = bam_plp_init(f, aux);
    bam_plp_set_maxcnt(p->plp, max_count);
}

#endif
