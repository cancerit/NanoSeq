#ifndef PILEUP_STATE_H_
#define PILEUP_STATE_H_

#include "aux.h"
#include "compressor.h"
#include "options.h"
#include "ref.h"

typedef struct pileup_state_t {
    aux_t *bulk_aux;
    aux_t *duplex_aux;
    Ref *ref;
    Options *opts;
    GzipCompressor *compressor;
} pileup_state_t;

#endif
