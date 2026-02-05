#ifndef POS_STATS_H_
#define POS_STATS_H_

typedef struct pos_final_stats_t {
    double asxs = 0.0;
    double nm = 0.0;
    double proper_pairs = 0.0;
    // Duplex-only
    double clip = 0.0;
} pos_final_stats_t;

#endif
