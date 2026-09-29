#ifndef POS_STATS_H_
#define POS_STATS_H_

typedef struct pos_final_stats_t {
    int asxs = 0;
    float nm = 0.0f;
    int proper_pairs = 0;
    // Duplex-only
    int clip = 0;
} pos_final_stats_t;

#endif
