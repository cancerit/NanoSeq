#ifndef PROBS_H_
#define PROBS_H_

#include <cstdint>
#include <cmath>

#include "constants.h"

#define ALT_BASES 3.0  // ALPH_LEN - 1
#define POWER 10.0

#define PROB_INDEX_ERR 0
#define PROB_INDEX_CORRECT 1
#define PROB_STRIDE 2

#define PROB_ARRAY_LENGTH (256 * PROB_STRIDE)

typedef struct probs_t {
  // pairs of probabilities (error and correct) for each quality score
  double probs[PROB_ARRAY_LENGTH];
} probs_t;

inline void probs_init(probs_t *p) {
  p->probs[PROB_INDEX_ERR] = 0.0;
  p->probs[PROB_INDEX_CORRECT] = 0.0;

  double p_error;
  for (int i = 0; i <= UINT8_MAX; ++i) {
    p_error = std::pow(POWER, (-i / POWER));
    p->probs[i * PROB_STRIDE + PROB_INDEX_ERR]     = std::log10(p_error);
    p->probs[i * PROB_STRIDE + PROB_INDEX_CORRECT] = std::log10((1.0 - p_error) / ALT_BASES);
  }
}

inline void probs_add_p_error(const probs_t *p, const int qual, int base_code, double *probs) {
    if (base_code != ALLELE_DEL) {
        const double pc = p->probs[qual * PROB_STRIDE + PROB_INDEX_CORRECT];
        double dp[ALLELE_COUNT - 1] = {pc, pc, pc, pc};
        dp[base_code] = p->probs[qual * PROB_STRIDE + PROB_INDEX_ERR];
        probs[ALLELE_A] += dp[ALLELE_A];
        probs[ALLELE_C] += dp[ALLELE_C];
        probs[ALLELE_G] += dp[ALLELE_G];
        probs[ALLELE_T] += dp[ALLELE_T];
    }
}

inline void finalise_consensus_quality_scores(double probs[ALPH_LEN]) {
    // log10 sum exp
    double maxp = probs[0];
    for (size_t i = 1; i < ALPH_LEN; ++i) {
        if (probs[i] > maxp) {
            maxp = probs[i];
        }
    }

    double sumexp = 0.0;
    for (size_t i = 0; i < ALPH_LEN; ++i) {
        sumexp += std::pow(POWER, probs[i] - maxp);
    }

    const double logsumexp = std::log10(sumexp) + maxp;
    // normalize sum log10 probs and convert to Phred based quality score
    for (size_t i = 0; i < ALPH_LEN; ++i) {
        probs[i] = (probs[i] - logsumexp) * -POWER;
    }
}

#endif
