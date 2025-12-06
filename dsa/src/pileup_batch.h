#ifndef PILEUP_BATCH_H_
#define PILEUP_BATCH_H_

#include <vector>
#include <assert.h>
#include "constants.h"

class PileupBatch {
  private:
    // const char region[MAX_REGION_STR_LENGTH];
    int32_t start;
    int32_t end;

  public:
    void Update(const int start, const int end);
    std::vector<uint8_t> mask;
};

#endif
