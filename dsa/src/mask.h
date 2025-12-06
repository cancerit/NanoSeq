#ifndef MASK_H_
#define MASK_H_

#include "htslib/tbx.h"
#include "htslib/kstring.h"
#include <vector>
#include "range.h"

class Mask {
  private:
    range_t range;
    std::vector<uint8_t> mask;

  public:
    void Reset(const range_t range);
    void Update(const range_t range, const uint8_t flag);
    uint8_t GetFlag(const int32_t pos);
};

#endif
