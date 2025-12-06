#ifndef MASK_H_
#define MASK_H_

#include "htslib/tbx.h"
#include "htslib/kstring.h"
#include <vector>

class Mask {
  private:
    int32_t start, end;
    std::vector<uint8_t> mask;

  public:
    void Reset(const int32_t start, const int32_t end);
    void Update(const int32_t start, const int32_t end, const uint8_t flag);
};

#endif
