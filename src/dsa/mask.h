#ifndef MASK_H_
#define MASK_H_

#include "slice.h"

class Mask {
  private:
    Slice<uint8_t> mask;

  public:
    uint64_t CountBytesSet();
    void Reset(const range_t range);
    void Update(const range_t range, const uint8_t flag);
    uint8_t GetFlag(const int32_t pos);
};

#endif
