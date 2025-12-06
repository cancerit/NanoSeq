#ifndef MASK_H_
#define MASK_H_

#include "htslib/tbx.h"
#include "htslib/kstring.h"

class Mask {
  private:
    uint8_t index;
    uint8_t flag;
    const char *bed_fp;
    htsFile *f;
    tbx_t *tbx;
    void UpdateMaskAt(const int start, const int pos, std::vector<uint8_t> mask);
    void UpdateMask(const int start, const int end, const int a, const int b, std::vector<uint8_t> mask);

  public:
    Mask() = default;
    Mask(const uint8_t index);
    void Init(const char *bed_fp);
    int LoadRange(const char *contig, const int start, const int end, std::vector<uint8_t> mask);
};

#endif
