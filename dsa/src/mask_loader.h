#ifndef MASK_LOADER_H_
#define MASK_LOADER_H_

#include "htslib/tbx.h"
#include "htslib/kstring.h"
#include "mask.h"

class MaskLoader {
  private:
    uint8_t index;
    uint8_t flag;
    const char *bed_fp;
    htsFile *f;
    tbx_t *tbx;

  public:
    MaskLoader() = default;
    MaskLoader(const uint8_t index);
    void Init(const char *bed_fp);
    void LoadMask(const char *contig, const int start, const int end, Mask &mask);
};

#endif
