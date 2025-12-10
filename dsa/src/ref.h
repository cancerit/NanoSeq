#ifndef REF_H_
#define REF_H_

#include "slice.h"
#include "htslib/faidx.h"

class Ref {
  private:
    Slice<char> seq;

  public:
    faidx_t *fai;
    void Init(const char *fai_fp);
    void Fetch(const char *contig, const range_t range);
    char *From(const int32_t pos);
    const std::string ToString();
    std::string_view GetTripletAround(const int32_t pos);
};

#endif
