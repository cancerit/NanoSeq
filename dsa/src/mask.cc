#include <format>
#include "utils.h"
#include "mask.h"

Mask::Mask(const uint8_t index) {
  this->index = index;
  this->flag = static_cast<uint8_t>(1) << this->index;
}

void Mask::Init(const char *bed_fp) {
    this->bed_fp = bed_fp;

    if (bed_fp == NULL || bed_fp[0] == '\0') {
      throw std::runtime_error("No BED file name provided!");
    }

    this->f = hts_open(bed_fp, "r");
    if (this->f == NULL) {
      throw std::runtime_error(
        std::format("Error: failed to open {}\n", bed_fp));
    }

    this->tbx = tbx_index_load(bed_fp);
    if (this->tbx == NULL) {
      throw std::runtime_error(
        std::format("Error: failed to open .tbi index of {}\n", bed_fp));
    }
}

inline void Mask::UpdateMaskAt(const int start, const int pos, std::vector<uint8_t> mask) {
  mask[pos - start] |= this->flag;
}

void Mask::UpdateMask(const int start, const int end, const int a, const int b, std::vector<uint8_t> mask) {
  if (a == b) {
    UpdateMaskAt(start, a, mask);
  } else {
    for (int i = a; i <= b; ++i) {
      UpdateMaskAt(start, a, mask);
    }
  }
}

int Mask::LoadRange(const char *contig, const int start, const int end, std::vector<uint8_t> mask) {
  char region[MAX_REGION_STR_LENGTH];
  get_region(contig, start, end, region);
  hts_itr_t *itr = tbx_itr_querys(this->tbx, region);
  if (itr) {
    kstring_t str;
    int32_t nfields;
    int32_t *fields;
    int a, b;
    while (tbx_itr_next(this->f, this->tbx, itr, &str) >= 0) {
      fields = ksplit(&str, 0, &nfields);
      a = std::stoi(&str.s[fields[1]]);
      b = std::stoi(&str.s[fields[2]]);
      this->UpdateMask(start, end, a, b, mask);
    }
  }
  tbx_itr_destroy(itr);
  return 0;
}
