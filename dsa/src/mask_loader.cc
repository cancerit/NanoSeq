#include <format>
#include "utils.h"
#include "mask_loader.h"

MaskLoader::MaskLoader(const uint8_t index) {
  this->index = index;
  this->flag = static_cast<uint8_t>(1) << this->index;
}

void MaskLoader::Init(const char *bed_fp) {
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

void MaskLoader::LoadMask(const char *contig, const int start, const int end, Mask &mask) {
  mask.Reset({start, end});

  char region[MAX_REGION_STR_LENGTH];
  get_region(contig, start, end, region);
  hts_itr_t *itr = tbx_itr_querys(this->tbx, region);
  range_t r;
  if (itr) {
    kstring_t str;
    int32_t nfields;
    int32_t *fields;
    while (tbx_itr_next(this->f, this->tbx, itr, &str) >= 0) {
      fields = ksplit(&str, 0, &nfields);
      r.start = std::stoi(&str.s[fields[1]]);
      r.end = std::stoi(&str.s[fields[2]]);
      mask.Update(r, this->flag);
    }
  }
  tbx_itr_destroy(itr);
}
