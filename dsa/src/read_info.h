#ifndef READ_INFO_H_
#define READ_INFO_H_

#include "htslib/sam.h"
#include <string.h>
#include <iostream>

#include <cassert>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <utility>

#include "duplex_tag_info.h"
#include "constants.h"

static inline int read_has_flag(const bam1_t *b, const uint16_t flag) {
    return (b->core.flag & flag) != 0;
}

static inline int get_strand_index(const bam1_t *b) {
  // ASSUMPTION: proper pair and strand have already been verified
  static_assert(STRAND_INDEX_FORWARD == 0);
  static_assert(STRAND_INDEX_REVERSE == 1);

  // TODO: optimise
  if (read_has_flag(b, BAM_FMREVERSE)) {
    return STRAND_INDEX_FORWARD;
  } else if (read_has_flag(b, BAM_FREVERSE)) {
    return STRAND_INDEX_REVERSE;
  } else {
    return STRAND_INDEX_IGNORE;
  }
}

static inline const uint8_t *get_tag(const bam1_t *read, const char *tag) {
    const uint8_t *t = bam_aux_get(read, tag);
    if (t == NULL) {
        throw std::runtime_error(std::format(
            "Missing {} tag!", tag));
    }
    return t;
}

static inline std::string get_duplex_id(const bam1_t *read) {
    return bam_aux2Z(get_tag(read, "RB"));
}

static inline int64_t get_int_tag(const bam1_t *read, const char *tag) {
    return bam_aux2i(get_tag(read, tag));
}

static inline int64_t get_as_minus_xs(const bam1_t *read) {
    const int64_t as = get_int_tag(read, "AS");
    const int64_t xs = get_int_tag(read, "XS");
    return as - xs;
}

static inline int read_is_in_proper_pair(const bam1_t *b) {
  return read_has_flag(b, BAM_FPROPER_PAIR);
}

static inline int get_read_type_index(const bam1_t *b) {
  // ASSUMPTION: that one and only one BAM_FREAD* flag is set to be verified earlier
  //  (this would be an unpredictable branch, while flag consistency verification is predictable)
  static_assert(READ_TYPE_INDEX_READ_1 == 0);
  static_assert(READ_TYPE_INDEX_READ_2 == 1);
  return read_has_flag(b, BAM_FREAD2);
}

static inline int64_t get_nm(const bam1_t *read) {
    return get_int_tag(read, "NM");
}

static inline int get_is_5p_clipped(const bam1_t *b, const int32_t strand) {
    const uint32_t *cigar = bam_get_cigar(b);
    // TODO: optimise
    switch (strand) {
        case STRAND_INDEX_FORWARD:
            return (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP);
        case STRAND_INDEX_REVERSE:
            return (bam_cigar_op(cigar[b->core.n_cigar - 1]) == BAM_CSOFT_CLIP);
        default:
            return 0;
    }
}

static inline int get_is_proper_pair(const bam1_t *b) {
    return read_has_flag(b, BAM_FPROPER_PAIR);
}

typedef struct ReadInfo {
    int32_t strand;
    int64_t as_xs;
    int64_t nm;
    uint8_t is_clipped;
    uint8_t is_proper_pair;
    std::string id;
    // duplex_tag_info idf;
} ReadInfo;

static void read_info_init(ReadInfo *ri, const bam1_t *read) {
    ri->id = get_duplex_id(read);
    ri->strand = get_strand_index(read);
    ri->as_xs = get_as_minus_xs(read);
    ri->nm = get_nm(read);
    ri->is_clipped = get_is_5p_clipped(read, ri->strand);
    ri->is_proper_pair = get_is_proper_pair(read);
}

#endif
