#ifndef READ_INFO_H_
#define READ_INFO_H_

#include "htslib/sam.h"
#include <format>
#include <string>
#include "constants.h"

typedef struct read_info_t {
    int32_t strand;
    int32_t read_index;

    int64_t asxs;
    int64_t nm;
    int64_t proper_pair;

    // Duplex-only
    int64_t is_5p_clipped;
} read_info_t;

inline int32_t read_info_get_r_type(const read_info_t *r) {
    return r->strand != STRAND_INDEX_IGNORE ? RTYPES[r->strand][r->read_index] : 0;
}

inline int read_has_flag(const bam1_t *b, const uint16_t flag) {
    return (b->core.flag & flag) != 0;
}

inline bool read_has_tag(const bam1_t *read, const char *tag) {
    return bam_aux_get(read, tag) != NULL;
}

inline int get_strand_index(const bam1_t *b) {
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

inline const uint8_t *get_tag(const bam1_t *read, const char *tag) {
    const uint8_t *t = bam_aux_get(read, tag);
    if (t == NULL) {
        throw std::runtime_error(std::format(
            "Missing {} tag!", tag));
    }
    return t;
}

inline std::string get_duplex_id(const bam1_t *read) {
    return bam_aux2Z(get_tag(read, "RB"));
}

inline int64_t get_int_tag(const bam1_t *read, const char *tag) {
    return bam_aux2i(get_tag(read, tag));
}

inline int64_t get_as_minus_xs(const bam1_t *read) {
    const int64_t as = get_int_tag(read, "AS");
    const int64_t xs = get_int_tag(read, "XS");
    return as - xs;
}

inline int64_t read_is_in_proper_pair(const bam1_t *b) {
  return static_cast<int64_t>(read_has_flag(b, BAM_FPROPER_PAIR));
}

inline int get_read_type_index(const bam1_t *b) {
  // ASSUMPTION: that one and only one BAM_FREAD* flag is set to be verified earlier
  //  (this would be an unpredictable branch, while flag consistency verification is predictable)
  static_assert(READ_TYPE_INDEX_READ_1 == 0);
  static_assert(READ_TYPE_INDEX_READ_2 == 1);
  return read_has_flag(b, BAM_FREAD2);
}

inline int64_t get_nm(const bam1_t *read) {
    return get_int_tag(read, "NM");
}

inline int get_is_5p_clipped(const bam1_t *b, const int32_t strand) {
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

inline int get_is_proper_pair(const bam1_t *b) {
    return read_has_flag(b, BAM_FPROPER_PAIR);
}

inline void read_info_init(read_info_t *r, const bam1_t *read) {
    r->strand = get_strand_index(read);
    r->read_index = get_read_type_index(read);

    r->asxs = get_as_minus_xs(read);
    r->nm = get_nm(read);
    r->proper_pair = read_is_in_proper_pair(read);
}

#endif
