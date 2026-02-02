#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// DNA alphabet length
#define ALPH_LEN 4

// TODO: verify!
#define BAM_NT_A 1
#define BAM_NT_C 2
#define BAM_NT_G 4
#define BAM_NT_T 8

#define ALLELE_DISCARDED 0
#define ALLELE_A 1
#define ALLELE_C 2
#define ALLELE_G 3
#define ALLELE_T 4
#define ALLELE_DEL 5
#define ALLELE_COUNT 6

static const int nt16_allele[16] = {
  [0] = ALLELE_DISCARDED,
  [BAM_NT_A] = ALLELE_A,
  [BAM_NT_C] = ALLELE_C,
  [3] = ALLELE_DISCARDED,
  [BAM_NT_G] = ALLELE_G,
  [5] = ALLELE_DISCARDED,
  [6] = ALLELE_DISCARDED,
  [7] = ALLELE_DISCARDED,
  [BAM_NT_T] = ALLELE_T,
  [9] = ALLELE_DISCARDED,
  [10] = ALLELE_DISCARDED,
  [11] = ALLELE_DISCARDED,
  [12] = ALLELE_DISCARDED,
  [13] = ALLELE_DISCARDED,
  [14] = ALLELE_DISCARDED,
  [15] = ALLELE_DISCARDED,
};

#define ALLELE_INDEX_A 0
#define ALLELE_INDEX_C 1
#define ALLELE_INDEX_G 2
#define ALLELE_INDEX_T 3

#define STRAND_COUNT 2
#define READ_TYPE_COUNT 2

static const int RTYPES[STRAND_COUNT][READ_TYPE_COUNT] = {
  {0, 1},
  {1, 0}
};

#define BAM_COUNT 2
#define MASK_COUNT 2

#define BULK_INDEX 0
#define DUPLEX_INDEX 1

// Bundle type indices
#define BUNDLE_TYPE_BULK   BULK_INDEX
#define BUNDLE_TYPE_DUPLEX DUPLEX_INDEX
#define BUNDLE_TYPES_COUNT 2

// Mask indices
#define MASK_INDEX_SNP   0
#define MASK_INDEX_NOISE 1
#define MASK_COUNT 2

#define MASK_FLAG_SNP   (1 << MASK_INDEX_SNP)
#define MASK_FLAG_NOISE (1 << MASK_INDEX_NOISE)

// For genomic region strings (contig:start-end)
#define MAX_REGION_STR_LENGTH 1024

#define READ_TYPE_INDEX_READ_1 0
#define READ_TYPE_INDEX_READ_2 1

#define STRAND_INDEX_FORWARD 0
#define STRAND_INDEX_REVERSE 1
#define STRAND_INDEX_IGNORE -1

// RTYPE is just the strand index for bulk
#define RTYPE_A 0
#define RTYPE_B 1
#define RTYPE_COUNT 2

#define MAX_DSA_LINE_LENGHT 4096

#define BED_INDEX_CONTIG 0
#define BED_INDEX_START 1
#define BED_INDEX_END 2

#endif
