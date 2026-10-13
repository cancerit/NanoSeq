/*########## LICENCE ##########
# Copyright (c) 2022, 2025, 2026 Genome Research Ltd
#
# Authors: Luca Barbon <lb29@sanger.ac.uk>
#
# This file is part of NanoSeq.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
##########################*/

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <stdint.h>

// DNA alphabet length
#define ALPH_LEN 4

// TODO: verify!
#define BAM_NT_A 1
#define BAM_NT_C 2
#define BAM_NT_G 4
#define BAM_NT_T 8

#define ALLELE_A 0
#define ALLELE_C 1
#define ALLELE_G 2
#define ALLELE_T 3
#define ALLELE_DEL 4
#define ALLELE_COUNT 5

#define ALLELE_INVALID 255

// ASSUMPTION: the 4-bit base has been validated as canonical beforehand
static const uint8_t canonical_nt16_minus_one_to_allele[8] = {
  [BAM_NT_A - 1] = ALLELE_A,
  [BAM_NT_C - 1] = ALLELE_C,
  [BAM_NT_G - 1] = ALLELE_G,
  [BAM_NT_T - 1] = ALLELE_T
};

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
