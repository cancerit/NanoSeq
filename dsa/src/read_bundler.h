/*########## LICENCE ##########
# Copyright (c) 2022, 2025 Genome Research Ltd
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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

#ifndef READ_BUNDLER_H_
#define READ_BUNDLER_H_

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

#include "limits.h"

#include "htslib/sam.h"
#include "options.h"

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

typedef struct {
  std::string id;
  int beg;
  int end;
  int strand;
  std::string fwd_bc;
  std::string rev_bc;
} identifier;

struct bundle {
  float dplx_depth[STRAND_COUNT][READ_TYPE_COUNT];
  uint64_t counts[RTYPE_COUNT][ALLELE_COUNT];
  std::vector<int> asxs[RTYPE_COUNT];
  std::vector<int> clip[RTYPE_COUNT];
  std::vector<int> nmms[RTYPE_COUNT];

  uint64_t rtype_ppair_counts[RTYPE_COUNT];
  uint64_t rtype_read_counts[RTYPE_COUNT];  // then divide ppair to get the averages

  // TODO: replace character key with index (?)
  std::vector<std::pair<char, int>> call[RTYPE_COUNT];
  std::vector<double> consensus[BUNDLE_TYPES_COUNT];
  identifier idf;
  int bundle_type;
};

// TODO: consider a more compact duplex ID as key
typedef std::map<std::string, bundle> bundles;
typedef std::vector<const bam_pileup1_t*> pileups;

static inline int read_has_flag(const bam1_t *b, const uint16_t flag) {
  return (b->core.flag & flag) != 0;
}

static inline int read_is_in_proper_pair(const bam1_t *b) {
  return read_has_flag(b, BAM_FPROPER_PAIR);
}

class ReadBundler {
  public:
    int pos;
    int offset;
    char* AuxTagToChar(bam1_t* b, const char* tag);
    int AuxTagToInt(bam1_t* b, const char* tag);
    int ASMinusXS(bam1_t* b);
    int IsFivePrimeClipped(bam1_t* b, int readstrand);
    bool IsTemplate(const int beg, const int end);
    std::pair<int, int> BaseAndQual(const bam_pileup1_t* p);
    bool BulkIsUsable(bam1_t *b);
    identifier DplxIdentifier(const bam_pileup1_t* p);
    void UpdateDplxBundle(identifier idf, bundle* bndl, const bam_pileup1_t* p);
    void UpdateBulkBundle(bundle* bndl, const bam_pileup1_t* p, int min_base_quality);
    void DplxConsensus(bundle* bndl);
    bundles DplxBundles(int pos, int offset, int min_dplx_depth, pileups plps);
    bundle BulkBundle(pileups plps, int min_base_quality);
};

#endif  // READ_BUNDLER_H_
