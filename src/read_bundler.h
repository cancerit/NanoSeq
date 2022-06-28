/*########## LICENCE ##########
# Copyright (c) 2022 Genome Research Ltd
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


typedef struct {
  std::string id;
  int beg;
  int end;
  int strand;
  std::string fwd_bc;
  std::string rev_bc;
} identifier;


struct bundle {
  std::map<int, std::map<int, float>> dplx_depth;
  std::map<int, float> bulk_depth;
  std::map<int, std::map<char, int>> counts;
  std::map<int, std::vector<int>> asxs;
  std::map<int, std::vector<int>> clip;
  std::map<int, std::vector<int>> nmms;
  std::map<int, std::vector<int>> ppair;
  std::map<int, std::vector<std::pair<char, int>>> call;
  std::map<int, std::vector<double>> consensus;
  identifier idf;
  int bundle_type;
};


typedef std::map<std::string, bundle> bundles;

typedef std::vector<const bam_pileup1_t*> pileups;


class ReadBundler {
 public:
    int pos;

    int offset;

    int AuxTagIsPresent(bam1_t* b, const char* tag);

    char* AuxTagToChar(bam1_t* b, const char* tag);

    int AuxTagToInt(bam1_t* b, const char* tag);

    int ASMinusXS(bam1_t* b);

    int IsFivePrimeClipped(bam1_t* b, int readstrand);

    bool IsTemplate(int beg, int end);

    std::string ReadOrientation(int readstrand, int readnumber);

    int ReadStrand(bam1_t* b);

    int ReadNumber(bam1_t* b);

    std::pair<char, int> BaseAndQual(const bam_pileup1_t* p);

    bool BulkIsUsable(bam1_t *b);

    int IsProperPair(bam1_t *b);

    identifier DplxIdentifier(const bam_pileup1_t* p);

    void UpdateDplxBundle(identifier idf, bundle* bndl,
      const bam_pileup1_t* p);

    void UpdateBulkBundle(bundle* bndl, const bam_pileup1_t* p, int min_base_quality);

    void DplxConsensus(bundle* bndl);

    bundles DplxBundles(int pos, int offset, int min_dplx_depth, pileups plps);

    bundle BulkBundle(pileups plps, int min_base_quality);

};

#endif  // READ_BUNDLER_H_
