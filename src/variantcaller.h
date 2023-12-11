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

#ifndef VARIANTCALLER_H_
#define VARIANTCALLE_H_

#include <algorithm>
#include <map>
#include <sstream>
#include <cassert>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <string.h>
#include <tuple>
#include <vector>
#include "limits.h"
#include "float.h"
#include <unistd.h>
#include "gzstream.h"

typedef struct {
  std::string chrom;
  int chrom_beg;
  std::string context;
  int snp;
  int shearwater;
  float bulk_asxs;
  float bulk_nm;
  int bfwd_A;
  int bfwd_C;
  int bfwd_G;
  int bfwd_T;
  int bfwd_I;
  int brev_A;
  int brev_C;
  int brev_G;
  int brev_T;
  int brev_I;
  int bp_beg;
  int bp_end;
  int bndl_type;
  float dplx_asxs;
  float dplx_clip;
  float dplx_nm;
  int f1r2_A;
  int f1r2_C;
  int f1r2_G;
  int f1r2_T;
  int f1r2_I;
  int f2r1_A;
  int f2r1_C;
  int f2r1_G;
  int f2r1_T;
  int f2r1_I;
  int f1r2_A_Q;
  int f1r2_C_Q;
  int f1r2_G_Q;
  int f1r2_T_Q;
  int f2r1_A_Q;
  int f2r1_C_Q;
  int f2r1_G_Q;
  int f2r1_T_Q;
  float bfwd_canonical;
  float brev_canonical;
  float f1r2_canonical;
  float f2r1_canonical;
  float bfwd_total;
  float brev_total;
  float bulk_ppair;
  float dplx_ppair;
  float f1r2_total;
  float f2r1_total;
  int left;
  int right;
  int min_qpos;
  char f1r2_call;
  char f2r1_call;
  char call;
  int isvariant;
  int ismasked;
  std::string pyrcontext;
  std::string dplx_barcode;

  // Filters
  bool dplx_clip_filter;
  bool alignment_score_filter;
  bool mismatch_filter;
  bool matched_normal_filter;
  bool duplex_filter;
  bool consensus_base_quality_filter;
  bool indel_filter;
  bool five_prime_trim_filter;
  bool three_prime_trim_filter;
  bool proper_pair_filter;
  bool vaf_filter;
  bool pass_all_filters;
} row_t;


class VariantCaller {
  public:
    igzstream gzin;
      
    std::ofstream fout;

    std::ofstream fout_coverage;

    std::ofstream fout_discarded;

    ogzstream gzout_coverage;

    const char* outfile;

    const char* outfile_coverage;

    const char* outfile_discarded;

    const char* bed;

    int asxs;

    int bulk;
    
    int bulk_total; // added by fa8

    int dplx;

    float clip;

    float frac;

    float indel;

    int nmms;

    float ppair;

    int qual;

    int readlen;

    int min_cycle;

    int max_cycle;

    float vaf;

    int coverage;

    std::map<char, std::map<int, std::map<int, int>>> call_by_qpos;

    std::map<std::string, std::map<int, int>> pyr_by_mask;

    std::map<int, std::map<int, std::map<int, std::map<int, int>>>> read_bundles;

    std::map<int, std::map<int, int>> burdens;

    row_t ParseRow(std::string line);

    void CallDuplex(row_t *row);

    int DplxClipFilter(row_t *row);

    int AlignmentScoreFilter(row_t *row);

    int MismatchFilter(row_t *row);

    int MatchedNormalFilter(row_t *row);

    int DuplexFilter(row_t *row);

    int ConsensusBaseQualityFilter(row_t *row);

    int IndelFilter(row_t *row);

    int FivePrimeTrimFilter(row_t *row);

    int ThreePrimeTrimFilter(row_t *row);

    int ProperPairFilter(row_t *row);

    void ApplyFilters(row_t *row);
    
    int VafFilter(row_t *row);

    int IsVariant(row_t *row);

    int IsMasked(row_t *row);

    std::string ReverseComplementContext(std::string seq);

    std::string PyrimidineContext(row_t *row);

    std::string PyrimidineSubstitution(row_t *row);

    std::string StrandContext(row_t *row);

    std::string StrandSubstitution(row_t *row);

    std::string StrandMismatch(row_t *row);

    bool ContextIsCanonical(row_t *row);

    void CollectMetrics();

    void WriteVariants(row_t *row);

    void WriteDiscardedVariants(row_t *row);

    void WriteMismatches(row_t *row);

    void WriteMetrics();

};

#endif  // VARIANTCALLER_H_
