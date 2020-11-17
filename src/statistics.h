/**   LICENCE
* Copyright (c) 2020 Genome Research Ltd.
* 
* Author: Cancer Genome Project <cgphelp@sanger.ac.uk>
* 
* This file is part of NanoSeq.
* 
* NanoSeq is free software: you can redistribute it and/or modify it under
* the terms of the GNU Affero General Public License as published by the Free
* Software Foundation; either version 3 of the License, or (at your option) any
* later version.
* 
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
* details.
* 
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef STATISTICS_H_
#define STATISTICS_H_
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <limits.h>
#include <map>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include "htslib/sam.h"

class Statistics {
 private:


 public:
   int min_dplx_depth;
 
   char* input_bam;

   char* filtered_bam;
 
   std::string dir;

   htsFile *in1;

   htsFile *in2;

   hts_idx_t *idx1;

   hts_idx_t *idx2;

   bam_hdr_t *head1;

   bam_hdr_t *head2;

   int reads;

   std::map<int, int> reads_per_bundle;

   std::map<int, int> reads_per_usable_bundle;

   std::map<float, int> tlens_per_usable_bundle;

   std::ofstream filter_metrics;

   std::ofstream bundle_metrics;

   std::ofstream tlen_metrics;

   void LoadFiles();

   int ReadStrand(bam1_t* b);

   int ReadNumber(bam1_t* b);

   std::vector<std::string> Tokenize(std::string str, char delimiter);

   bool BamIsPreProcessed(bam_hdr_t *head);

   float MeanOfOneVector(std::vector<int> v1);

   void CollectMetrics();

   void WriteOut();

};

#endif  // STATISTICS_H_
