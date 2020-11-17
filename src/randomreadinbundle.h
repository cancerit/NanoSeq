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

#ifndef RANDOMIZE_DUPLICATES_H_
#define RANDOMIZE_DUPLICATES_H_

#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cassert>
#include "header.h"
#include "htslib/sam.h"

class RandomlySelectOneReadPerBundle {
 private:


 public:
   char* infile;

   char* outfile;

   htsFile *in;

   htsFile *out;

   hts_idx_t *idx;

   bam_hdr_t *head;

   int min_reads;

   void LoadFiles();

   bool BamIsCorrectlyPreprocessed(bam_hdr_t *head);

   void UpdateHeader(int argc, char **argv);

   int ReadStrand(bam1_t* b);

   void WriteOut(bam1_t* b);

   void SelectReads();

};

#endif  // RANDOMIZE_DUPLICATES_H_
