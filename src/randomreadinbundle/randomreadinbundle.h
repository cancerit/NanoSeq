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
