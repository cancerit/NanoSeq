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

/**
  Tabix BED files prepared using:
    sort -k 1,1 -k 2,2n -k 3,3n in.bed | bgzip -c > out.sorted.bed.gz
    tabix -pbed out.sorted.bed.gz

**/

#include "bed_reader.h"

#define KS_SEP_LINE  2


// Load intervals from tabix BED file
void Bed::Load(const char* bed_filename, const char* rname, int beg, int end, ogzstream & gzout) {
  htsFile *fp = hts_open(bed_filename, "r");
  if (fp == NULL) {
    std::stringstream er;
    er << "Error: failed to open ";
    er << bed_filename;
    er << std::endl;
    throw std::runtime_error(er.str());
  }
  tbx_t* tbx = tbx_index_load(bed_filename);
  if (!tbx) {
    std::stringstream er;
    er << "Error: failed to open .tbi index of ";
    er << bed_filename;
    er << std::endl;
    throw std::runtime_error(er.str());
  }
  char region[65536];
  sprintf(region, "%s:%d-%d", rname, beg, end);
  hts_itr_t* itr = tbx_itr_querys(tbx, region);
  if (itr) {
    kstring_t str;
    int32_t nfields;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
      int32_t* fields = ksplit(&str, 0, &nfields);
      int beg = std::stoi(&str.s[fields[1]]);
      int end = std::stoi(&str.s[fields[2]]);
      intervals.push_back(std::make_pair(beg, end));
    }
  }
  tbx_itr_destroy(itr);
  gzout << "# ";
  gzout << intervals.size();
  gzout << " intervals added from ";
  gzout << bed_filename;
  gzout << " ";
  gzout << region;
  gzout << std::endl;
}


// Returns true if position intersects with BED position, otherwise false
bool Bed::Intersects(int pos) {
  bool ret = false;
  while (intervals.size() > 0) {
    std::pair<int, int> interval = intervals[0];
    int beg = interval.first;
    int end = interval.second;
    if (pos < beg) {
      ret = false;
      break;
    } else if (pos >= end) {
      intervals.pop_front();
    } else {
      assert((pos >= beg) && (pos < end));
      ret = true;
      break;
    }
  }
  return ret;
}
