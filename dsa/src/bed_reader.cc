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

/**
  Tabix BED files prepared using:
    sort -k 1,1 -k 2,2n -k 3,3n in.bed | bgzip -c > out.sorted.bed.gz
    tabix -pbed out.sorted.bed.gz

**/

#include "bed_reader.h"
#include "bedtk_lite.h"

#define KS_SEP_LINE  2

int Bed::Load(const char *bed_filename, ogzstream &gzout, const bool out2stdout) {
  /*
  Changed: loading the full BED file to support multiple regions; no need for the tabix index.

  TODO:
  - throw exceptions instead of using return codes?
  - consider logging to stderr instead of gzout (verify expectations)
  */

  if (!bed_filename || bed_filename[0] == '\0') {
    return 1;
  }
  // Load the intervals
  this->intervals = read_bed3(bed_filename);
  if (!intervals) {
    return 1;
  }
  // Index the intervals
  cr_index(this->intervals);
  if (!out2stdout) {
    gzout << "# ";
    gzout << intervals->n_r;
    gzout << " intervals added from ";
    gzout << bed_filename;
    gzout << std::endl;
  }
  return 0;
}

// Returns true if position intersects with BED position, otherwise false
bool Bed::Intersects(const char *contig, const int pos) {
  const int64_t n_b = cr_overlap(this->intervals, contig, pos, pos, &this->b, &this->m_b);
  return n_b != 0;
}
