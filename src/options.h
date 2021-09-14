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

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <array>

struct Options {
  std::array<const char*, 2> bams;
  std::array<const char*, 2> beds;
  const char* fasta;
  const char* rname;
  char* oname;
  int min_base_quality;
  int min_mapQ;
  int max_plp_depth;
  int min_dplx_depth;
  int max_dplx_depth;
  int offset;
  int beg;
  int end;
  bool out2stdout;
  bool doTests;
};

#endif  // OPTIONS_H_
