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

#ifndef BED_H_
#define BED_H_

#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "htslib/tbx.h"
#include "htslib/kstring.h"
#include "gzstream.h"

class Bed {
 public:
   const char* rname;

   std::deque<std::pair<int, int>> intervals;

   void Load(const char* bed_filename, const char* rname, int beg, int end, ogzstream & gzout, bool out2stdout);

   bool Intersects(int pos);

};

#endif  // BED_H_
