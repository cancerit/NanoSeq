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

#ifndef WRITEOUT_H_
#define WRITEOUT_H_

#include <algorithm>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <numeric>
#include "gzstream.h"
#include "read_bundler.h"


class WriteOut {
 public:
    Options *opts;

    int Round(double d);

    std::string BaseQuals(bundle *bin);

    std::string Counts(bundle *bin);

    float MeanOfOneVector(std::vector<int> v1);

    float MeanOfTwoVectors(std::vector<int> v1, std::vector<int> v2);

    std::string BulkCovariates(bundle *bin);

    std::string DplxCovariates(bundle *bin);

    std::string ProperPairs(bundle *bin);

    std::string IdentifierString(bundle *bin);

    void WriteRows(bundle bulk, bundles dplx, std::string posn, ogzstream & gzout);
};

#endif  // WRITEOUT_H_
