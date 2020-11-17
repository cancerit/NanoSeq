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

#ifndef PILEUP_H_
#define PILEUP_H_

#include <limits.h>
#include <array>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <iostream>
#include <string>
#include "gzstream.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "bed_reader.h"
#include "options.h"
#include "read_bundler.h"
#include "writeout.h"


typedef struct {
  htsFile* fp;
  hts_itr_t* iter;
  int min_mapQ;
  sam_hdr_t* head;
} aux_t;


class Pileup {
 private:
    Options *opts;
    faidx_t* fai;
    Bed mask;
    Bed snp;
    int n;
    int tid;
    aux_t **data;
    bam_mplp_t mplp;
    int *n_plp;
    const bam_pileup1_t **plp;
    ogzstream  gzout;


 public:
    void Initiate(Options *options);

    std::string Header();

    std::string PositionString(int pos);

    char* GetTrinucleotideContext(int pos);

    void MultiplePileup();
};

#endif  // PILEUP_H_
