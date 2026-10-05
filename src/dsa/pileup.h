/*########## LICENCE ##########
# Copyright (c) 2022, 2026 Genome Research Ltd
#
# Author: Luca Barbon <lb29@sanger.ac.uk>
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


#ifndef PILEUP_H_
#define PILEUP_H_

#include <limits.h>
#include <string>
#include <vector>
#include <string>
#include "htslib/sam.h"
#include "mask_loader.h"
#include "options.h"
#include "constants.h"
#include "range.h"
#include "ref.h"
#include "aux.h"

class Pileup {
  private:
    Options *opts;
    Ref ref;

    // BAI/CRAI indices for sample and normal
    hts_idx_t *indices[BAM_COUNT];

    const char *regions;  // Regions to process
    MaskLoader masks[MASK_COUNT];

    aux_t data[BUNDLE_TYPES_COUNT];
    std::vector<genomic_region_t> ranges;
    int GetTID(const char *contig);
    const char *GetContig(const int32_t tid);
    void LoadRanges();

    sam_hdr_t *GetHeader(const int i);
    sam_hdr_t *GetBulkHeader();
    sam_hdr_t *GetDuplexHeader();

  public:
    Pileup();
    void DestroyIterators();
    void Initiate(Options *options);
    void InitIterators(const genomic_region_t *r);
    std::string Header();
    void MultiplePileup();
};

#endif  // PILEUP_H_
