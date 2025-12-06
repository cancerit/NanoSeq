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
#include "mask.h"
#include "options.h"
#include "read_bundler.h"
#include "writeout.h"
#include "constants.h"

typedef struct {
  htsFile* fp;
  hts_itr_t* iter;
  int min_mapQ;
  int duplex;
  sam_hdr_t* head;
} aux_t;

typedef struct {
  int32_t tid;
  int32_t start;
  int32_t end;
} range_t;

class Pileup {
  private:
    Options *opts;
    faidx_t *fai;

    // BAI/CRAI indices for sample and normal
    hts_idx_t *indices[BAM_COUNT];

    const char *regions;  // Regions to process
    Mask masks[MASK_COUNT];

    aux_t **data;
    bam_mplp_t mplp;
    ogzstream gzout;

  public:
    Pileup();
    void Initiate(Options *options);
    void InitIterators(const range_t *r);
    std::string Header();
    std::string PositionString(const char *contig, const int pos);
    void MultiplePileup();
};

/*
int xy(const char *fp) {
  gzFile f = gzopen(fp, "r");
  if (f == NULL) {
      fprintf(stderr, "\nFailed to open file '%s'!\n", fp);
      return 1;
  }

  kstring_t str;
  uint64_t total = 0;
  kstream_t *ks = ks_init(f);
  char *contig, *rest;
  int32_t start, end;
  while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
    total++;

    contig = parse_bed3b(str.s, &start, &end, &rest);
    if (contig == NULL) {
        fprintf(stderr, "\nContig not found!\n");
        return 1;
    }
  }
  return 0;
}
*/

#endif  // PILEUP_H_
