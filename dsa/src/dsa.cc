/*########## LICENCE ##########
# Copyright (c) 2022, 2025, 2026 Genome Research Ltd
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


#include <unistd.h>
#include "pileup.h"
#include "options.h"

#define MIN_MAPQ 0
#define MIN_BASE_QUALITY 30
#define MIN_DEPTH_DEFAULT 2
#define COMPRESSION_LEVEL_DEFAULT 2  // as in the current bgzip call

void Usage() {
  fprintf(stderr, "\nUsage:\n");
  fprintf(stderr, "\t-A\tBulk BAM/CRAM file name\n");
  fprintf(stderr, "\t-B\tDuplex BAM/CRAM file name\n");
  fprintf(stderr, "\t-I\tRegions BED file name\n");
  fprintf(stderr, "\t-C\tSNP BED file name\n");
  fprintf(stderr, "\t-D\tMask BED file name\n");
  fprintf(stderr, "\t-R\tReference sequence file (faidx indexed)\n");
  fprintf(stderr, "\t-Q\tMinimum base quality for bulk sequencing (default %d)\n", MIN_BASE_QUALITY);
  fprintf(stderr, "\t-M\tRemove duplex reads w/ MAPQ smaller than this (default %d)\n", MIN_MAPQ);
  fprintf(stderr, "\t-d\tMinimum duplex depth (default %d)\n", MIN_DEPTH_DEFAULT);
  fprintf(stderr, "\t-O\tOutput file\n");
  fprintf(stderr, "\t-x\tCompression level (default %d)\n", COMPRESSION_LEVEL_DEFAULT);
  fprintf(stderr, "\t-h\tHelp\n");
}

static int SetupOptions(int argc, char **argv, Options *opts) {
  opts->max_plp_depth    = 20000000;
  opts->min_dplx_depth   = MIN_DEPTH_DEFAULT;
  opts->offset           = 1;  // Used to correct genomic positions when comparing to duplex boundaries
  opts->min_base_quality = MIN_BASE_QUALITY;
  opts->min_mapQ         = MIN_MAPQ;
  opts->doTests          = true;
  opts->compression_level = COMPRESSION_LEVEL_DEFAULT;
  opts->debug_mode = false;
  int opt = 0;

  // TODO: make output file mandatory for now?

  while ((opt = getopt(argc, argv, "A:B:I:C:D:R:Q:M:d:O:x:thv")) >= 0) {
    switch (opt) {
      case 'A':
        opts->bams[0] = optarg;
        break;
      case 'B':
        opts->bams[1] = optarg;
        break;
      case 'I':
        opts->ranges_bed = optarg;
        break;
      case 'C':
        opts->beds[MASK_INDEX_SNP] = optarg;
        break;
      case 'D':
        opts->beds[MASK_INDEX_NOISE] = optarg;
        break;
      case 'R':
        opts->fasta = optarg;
        break;
      case 'Q':
        opts->min_base_quality = std::stoi(optarg);
        break;
      case 'M':
        opts->min_mapQ = std::stoi(optarg);
        break;
      case 'd':
        opts->min_dplx_depth = std::stoi(optarg);
        break;
      case 'O':
        opts->oname = optarg;
        break;
      case 'x':
        opts->compression_level = std::stoi(optarg);
        break;
      case 't':
        opts->doTests = false;
        break;
      case 'v':
          opts->debug_mode = true;
          break;
      case 'h':
        Usage();
        exit(0);
      default:
        break;
    }
  }

  return options_validate(opts);
}

int main(int argc, char **argv) {
  Options opts = {};
  if (SetupOptions(argc, argv, &opts)) {
      return 1;
  }

  Pileup pileup;
  pileup.Initiate(&opts);
  pileup.MultiplePileup();
  return 0;
}
