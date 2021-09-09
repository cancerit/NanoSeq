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

#include "./pileup.h"
#include "./options.h"
char buffer[400];

void Usage() {
  fprintf(stderr, "\nUsage:\n");
  fprintf(stderr, "\t-A\tBulk BAM file name\n");
  fprintf(stderr, "\t-B\tDuplex BAM file name\n");
  fprintf(stderr, "\t-C\tSNP BED file name\n");
  fprintf(stderr, "\t-D\tMask BED file name\n");
  fprintf(stderr, "\t-R\tReference sequence file (faidx indexed)\n");
  fprintf(stderr, "\t-Q\tMinimum base quality for bulk sequencing (def 30)\n");
  fprintf(stderr, "\t-M\tRemove duplex reads w/ MAPQ smaller than this (def 0)\n"); 
  fprintf(stderr, "\t-r\tReference or contig name\n");
  fprintf(stderr, "\t-b\tStart coordinate\n");
  fprintf(stderr, "\t-e\tEnd coordinate\n");
  fprintf(stderr, "\t-d\tMinimum duplex depth (default 2)\n");
  fprintf(stderr, "\t-O\tOutput file\n");
  fprintf(stderr, "\t-h\tHelp\n");
}


static void SetupOptions(int argc, char **argv, Options *opts) {
  opts->max_plp_depth    = 1000000;
  opts->min_dplx_depth   = 2;
  opts->offset           = 1;
  opts->min_base_quality = 30;
  opts->min_mapQ         = 0;
  opts->out2stdout       = true;
  opts->testBulk         = true;
  char suffix[] = ".gz";
  int opt = 0;
  while ((opt = getopt(argc, argv, "A:B:C:D:R:Q:M:r:b:e:d:O:th")) >= 0) {
    switch (opt) {
      case 'A':
        opts->bams[0] = optarg;
        break;
      case 'B':
        opts->bams[1] = optarg;
        break;
      case 'C':
        opts->beds[0] = optarg;
        break;
      case 'D':
        opts->beds[1] = optarg;
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
      case 'r':
        opts->rname = optarg;
        break;
      case 'b':
        opts->beg = std::stoi(optarg);
        break;
      case 'e':
        opts->end = std::stoi(optarg);
        break;
      case 'd':
        opts->min_dplx_depth = std::stoi(optarg);
        break;
      case 'O':
        strcpy(buffer,optarg);
        strcat(buffer,suffix);
        opts->oname = buffer;
        opts->out2stdout = false;
        break;
      case 't':
        opts->testBulk = false;
        break;
      case 'h':
        Usage();
        exit(0);
      default:
        break;
    }
  }
}


int main(int argc, char **argv) {
  Options opts;
  SetupOptions(argc, argv, &opts);
  Pileup pileup;
  pileup.Initiate(&opts);
  pileup.MultiplePileup();
  return 0;
}
