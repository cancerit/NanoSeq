/*########## LICENCE ##########
# Copyright (c) 2022 Genome Research Ltd
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
#include "bamaddreadbundles.h"


void BamAddReadBundles::LoadFiles() {
  this->in = hts_open(this->infile, "r");
  if (this->in == 0) {
    std::stringstream er;
    er << "Fail to open input BAM/CRAM file ";
    er << this->infile;
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
  this->head = sam_hdr_read(this->in);
  if (this->head == NULL || this->head->n_targets == 0) {
    std::stringstream er;
    er << "Error: BAM/CRAM file does not have header.";
    er << std::endl;
    throw std::runtime_error(er.str());
  }
  
  char * mode = sam_open_mode_opts( this->outfile,"w",NULL);
  if ( mode == NULL ) {
      std::stringstream er;
      er << "Error : Output extension must be .bam or .cram";
      er << std::endl;
      throw std::runtime_error(er.str());
  }
  this->out = hts_open(this->outfile, mode);

  if (!this->out) {
    std::stringstream er;
    er << "Error: failed to open ";
    er << this->outfile;
    er << " for output.";
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
}


void BamAddReadBundles::UpdateHeader(int argc, char **argv) {
  // command line argument
  std::stringstream ss;
  int i;
  for (i = 0; i < (argc-1); ++i) {
    ss << argv[i];
    ss << " ";
  }
  ss << argv[i];
  sam_hdr_t *sh = sam_hdr_parse(this->head->l_text,this->head->text);
  if (  sh == NULL ) exit(1);
  char pg[]   = "bamaddreadbundles";
  int rco = sam_hdr_add_pg(sh, pg, "VN", "1.0", "CL", ss.str().c_str(), NULL);
  if ( rco < 0 ) exit(1);
  free(this->head->text);
  rco = sam_hdr_rebuild(sh);
  if ( rco < 0 ) exit(1);
  this->head->text   = strdup(sam_hdr_str(sh));
  this->head->l_text = sam_hdr_length(sh);
  sam_hdr_destroy(sh);
  if (sam_hdr_write(this->out, this->head) < 0) {
      std::stringstream er;
      er << "Failed to write header";
      er << std::endl;
      throw std::runtime_error(er.str());
  }
}


bool BamAddReadBundles::HasAux(bam1_t* b, const char* tag) {
  return (bam_aux_get(b, tag))? true : false;
}


bool BamAddReadBundles::ReadIsUsable(bam1_t* b) {
  int suppl     = (b->core.flag & BAM_FSUPPLEMENTARY)? 1 : 0;
  int qcfail    = (b->core.flag & BAM_FQCFAIL)? 1 : 0;
  int unmapped  = (b->core.flag & BAM_FUNMAP)? 1 : 0;
  int secondary = (b->core.flag & BAM_FSECONDARY)? 1 : 0;
  int has_od    = (BamAddReadBundles::HasAux(b, "od"))? 1 : 0;
  int paired    = (b->core.flag & BAM_FPROPER_PAIR)? 1 : 0;
  int has_rc    = (BamAddReadBundles::HasAux(b, "rc"))? 1 : 0;
  int has_mc    = (BamAddReadBundles::HasAux(b, "mc"))? 1 : 0;
  int has_rb    = (BamAddReadBundles::HasAux(b, "rb"))? 1 : 0;
  int has_mb    = (BamAddReadBundles::HasAux(b, "mb"))? 1 : 0;
  if (((suppl + qcfail + unmapped + secondary + has_od) == 0) &&
      ((paired  + has_rc + has_mc + has_rb + has_mb) == 5)) {
    return true;
  }
  return false;
}

bool BamAddReadBundles::ReadIsWritable(bam1_t* b) {
  int qcfail    = (b->core.flag & BAM_FQCFAIL)? 1 : 0;
  int has_RB    = (BamAddReadBundles::HasAux(b, "RB"))? 1 : 0;
  int has_od    = (BamAddReadBundles::HasAux(b, "od"))? 1 : 0;
  int has_rb    = (BamAddReadBundles::HasAux(b, "rb"))? 1 : 0;
  int has_mb    = (BamAddReadBundles::HasAux(b, "mb"))? 1 : 0;
  if ( has_RB == 1) {
    return true;
  }
  if ((( has_rb + has_mb) == 2) && ((qcfail + has_od) == 0 )) {
    return true;
  }
  return false;
}

void BamAddReadBundles::AddAuxTags(bam1_t* b) {
  int rc     = bam_aux2i(bam_aux_get(b, "rc"));
  int mc     = bam_aux2i(bam_aux_get(b, "mc"));
  char* rb   = bam_aux2Z(bam_aux_get(b, "rb"));
  char* mb   = bam_aux2Z(bam_aux_get(b, "mb"));
  int strand = (b->core.flag & BAM_FREVERSE)? 1: 0;
  std::stringstream ss;
  ss << this->head->target_name[b->core.tid];
  ss << ",";
  ss << std::min(mc, rc);
  ss << ",";
  ss << std::max(mc, rc);
  ss << ",";
  if (strand == 0) {
    ss << rb;
    ss << ",";
    ss << mb;
  } else {
    ss << mb;
    ss << ",";
    ss << rb;
  }
  std::string str = ss.str();
  const char* cstr = str.c_str();
  int len = str.length() + 1;
  uint8_t *data = const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(cstr));
  int rco = bam_aux_append(b, "RB", 'Z', len, data);
  if ( rco < 0 ) exit(1);
}


void BamAddReadBundles::DelAuxTags(bam1_t* b) {
  uint8_t* t_tag;
  if (( t_tag = bam_aux_get(b, "mc")) != NULL  ) {
    bam_aux_del(b, t_tag);
  }
  if (( t_tag = bam_aux_get(b, "rc")) != NULL  ) {
    bam_aux_del(b, t_tag);
  }
  if (( t_tag = bam_aux_get(b, "mb")) != NULL  ) {
    bam_aux_del(b, t_tag);
  }
  if (( t_tag = bam_aux_get(b, "rb")) != NULL  ) {
    bam_aux_del(b, t_tag);
  }
  if (( t_tag = bam_aux_get(b, "MQ")) != NULL  ) {
    bam_aux_del(b, t_tag);
  }
  if (( t_tag = bam_aux_get(b, "ms")) != NULL  ) {
    bam_aux_del(b, t_tag);
  }
  if (( t_tag = bam_aux_get(b, "MC")) != NULL  ) {
    bam_aux_del(b, t_tag);
  }
}


void BamAddReadBundles::WriteOut(bam1_t* b) {
  if ((sam_write1(this->out, this->head, b) < 0)) {
    std::stringstream er;
    er << "Error: failed to write record.";
    er << std::endl;
    throw std::runtime_error(er.str());
  }
}


void BamAddReadBundles::FilterAndTagReads() {
  bam1_t *b = bam_init1();
  int ret;
  while ( 1 ) {
    ret = sam_read1(this->in, this->head, b);
    if ( ret == -1 ) break;
    if ( ret < -1 ) {
      std::stringstream er;
      er << "Error: failure while reading input BAM";
      er << std::endl;
      throw std::runtime_error(er.str());
    }
    if (BamAddReadBundles::ReadIsUsable(b) == true) {
      BamAddReadBundles::AddAuxTags(b);
      BamAddReadBundles::DelAuxTags(b);
    }
    if ( BamAddReadBundles::ReadIsWritable(b) == true ) {
      BamAddReadBundles::WriteOut(b);
    }
  }
  bam_destroy1(b);
}


void BamAddReadBundles::CleanUp() {
  bam_hdr_destroy(this->head);
  hts_close(this->in);
  if (sam_close(this->out) < 0) {
    std::stringstream er;
    er << "Error closing output file.";
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
}


void Usage() {
  fprintf(stderr, "\nUsage:\n");
  fprintf(stderr, "\t-I\tInput BAM/CRAM file name\n");
  fprintf(stderr, "\t-O\tOutput BAM/CRAM file name\n");
  fprintf(stderr, "\t-h\tHelp\n");
}


int main(int argc, char **argv) {
  BamAddReadBundles barb;
  barb.infile = NULL;
  barb.outfile = NULL;
  int opt = 0;
  while ((opt = getopt(argc, argv, "I:O:h")) >= 0) {
    switch (opt) {
      case 'I':
        barb.infile = optarg;
        break;
      case 'O':
        barb.outfile = optarg;
        break;
      case 'h':
        Usage();
        exit(0);
      default:
        break;
    }
  }
  if (barb.infile == NULL) {
    std::stringstream er;
    er << "Error: no input file specified";
    throw std::invalid_argument(er.str());
  }
  if (barb.outfile == NULL) {
    std::stringstream er;
    er << "Error: no output file specified";
    throw std::invalid_argument(er.str());
  }
  barb.LoadFiles();
  barb.UpdateHeader(argc, argv);
  barb.FilterAndTagReads();
  barb.CleanUp();
  return 0;
}
