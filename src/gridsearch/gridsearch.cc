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


/*
   Post-processing filter for BED files output by caller
   
   To run:
   ./gridsearch -B /lustre/scratch119/casm/team78/ro4/drseq/54-HpyCH4V-AluI/tables/variants.bed
*/


#include "gridsearch.h"


char ALPH[4] =
  {'A', 'C', 'G', 'T'};

std::map<char, int> INDEX =
  {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};


row_t GridSearch::ParseRow(std::string line) {
  std::string token;
  std::vector<std::string> tokens;
  std::istringstream iss(line);
  while (std::getline(iss, token, '\t')) {
    tokens.push_back(token);
  }
  row_t row;
  row.chrom          = tokens[0];
  row.chrom_beg      = std::stoi(tokens[1]);
  row.context        = tokens[3];
  row.snp            = std::stoi(tokens[4]);
  row.shearwater     = std::stoi(tokens[5]);
  row.bulk_asxs      = std::stof(tokens[6]);
  row.bulk_nm        = std::stof(tokens[7]);
  row.bfwd_A         = std::stoi(tokens[8]);
  row.bfwd_C         = std::stoi(tokens[9]);
  row.bfwd_G         = std::stoi(tokens[10]);
  row.bfwd_T         = std::stoi(tokens[11]);
  row.bfwd_I         = std::stoi(tokens[12]);
  row.brev_A         = std::stoi(tokens[13]);
  row.brev_C         = std::stoi(tokens[14]);
  row.brev_G         = std::stoi(tokens[15]);
  row.brev_T         = std::stoi(tokens[16]);
  row.brev_I         = std::stoi(tokens[17]);
  row.bp_beg         = std::stoi(tokens[18]);
  row.bp_end         = std::stoi(tokens[19]);
  row.bndl_type      = std::stoi(tokens[21]);
  row.dplx_asxs      = stof(tokens[22]);
  row.dplx_clip      = stof(tokens[23]);
  row.dplx_nm        = stof(tokens[24]);
  row.f1r2_A         = std::stoi(tokens[25]);
  row.f1r2_C         = std::stoi(tokens[26]);
  row.f1r2_G         = std::stoi(tokens[27]);
  row.f1r2_T         = std::stoi(tokens[28]);
  row.f1r2_I         = std::stoi(tokens[29]);
  row.f2r1_A         = std::stoi(tokens[30]);
  row.f2r1_C         = std::stoi(tokens[31]);
  row.f2r1_G         = std::stoi(tokens[32]);
  row.f2r1_T         = std::stoi(tokens[33]);
  row.f2r1_I         = std::stoi(tokens[34]);
  row.f1r2_A_Q       = std::stoi(tokens[35]);
  row.f1r2_C_Q       = std::stoi(tokens[36]);
  row.f1r2_G_Q       = std::stoi(tokens[37]);
  row.f1r2_T_Q       = std::stoi(tokens[38]);
  row.f2r1_A_Q       = std::stoi(tokens[39]);
  row.f2r1_C_Q       = std::stoi(tokens[40]);
  row.f2r1_G_Q       = std::stoi(tokens[41]);
  row.f2r1_T_Q       = std::stoi(tokens[42]);
  row.bulk_ppair     = std::stof(tokens[43]);
  row.dplx_ppair     = std::stof(tokens[44]);
  row.bfwd_canonical = row.bfwd_A + row.bfwd_C + row.bfwd_G + row.bfwd_T;
  row.brev_canonical = row.brev_A + row.brev_C + row.brev_G + row.brev_T;
  row.f1r2_canonical = row.f1r2_A + row.f1r2_C + row.f1r2_G + row.f1r2_T;
  row.f2r1_canonical = row.f2r1_A + row.f2r1_C + row.f2r1_G + row.f2r1_T;
  row.bfwd_total     = row.bfwd_canonical + row.bfwd_I;
  row.brev_total     = row.brev_canonical + row.brev_I;
  row.f1r2_total     = row.f1r2_canonical + row.f1r2_I;
  row.f2r1_total     = row.f2r1_canonical + row.f2r1_I;
  row.left           = row.chrom_beg - row.bp_beg + 1;
  row.right          = row.bp_end - row.chrom_beg - 1;
  return row;
}


bool GridSearch::IsMasked(row_t *row) {
  if ((row->snp + row->shearwater) == 0) {
    return false;
  } else {
    return true;
  }
}


void GridSearch::CallDuplex(row_t *row) {
  std::vector<int> f1r2 = {row->f1r2_A, row->f1r2_C, row->f1r2_G, row->f1r2_T};
  std::vector<int> f2r1 = {row->f2r1_A, row->f2r1_C, row->f2r1_G, row->f2r1_T};
  row->f1r2_call = 'N';
  row->f2r1_call = 'N';
  for (int i = 0; i < 4; i++) { 
    if ((f1r2[i]/row->f1r2_canonical) >= this->frac) {
      row->f1r2_call = ALPH[i];
    }
    if ((f2r1[i]/row->f2r1_canonical) >= this->frac) {
      row->f2r1_call = ALPH[i];
    }
  }
}


// mean 5’ clipped bases on duplex reads { 0 : 'no.filter', 1 : '<=0.1', 2 : '==0' }
int GridSearch::DplxClipFilter(row_t *row) {
  int t = 0;
  std::vector<float> parameters = { 0.1, 0 }; 
  for (int j = 0; j < parameters.size(); j++) {
    if (row->dplx_clip <= parameters[j]) {
      t++;
    }
  }
  return t;
}


// AS-XS { 0 : 'no.filter', 1 : '>=0', 2 : '>=50', 3 : '>=100' }
int GridSearch::AlignmentScoreFilter(row_t *row) {
  int t = 0;
  std::vector<int> parameters = { 0, 50, 100 };
  for (int j = 0; j < parameters.size(); j++) {
    if ((row->bulk_asxs >= parameters[j]) &&
        (row->dplx_asxs >= parameters[j])) { 
      t++;
    }
  }
  return t;
}


// NM { 0 : 'no.filter', 1 : '<=5', 2 : '<=3', 3 : '<=2' }
int GridSearch::MismatchFilter(row_t *row) {
  int t = 0;
  std::vector<int> parameters = { 5, 3, 2 };
  for (int j = 0; j < parameters.size(); j++) {
    if ((row->bulk_nm <= parameters[j]) &&
        (row->dplx_nm <= parameters[j])) {
      t++;
    }  
  }
  return t;
}


// matched-normal reads per strand { 0 : 'no.filter', 1 : '>=5', 2 : '>=10'}
int GridSearch::MatchedNormalFilter(row_t *row) {
  int t = 0;
  std::vector<int> parameters = { 5, 10 };
  for (int j = 0; j < parameters.size(); j++) {
    if ((row->bfwd_total >= parameters[j]) &&
        (row->brev_total >= parameters[j])) {
      t++;
    }
  }
  return t;
}


// duplex reads per strand { 0 : '>=2', 1 : '>=3' }
int GridSearch::DuplexFilter(row_t *row) {
  int t = 0;
  std::vector<int> parameters = { 3 };
  for (int j = 0; j < parameters.size(); j++) {
    if ((row->f1r2_total >= parameters[j]) &&
        (row->f2r1_total >= parameters[j])) {
      t++;
    }
  }
  return t;
}


// duplex consensus base quality { 0 : 'no.filter', 1 : '>=30', 2 : '>=60', 3 : '>=90' }
int GridSearch::ConsensusBaseQualityFilter(row_t *row) {
  std::vector<int> f1r2 = { row->f1r2_A_Q, row->f1r2_C_Q, row->f1r2_G_Q, row->f1r2_T_Q };
  std::vector<int> f2r1 = { row->f2r1_A_Q, row->f2r1_C_Q, row->f2r1_G_Q, row->f2r1_T_Q };
  int f1r2_qual = f1r2[INDEX[row->f1r2_call]];
  int f2r1_qual = f2r1[INDEX[row->f2r1_call]];
  int t = 0;
  std::vector<int> parameters = { 30, 60, 90 };
  for (int j = 0; j < parameters.size(); j++) {
    if ((f1r2_qual >= parameters[j]) &&
        (f2r1_qual >= parameters[j])) {
      t++;
    }
  }
  return t;
}


// proportion of reads with an indel { 0 : 'no.filter', 1 : '<=0.05', 2 : '==0' }
int GridSearch::IndelFilter(row_t *row) {
  int t = 0;
  std::vector<float> parameters = { 0.05, 0.0 }; 
  for (int j = 0; j < parameters.size(); j++) {
    if (((row->bfwd_I/row->bfwd_total)  <= parameters[j]) &&
        ((row->brev_I/row->brev_total)  <= parameters[j]) &&
        ((row->f1r2_I/row->f1r2_total) <= parameters[j]) &&
        ((row->f2r1_I/row->f2r1_total) <= parameters[j])) {
      t++;
    }
  }
  return t;
}


// distance from 5' termini { 0 : 'no.filter', 1 : '>=5', 2 : '>=10' }
int GridSearch::FivePrimeTrimFilter(row_t *row) {
  int t = 0;
  std::vector<int> parameters = { 5, 10 };
  for (int j = 0; j < parameters.size(); j++) {
    if ((row->bndl_type == 1) &&
        (row->left >= parameters[j])) {
      t++;
    }
    if ((row->bndl_type == 2) &&
        (row->right >= parameters[j])) {
      t++;
    }
    if ((row->bndl_type == 3) &&
        (row->left >= parameters[j]) && 
        (row->right >= parameters[j])) {
      t++;
    }
  }
  return t;
}


// distance from read 3' { 0 : no filter, 1 : >=5, 2 : >= 10 }
int GridSearch::ThreePrimeTrimFilter(row_t *row) {
  int t = 0;
  std::vector<int> parameters = { 5, 10 };
  for (int j = 0; j < parameters.size(); j++) {
    int max = this->readlen - parameters[j];
    if ((row->bndl_type == 1) &&
        (row->left <= max)) {
      t++;
    }
    if ((row->bndl_type == 2) &&
        (row->right <= max)) {
      t++;
    } 
    if ((row->bndl_type == 3) &&
        (row->left <= max) &&
        (row->right <= max)) {
      t++;
    }
  }
  return t;
}


// proportion of matched-normal reads that are properly-paired { 0 : 'no.filter', 1 : '>=0.95', 2 : '==1' }
int GridSearch::ProperPairFilter(row_t *row) {
  int t = 0;
  std::vector<float> parameters = { 0.95, 1.0 };
  for (int j = 0; j < parameters.size(); j++) {
    if ((row->bulk_ppair >= parameters[j]) &&
        (row->dplx_ppair >= parameters[j])) {
      t++;
    }
  }
  return t;
}


int GridSearch::IsVariant(row_t *row) {
  std::vector<int> bfwd = {row->bfwd_A, row->bfwd_C, row->bfwd_G, row->bfwd_T};
  std::vector<int> brev = {row->brev_A, row->brev_C, row->brev_G, row->brev_T};
  int i = INDEX[row->call];
  if (((bfwd[i]/row->bfwd_canonical) <= this->vaf) &&
      ((brev[i]/row->brev_canonical) <= this->vaf) &&
      (row->context[1] != row->call)) {
    return 1;
  } else {
    return 0;
  }
}


std::string GridSearch::Filter(row_t *row) {
  int t1  = GridSearch::DplxClipFilter(row);
  int t2  = GridSearch::AlignmentScoreFilter(row);
  int t3  = GridSearch::MismatchFilter(row);
  int t4  = GridSearch::MatchedNormalFilter(row);
  int t5  = GridSearch::DuplexFilter(row);
  int t6  = GridSearch::ConsensusBaseQualityFilter(row);
  int t7  = GridSearch::IndelFilter(row);
  int t8  = GridSearch::FivePrimeTrimFilter(row);
  int t9  = GridSearch::ThreePrimeTrimFilter(row);
  int t10 = GridSearch::ProperPairFilter(row);
  int t11 = GridSearch::IsVariant(row);
  std::stringstream ss;
  ss << t1;
  ss << ",";
  ss << t2;
  ss << ",";
  ss << t3;
  ss << ",";
  ss << t4;
  ss << ",";
  ss << t5;
  ss << ",";
  ss << t6;
  ss << ",";
  ss << t7;
  ss << ",";
  ss << t8;
  ss << ",";
  ss << t9;
  ss << ",";
  ss << t10;
  ss << ",";
  ss << t11;
  return ss.str();
}


void GridSearch::CollectMetrics() {
  this->gzin.open(this->bed);
  if ( ! this->gzin) {
    std::stringstream er;
    er << "Error: failed to open bed file " << this->bed << std::endl;
    throw std::runtime_error(er.str());
  }
  std::string line;
  while (std::getline(this->gzin, line)) {
    if (line[0] != '#') {
      row_t row = GridSearch::ParseRow(line);
      if ((GridSearch::IsMasked(&row) == true) ||
          (row.bfwd_canonical == 0) ||
          (row.brev_canonical == 0) ||
          (row.f1r2_canonical == 0) ||
          (row.f2r1_canonical == 0)) {
        continue;
      }
      GridSearch::CallDuplex(&row);
      if ((row.f1r2_call == 'N') ||
          (row.f2r1_call == 'N') || 
          (row.f1r2_call != row.f2r1_call)) {
        continue;
      }
      row.call = row.f1r2_call;
      this->counts[GridSearch::Filter(&row)]++;
    }
  }
}


void GridSearch::WriteMetrics() {
  for (auto f1 : this->counts) {
    std::cout << f1.first;
    std::cout << ",";
    std::cout << f1.second;
    std::cout << std::endl;
  }
}


void Usage() {
  fprintf(stderr, "\nUsage:\n");
  fprintf(stderr, "\t-B\tBED input file name\n");
  fprintf(stderr, "\t-h\tHelp\n\n");
}


void Options(int argc, char **argv, GridSearch *vc) {
  vc->bed     = NULL;
  vc->frac    = 0.9;
  vc->readlen = 151;
  vc->vaf     = 0.01;
  int opt = 0;
  while ((opt = getopt(argc, argv, "B:f:r:v:h")) >= 0) {
    switch (opt) {
      case 'B':
        vc->bed = optarg;
        break;
      case 'f':
        vc->frac = std::stof(optarg);
        break;
      case 'r':
        vc->readlen = std::stoi(optarg);
        break;
      case 'v':
        vc->vaf = std::stof(optarg);
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
  GridSearch vc;
  Options(argc, argv, &vc);
  vc.CollectMetrics();
  vc.WriteMetrics();
  return 0;
}

