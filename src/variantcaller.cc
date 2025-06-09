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

#include "variantcaller.h"

char ALPH[4] =
  {'A', 'C', 'G', 'T'};

std::map<char, int> INDEX =
  {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

std::map<char, char> COMPLEMENT = 
  {{'A','T'}, {'C','G'}, {'G','C'}, {'T','A'}};


row_t VariantCaller::ParseRow(std::string line) {
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
  row.dplx_barcode   = tokens[20];
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
  row.min_qpos       = std::min(row.left, row.right);
  return row;
}


void VariantCaller::CallDuplex(row_t *row) {
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


int VariantCaller::DplxClipFilter(row_t *row) {
  if (row->dplx_clip <= this->clip) {
    return 1;
  } else {
    return 0;
  }
}


int VariantCaller::AlignmentScoreFilter(row_t *row) {
  if ((row->bulk_asxs >= this->asxs) &&
      (row->dplx_asxs >= this->asxs)) { 
    return 1;
  } else {
    return 0;
  }
}


int VariantCaller::MismatchFilter(row_t *row) {
 // fa8 method to check the NM filter is ok after counting the mutation
  if(row->isvariant) {
    if ((row->bulk_nm <= this->nmms) &&
        (row->dplx_nm <= this->nmms+1)) { //if it is variant we should not count the mutation as a mismatch
      return 1;
    } else {
      return 0;
    }
  } else {
    if ((row->bulk_nm <= this->nmms) &&
        (row->dplx_nm <= this->nmms)) {
      return 1;
    } else {
      return 0;
    }
  }
}


int VariantCaller::MatchedNormalFilter(row_t *row) {
  if ((row->bfwd_total >= this->bulk) &&
      (row->brev_total >= this->bulk) &&
      (row->brev_total + row->bfwd_total >= this->bulk_total)) { // Modified by fa8, to include bulk_total condition
    return 1;
  } else {
    return 0;
  }
}


int VariantCaller::DuplexFilter(row_t *row) {
  if ((row->f1r2_total >= this->dplx) &&
      (row->f2r1_total >= this->dplx)) {
    return 1;
  } else {
    return 0;
  }
}


int VariantCaller::ConsensusBaseQualityFilter(row_t *row) {
  std::vector<int> f1r2 = { row->f1r2_A_Q, row->f1r2_C_Q, row->f1r2_G_Q, row->f1r2_T_Q };
  std::vector<int> f2r1 = { row->f2r1_A_Q, row->f2r1_C_Q, row->f2r1_G_Q, row->f2r1_T_Q };
  int f1r2_qual = f1r2[INDEX[row->f1r2_call]];
  int f2r1_qual = f2r1[INDEX[row->f2r1_call]];
  if ((f1r2_qual >= this->qual) &&
      (f2r1_qual >= this->qual)) {
    return 1;
  } else {
    return 0;
  }
}


int VariantCaller::IndelFilter(row_t *row) {
  if(row->bfwd_total > 0 && row->bfwd_I/row->bfwd_total > this->indel) {
  	return 0;
  }
  if(row->brev_total > 0 && row->brev_I/row->brev_total > this->indel) {
  	return 0;
  }
  if(row->f1r2_total > 0 && row->f1r2_I/row->f1r2_total > this->indel) {
  	return 0;
  }
  if(row->f2r1_total > 0 && row->f2r1_I/row->f2r1_total > this->indel) {
  	return 0;
  }
  return 1;
}


int VariantCaller::FivePrimeTrimFilter(row_t *row) {
  if (row->bndl_type == 1) {
    if (row->left >= this->min_cycle) {
      return 1;
    } else {
      return 0;
    }
  } else if (row->bndl_type == 2) {
    if (row->right >= this->min_cycle) {
      return 1;
    } else {
      return 0;
    }
  } else {
    assert(row->bndl_type == 3);
    if ((row->left >= this->min_cycle) && 
        (row->right >= this->min_cycle)) {
      return 1;
    } else {
      return 0;
    }
  }
}


int VariantCaller::ThreePrimeTrimFilter(row_t *row) {
  int max = this->readlen - this->max_cycle;
  if (row->bndl_type == 1) {
    if (row->left <= max) {
      return 1;
    } else {
      return 0;
    }
  } else if (row->bndl_type == 2) {
    if (row->right <= max) {
      return 1;
    } else {
      return 0;
    }
  } else {
    assert(row->bndl_type == 3);
    if ((row->left <= max) &&
        (row->right <= max)) {
      return 1;
    } else {
      return 0;
    }
  }
}


int VariantCaller::ProperPairFilter(row_t *row) {
  if ((row->bulk_ppair >= this->ppair) &&
      (row->dplx_ppair >= this->ppair)) {
    return 1;
  } else {
    return 0;
  }
}


void VariantCaller::ApplyFilters(row_t *row) {
  int t1  = VariantCaller::DplxClipFilter(row);
  int t2  = VariantCaller::AlignmentScoreFilter(row);
  int t3  = VariantCaller::MismatchFilter(row);
  int t4  = VariantCaller::MatchedNormalFilter(row);
  int t5  = VariantCaller::DuplexFilter(row);
  int t6  = VariantCaller::ConsensusBaseQualityFilter(row);
  int t7  = VariantCaller::IndelFilter(row);
  int t8  = VariantCaller::FivePrimeTrimFilter(row);
  int t9  = VariantCaller::ThreePrimeTrimFilter(row);
  int t10 = VariantCaller::ProperPairFilter(row);
  bool pass = (t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10) == 10;

  // Store filter results in row
  row->dplx_clip_filter = t1;
  row->alignment_score_filter = t2;
  row->mismatch_filter = t3;
  row->matched_normal_filter = t4;
  row->duplex_filter = t5;
  row->consensus_base_quality_filter = t6;
  row->indel_filter = t7;
  row->five_prime_trim_filter = t8;
  row->three_prime_trim_filter = t9;
  row->proper_pair_filter = t10;
  row->pass_all_filters = pass;
}

int VariantCaller::VafFilter(row_t *row) {
  std::vector<int> bfwd = {row->bfwd_A, row->bfwd_C, row->bfwd_G, row->bfwd_T};
  std::vector<int> brev = {row->brev_A, row->brev_C, row->brev_G, row->brev_T};
  if(row->f1r2_call == row->f2r1_call) { // Check if both strands consistent
    row->call = row->f1r2_call;
  }
  int i = INDEX[row->call];
  if(row->bfwd_canonical > 0 && bfwd[i]/row->bfwd_canonical > this->vaf ) {
    return 0;
  } else if(row->brev_canonical > 0 && brev[i]/row->brev_canonical > this->vaf ) {
    return 0;
  } else {
	  return 1;
  }
}

int VariantCaller::IsVariant(row_t *row) {
  if(row->context[1] != row->call) {
    return 1;  
  } else {
    return 0;
  }
}


int VariantCaller::IsMasked(row_t *row) {
  if ((row->snp + row->shearwater) == 0) {
    return 0;
  } else {
    return 1;
  }
}


std::string VariantCaller::ReverseComplementContext(std::string context) {
  std::stringstream ss;
  ss << COMPLEMENT[context[2]];
  ss << COMPLEMENT[context[1]];
  ss << COMPLEMENT[context[0]];
  return ss.str();
}


std::string VariantCaller::PyrimidineContext(row_t *row) {
  if ((row->context[1] == 'A') ||
      (row->context[1] == 'G')) {
    return VariantCaller::ReverseComplementContext(row->context);
  } else {
    return row->context;
  }
}


std::string VariantCaller::PyrimidineSubstitution(row_t *row) {
  if ((row->context[1] == 'A') ||
      (row->context[1] == 'G')) {
    std::stringstream ss;
    ss << VariantCaller::ReverseComplementContext(row->context);
    ss << ">";
    ss << COMPLEMENT[row->call];
    return ss.str();
  } else {
    std::stringstream ss;
    ss << row->context;
    ss << ">";
    ss << row->call;
    return ss.str();
  }
}


std::string VariantCaller::StrandContext(row_t *row) {
  if (row->left < row->right) {
    return row->context;
  } else if (row->right < row->left) {
    return VariantCaller::ReverseComplementContext(row->context);
  } else {
    return "UNKNOWN";
  }
}


std::string VariantCaller::StrandSubstitution(row_t *row) {
  if (row->left < row->right) {
    std::stringstream ss;
    ss << row->context;
    ss << ">";
    ss << row->call;
    return ss.str();
  } else if (row->right < row->left) {
    std::stringstream ss;
    ss << VariantCaller::ReverseComplementContext(row->context);
    ss << ">";
    ss << COMPLEMENT[row->call];
    return ss.str();
  } else {
    return "UNKNOWN";
  }
}


std::string VariantCaller::StrandMismatch(row_t *row) {
  std::string mismatch;
  if (row->f1r2_call == row->context[1]) {
    // error on f2r1 strand
    std::stringstream ss;
    ss << VariantCaller::ReverseComplementContext(row->context);
    ss << ">";
    ss << COMPLEMENT[row->f2r1_call];
    mismatch = ss.str();
  } else if (row->f2r1_call == row->context[1]) {
    // error on f1r2 strand
    std::stringstream ss;
    ss << row->context;
    ss << ">";
    ss << row->f1r2_call;
    mismatch = ss.str();
  } else {
    mismatch = "UNKNOWN";
  }
  return mismatch;
}


bool VariantCaller::ContextIsCanonical(row_t *row) {
  int t = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      if (ALPH[i] == row->context[j]) {
        t++;
      }
    }
  }
  return (t == 3)? true : false;
}


void VariantCaller::CollectMetrics() {
  this->gzin.open(this->bed);
  if ( ! this->gzin ) {
    std::stringstream er;
    er << "Error: failed to open input file " << this->bed << std::endl;
    throw std::runtime_error(er.str());
  }
  this->coverage = 0;
  int curr = -1;
  std::string line;
  row_t lastRow;
  bool first = true;
  int cov = 0;
  while (std::getline(this->gzin, line)) {
    if (line[0] != '#') {
      row_t row = VariantCaller::ParseRow(line);
      // Modified by fa8 to not check the bulk
      if ((row.f1r2_canonical == 0) ||
          (row.f2r1_canonical == 0)) {
        continue;
      }
      VariantCaller::CallDuplex(&row);
      if ((row.f1r2_call == 'N') ||
          (row.f2r1_call == 'N')) {
        continue;
      }
      if (VariantCaller::ContextIsCanonical(&row) == false) {
        continue;
      }

      if (row.f1r2_call == row.f2r1_call) { // Check if both strands consistent
        //fa8: moved this so we know if it is variant or not before applying the NM filter:
        row.call = row.f1r2_call;
        row.isvariant  = VariantCaller::IsVariant(&row);
        row.ismasked   = VariantCaller::IsMasked(&row);
        row.pyrcontext = VariantCaller::PyrimidineContext(&row);
        VariantCaller::ApplyFilters(&row); // Apply filters to row
        if(row.isvariant && row.f1r2_call != row.context[1] && row.f2r1_call != row.context[1]) { // fa8: these conditions are redundant
	        row.vaf_filter = VariantCaller::VafFilter(&row); // fa8: This one has to go separately 
    	                                                     // from the other filters
      	} else {
      		row.vaf_filter = 1;
      	}
        if (row.pass_all_filters && row.vaf_filter) {
          this->burdens[row.ismasked][row.isvariant]++;
          this->call_by_qpos[row.call][row.min_qpos][row.ismasked]++;
          this->pyr_by_mask[row.pyrcontext][row.ismasked]++;
          this->read_bundles[row.f1r2_total][row.f2r1_total][row.ismasked][row.isvariant]++;
          if (row.isvariant) {
            VariantCaller::WriteVariants(&row);
          }
          //fa8: write the information needed for analysing coverage:
          // only for not masked sites!
          if(this->outfile_coverage != NULL) {
            if(row.shearwater == 0 && row.snp == 0) {
              if ( first ) {
                cov = 1;
                first = false;
              } else if ( row.chrom == lastRow.chrom && row.chrom_beg ==  lastRow.chrom_beg && row.context == lastRow.context) {
                cov++;
              } else {
                this->gzout_coverage << lastRow.chrom;
                this->gzout_coverage << "\t";
                this->gzout_coverage << lastRow.chrom_beg;
                this->gzout_coverage << "\t";
                this->gzout_coverage << lastRow.chrom_beg+1;
                this->gzout_coverage << "\t";
                this->gzout_coverage << lastRow.context << ";" << lastRow.context[1] << ";" << cov;
                this->gzout_coverage << std::endl;
                cov = 1;
              }
              lastRow = row;
            }
          }
          if (row.chrom_beg != curr) {
            this->coverage++;
            curr = row.chrom_beg;
          }
        }  else {  // retain variants that fail filters
          if (row.isvariant and this->outfile_discarded != NULL) {
            VariantCaller::WriteDiscardedVariants(&row);
          }
        }
      } else { // this is for strand-specific errors (where f1r2 != f2r1)
        // row.isvariant = 0;
        VariantCaller::ApplyFilters(&row); // Apply filters to row
        if (row.pass_all_filters) {
          VariantCaller::WriteMismatches(&row);
        }
      }
    }
  }
  //write last line to coverage file
  if(this->outfile_coverage != NULL and curr != -1 and cov > 0 ) {
    this->gzout_coverage << lastRow.chrom;
    this->gzout_coverage << "\t";
    this->gzout_coverage << lastRow.chrom_beg;
    this->gzout_coverage << "\t";
    this->gzout_coverage << lastRow.chrom_beg+1;
    this->gzout_coverage << "\t";
    this->gzout_coverage << lastRow.context << ";" << lastRow.context[1] << ";" << cov;
    this->gzout_coverage << std::endl;
  }
}


void VariantCaller::WriteVariants(row_t *row) {
  this->fout << "Variants";
  this->fout << "\t";
  this->fout << row->chrom;
  this->fout << "\t";
  this->fout << row->chrom_beg;
  this->fout << "\t";
  this->fout << row->context;
  this->fout << "\t";
  this->fout << row->snp;
  this->fout << "\t";
  this->fout << row->shearwater;
  this->fout << "\t";
  this->fout << row->bulk_asxs;
  this->fout << "\t";
  this->fout << row->bulk_nm;
  this->fout << "\t";
  this->fout << row->bfwd_A;
  this->fout << "\t";
  this->fout << row->bfwd_C;
  this->fout << "\t";
  this->fout << row->bfwd_G;
  this->fout << "\t";
  this->fout << row->bfwd_T;
  this->fout << "\t";
  this->fout << row->bfwd_I;
  this->fout << "\t";
  this->fout << row->brev_A;
  this->fout << "\t";
  this->fout << row->brev_C;
  this->fout << "\t";
  this->fout << row->brev_G;
  this->fout << "\t";
  this->fout << row->brev_T;
  this->fout << "\t";
  this->fout << row->brev_I;
  this->fout << "\t";
  this->fout << row->bp_beg;
  this->fout << "\t";
  this->fout << row->bp_end;
  this->fout << "\t";
  this->fout << row->bndl_type;
  this->fout << "\t";
  this->fout << row->dplx_asxs;
  this->fout << "\t";
  this->fout << row->dplx_clip;
  this->fout << "\t";
  this->fout << row->dplx_nm;
  this->fout << "\t";
  this->fout << row->f1r2_A;
  this->fout << "\t";
  this->fout << row->f1r2_C;
  this->fout << "\t";
  this->fout << row->f1r2_G;
  this->fout << "\t";
  this->fout << row->f1r2_T;
  this->fout << "\t";
  this->fout << row->f1r2_I;
  this->fout << "\t";
  this->fout << row->f2r1_A;
  this->fout << "\t";
  this->fout << row->f2r1_C;
  this->fout << "\t";
  this->fout << row->f2r1_G;
  this->fout << "\t";
  this->fout << row->f2r1_T;
  this->fout << "\t";
  this->fout << row->f2r1_I;
  this->fout << "\t";
  this->fout << row->f1r2_A_Q;
  this->fout << "\t";
  this->fout << row->f1r2_C_Q;
  this->fout << "\t";
  this->fout << row->f1r2_G_Q;
  this->fout << "\t";
  this->fout << row->f1r2_T_Q;
  this->fout << "\t";
  this->fout << row->f2r1_A_Q;
  this->fout << "\t";
  this->fout << row->f2r1_C_Q;
  this->fout << "\t";
  this->fout << row->f2r1_G_Q;
  this->fout << "\t";
  this->fout << row->f2r1_T_Q;
  this->fout << "\t";
  this->fout << row->bfwd_total;
  this->fout << "\t";
  this->fout << row->brev_total;
  this->fout << "\t";
  this->fout << row->f1r2_total;
  this->fout << "\t";
  this->fout << row->f2r1_total;
  this->fout << "\t";
  this->fout << row->left;
  this->fout << "\t";
  this->fout << row->right;
  this->fout << "\t";
  this->fout << row->min_qpos;
  this->fout << "\t";
  this->fout << row->call;
  this->fout << "\t";
  this->fout << row->isvariant;
  this->fout << "\t";
  this->fout << row->pyrcontext;
  this->fout << "\t";
  this->fout << VariantCaller::StrandContext(row);
  this->fout << "\t";
  this->fout << VariantCaller::PyrimidineSubstitution(row);
  this->fout << "\t";
  this->fout << VariantCaller::StrandSubstitution(row);
  this->fout << "\t";
  this->fout << row->ismasked;
  this->fout << "\t";
  this->fout << row->dplx_barcode;
  this->fout << std::endl;
}

// For the moment this is a copy of the above function
// Should write a general WriteVariant function and change the outstream for the specific destination
void VariantCaller::WriteDiscardedVariants(row_t *row) {
  this->fout_discarded << "DiscardedVariants";
  this->fout_discarded << "\t";
  this->fout_discarded << row->chrom;
  this->fout_discarded << "\t";
  this->fout_discarded << row->chrom_beg;
  this->fout_discarded << "\t";
  this->fout_discarded << row->context;
  this->fout_discarded << "\t";
  this->fout_discarded << row->snp;
  this->fout_discarded << "\t";
  this->fout_discarded << row->shearwater;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bulk_asxs;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bulk_nm;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bfwd_A;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bfwd_C;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bfwd_G;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bfwd_T;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bfwd_I;
  this->fout_discarded << "\t";
  this->fout_discarded << row->brev_A;
  this->fout_discarded << "\t";
  this->fout_discarded << row->brev_C;
  this->fout_discarded << "\t";
  this->fout_discarded << row->brev_G;
  this->fout_discarded << "\t";
  this->fout_discarded << row->brev_T;
  this->fout_discarded << "\t";
  this->fout_discarded << row->brev_I;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bp_beg;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bp_end;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bndl_type;
  this->fout_discarded << "\t";
  this->fout_discarded << row->dplx_asxs;
  this->fout_discarded << "\t";
  this->fout_discarded << row->dplx_clip;
  this->fout_discarded << "\t";
  this->fout_discarded << row->dplx_nm;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_A;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_C;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_G;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_T;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_I;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_A;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_C;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_G;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_T;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_I;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_A_Q;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_C_Q;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_G_Q;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_T_Q;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_A_Q;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_C_Q;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_G_Q;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_T_Q;
  this->fout_discarded << "\t";
  this->fout_discarded << row->bfwd_total;
  this->fout_discarded << "\t";
  this->fout_discarded << row->brev_total;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f1r2_total;
  this->fout_discarded << "\t";
  this->fout_discarded << row->f2r1_total;
  this->fout_discarded << "\t";
  this->fout_discarded << row->left;
  this->fout_discarded << "\t";
  this->fout_discarded << row->right;
  this->fout_discarded << "\t";
  this->fout_discarded << row->min_qpos;
  this->fout_discarded << "\t";
  this->fout_discarded << row->call;
  this->fout_discarded << "\t";
  this->fout_discarded << row->isvariant;
  this->fout_discarded << "\t";
  this->fout_discarded << row->pyrcontext;
  this->fout_discarded << "\t";
  this->fout_discarded << VariantCaller::StrandContext(row);
  this->fout_discarded << "\t";
  this->fout_discarded << VariantCaller::PyrimidineSubstitution(row);
  this->fout_discarded << "\t";
  this->fout_discarded << VariantCaller::StrandSubstitution(row);
  this->fout_discarded << "\t";
  this->fout_discarded << row->ismasked;
  this->fout_discarded << "\t";
  this->fout_discarded << row->dplx_barcode;
  this->fout_discarded << "\t";
  this->fout_discarded << row->dplx_clip_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->alignment_score_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->mismatch_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->matched_normal_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->duplex_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->consensus_base_quality_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->indel_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->five_prime_trim_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->three_prime_trim_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->proper_pair_filter;
  this->fout_discarded << "\t";
  this->fout_discarded << row->vaf_filter;
  this->fout_discarded << std::endl;
}


void VariantCaller::WriteMismatches(row_t *row) {
  this->fout << "Mismatches";
  this->fout << "\t";
  this->fout << row->chrom;
  this->fout << "\t";
  this->fout << row->chrom_beg;
  this->fout << "\t";
  this->fout << row->context;
  this->fout << "\t";
  this->fout << row->snp;
  this->fout << "\t";
  this->fout << row->shearwater;
  this->fout << "\t";
  this->fout << row->bulk_asxs;
  this->fout << "\t";
  this->fout << row->bulk_nm;
  this->fout << "\t";
  this->fout << row->bfwd_A;
  this->fout << "\t";
  this->fout << row->bfwd_C;
  this->fout << "\t";
  this->fout << row->bfwd_G;
  this->fout << "\t";
  this->fout << row->bfwd_T;
  this->fout << "\t";
  this->fout << row->bfwd_I;
  this->fout << "\t";
  this->fout << row->brev_A;
  this->fout << "\t";
  this->fout << row->brev_C;
  this->fout << "\t";
  this->fout << row->brev_G;
  this->fout << "\t";
  this->fout << row->brev_T;
  this->fout << "\t";
  this->fout << row->brev_I;
  this->fout << "\t";
  this->fout << row->bp_beg;
  this->fout << "\t";
  this->fout << row->bp_end;
  this->fout << "\t";
  this->fout << row->dplx_asxs;
  this->fout << "\t";
  this->fout << row->dplx_clip;
  this->fout << "\t";
  this->fout << row->dplx_nm;
  this->fout << "\t";
  this->fout << row->f1r2_A;
  this->fout << "\t";
  this->fout << row->f1r2_C;
  this->fout << "\t";
  this->fout << row->f1r2_G;
  this->fout << "\t";
  this->fout << row->f1r2_T;
  this->fout << "\t";
  this->fout << row->f1r2_I;
  this->fout << "\t";
  this->fout << row->f2r1_A;
  this->fout << "\t";
  this->fout << row->f2r1_C;
  this->fout << "\t";
  this->fout << row->f2r1_G;
  this->fout << "\t";
  this->fout << row->f2r1_T;
  this->fout << "\t";
  this->fout << row->f2r1_I;
  this->fout << "\t";
  this->fout << row->f1r2_A_Q;
  this->fout << "\t";
  this->fout << row->f1r2_C_Q;
  this->fout << "\t";
  this->fout << row->f1r2_G_Q;
  this->fout << "\t";
  this->fout << row->f1r2_T_Q;
  this->fout << "\t";
  this->fout << row->f2r1_A_Q;
  this->fout << "\t";
  this->fout << row->f2r1_C_Q;
  this->fout << "\t";
  this->fout << row->f2r1_G_Q;
  this->fout << "\t";
  this->fout << row->f2r1_T_Q;
  this->fout << "\t";
  this->fout << row->bfwd_total;
  this->fout << "\t";
  this->fout << row->brev_total;
  this->fout << "\t";
  this->fout << row->f1r2_total;
  this->fout << "\t";
  this->fout << row->f2r1_total;
  this->fout << "\t";
  this->fout << row->left;
  this->fout << "\t";
  this->fout << row->right;
  this->fout << "\t";
  this->fout << row->min_qpos;
  this->fout << "\t";
  this->fout << VariantCaller::StrandMismatch(row);;
  this->fout << "\t";
  this->fout << row->ismasked;
  this->fout << "\t";
  this->fout << row->dplx_barcode;
  this->fout << std::endl;
}


void VariantCaller::WriteMetrics() {
  // call versus query position
  for (auto c1 : this->call_by_qpos) {
    for (auto c2 : c1.second) {
      for (auto c3 : c2.second) {
        this->fout << "CallVsQpos";
        this->fout << "\t";
        this->fout << c1.first;
        this->fout << "\t";
        this->fout << c2.first;
        this->fout << "\t";
        this->fout << c3.first;
        this->fout << "\t";
        this->fout << c3.second;
        this->fout << std::endl;
      }
    }
  }
  // pyrimidine context versus mask
  for (auto p1 : this->pyr_by_mask) {
    for (auto p2 : p1.second) {
      this->fout << "PyrVsMask";
      this->fout << "\t";
      this->fout << p1.first;
      this->fout << "\t";
      this->fout << p2.first;
      this->fout << "\t";
      this->fout << p2.second;
      this->fout << std::endl;
    }
  }
  // read bundle sizes
  for (auto r1 : this->read_bundles) {
    for (auto r2 : r1.second) {
      for (auto r3 : r2.second) {
        for (auto r4 : r3.second) {
          this->fout << "ReadBundles";
          this->fout << "\t";
          this->fout << r1.first;
          this->fout << "\t";
          this->fout << r2.first;
          this->fout << "\t";
          this->fout << r3.first;
          this->fout << "\t";
          this->fout << r4.first;
          this->fout << "\t";
          this->fout << r4.second;
          this->fout << std::endl;
        }
      }
    }
  }
  // mutation burdens
  for (auto b1 : this->burdens) {
    for (auto b2 : b1.second) {
      this->fout << "Burdens";
      this->fout << "\t";
      this->fout << b1.first;
      this->fout << "\t";
      this->fout << b2.first;
      this->fout << "\t";
      this->fout << b2.second;
      this->fout << std::endl;
    }
  }
  // coverage
  this->fout << "Coverage";
  this->fout << "\t";
  this->fout << this->coverage;
  this->fout << std::endl;
  this->fout.close();
}


void Usage() {
  fprintf(stderr, "\nUsage:\n");
  fprintf(stderr, "\t-B\tBED input file name\n");
  fprintf(stderr, "\t-a\tminimum AS-XS\n");
  fprintf(stderr, "\t-b\tminimum number of bulk reads per strand\n");
  fprintf(stderr, "\t-z\tminimum number of bulk reads in total\n");
  fprintf(stderr, "\t-c\tmaximum fraction of clips\n");
  fprintf(stderr, "\t-d\tminimum number of dplx reads per strand\n");
  fprintf(stderr, "\t-f\tminimum fraction of reads for consensus\n");
  fprintf(stderr, "\t-i\tmaximum fraction of reads with an indel\n");
  fprintf(stderr, "\t-n\tmaximun number of mismatches\n");
  fprintf(stderr, "\t-p\tminimum fraction of reads that are proper-pairs\n");
  fprintf(stderr, "\t-q\tminimum consensus base quality\n");
  fprintf(stderr, "\t-r\tread length (after 5' trimming)\n");
  fprintf(stderr, "\t-m\tminimum cycle number\n");
  fprintf(stderr, "\t-x\tmaximum cycle number\n");
  fprintf(stderr, "\t-v\tmaximum bulk VAF\n");
  fprintf(stderr, "\t-O\tprefix of the output file\n");
  fprintf(stderr, "\t-D\tprefix of the discarded variant file [optional]\n");
  fprintf(stderr, "\t-U\tcoverage output file [optional]\n");
  fprintf(stderr, "\t-h\tHelp\n\n");
}


void Options(int argc, char **argv, VariantCaller *vc) {
  vc->bed               = NULL;          // permissive -> strict
  vc->asxs              = 50;
  vc->bulk              = 5;
  vc->bulk_total        = 10;
  vc->clip              = 0.1;           // 1 to 0
  vc->dplx              = 2;
  vc->frac              = 0.9;           // 0 to 1
  vc->indel             = 1;             // 1 to 0
  vc->nmms              = 3;
  vc->ppair             = 0.95;          // 0 to 1
  vc->qual              = 45;
  vc->readlen           = 146;
  vc->min_cycle         = 10;
  vc->max_cycle         = 10;
  vc->vaf               = 0.01;          // 1 to 0
  vc->outfile           = NULL;
  vc->outfile_discarded = NULL;
  int opt = 0;
  char suffix[] = ".gz";
  char buffer[600];
  while ((opt = getopt(argc, argv, "B:a:b:z:c:d:f:i:m:n:p:q:r:v:x:O:D:U:h")) >= 0) {
    switch (opt) {
      case 'B':
        vc->bed = optarg;
        break;
      case 'a':
        vc->asxs = std::stoi(optarg);
        break;
      case 'b':
        vc->bulk = std::stoi(optarg);
        break;
      //fa8:
      case 'z':
        vc->bulk_total = std::stoi(optarg);
        break;
      //
      case 'c':
        vc->clip = std::stof(optarg);
        break;
      case 'd':
        vc->dplx = std::stoi(optarg);
        break;
      case 'f':
        vc->frac = std::stof(optarg);
        break;
      case 'i':
        vc->indel = std::stof(optarg);
        break;
      case 'm':
        vc->min_cycle = std::stoi(optarg);
        break;
      case 'n':
        vc->nmms = std::stoi(optarg);
        break;
      case 'p':
        vc->ppair = std::stof(optarg);
        break;
      case 'q':
        vc->qual = std::stoi(optarg);
        break;
      case 'r':
        vc->readlen = std::stoi(optarg);
        break;
      case 'v':
        vc->vaf = std::stof(optarg);
        break;
      case 'x':
        vc->max_cycle = std::stoi(optarg);
        break;
      case 'O':
        vc->outfile = optarg;
        break;
      case 'D':
        vc->outfile_discarded = optarg;
        break;
      case 'U': //coverage
        strcpy(buffer,optarg);
        strcat(buffer,suffix);
        vc->outfile_coverage = buffer;
        break;
      case 'h':
        Usage();
        exit(0);
      default:
        break;
    }
  }
  if (vc->frac <= 0.5) {
    std::stringstream er;
    er << "Error: minimum fraction of reads to form consensus must be > 0.5";
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
  vc->fout.open(vc->outfile);
  if ( ! vc->fout.is_open() ) {
    std::stringstream er;
    er << "Error: cannot write output file ";
    er << vc->outfile;
    er << std::endl;
    throw std::runtime_error(er.str());
  }
  if (vc->outfile_discarded != NULL) { 
    vc->fout_discarded.open(vc->outfile_discarded);
    if ( ! vc->fout_discarded.is_open() ) {
      std::stringstream er;
      er << "Error: cannot write discarded_output file ";
      er << vc->outfile_discarded;
      er << std::endl;
      throw std::runtime_error(er.str());
    }
  }

  // Record the command
  int i;
  vc->fout <<"# ";
  for (i=0; i < argc; i++)
    {
    vc->fout << argv[i];
    vc->fout << " ";
    }
  vc->fout << std::endl;

  if (vc->outfile_discarded != NULL) { // If discarded_outfile exists, record command there too
    vc->fout_discarded << "# ";
    for (i=0; i < argc; i++) {
      vc->fout_discarded << argv[i];
      vc->fout_discarded << " ";
    }
    vc->fout_discarded << std::endl;
  }

  if(vc->outfile_coverage != NULL) {
    vc->gzout_coverage.open(vc->outfile_coverage);
  }
}


int main(int argc, char **argv) {
  VariantCaller vc;
  Options(argc, argv, &vc);
  vc.CollectMetrics();
  vc.WriteMetrics();
  return 0;
}
