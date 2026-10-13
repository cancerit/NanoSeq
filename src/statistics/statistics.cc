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

Run:
  ./statistics -O original.bam -F filtered.bam -M 2

*/


#include "statistics.h"


std::vector<std::string> RNAMES =
  { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
    "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" };


int Statistics::ReadStrand(bam1_t* b) {
  if (b->core.flag & BAM_FMREVERSE) {
    return 0;
  } else if (b->core.flag & BAM_FREVERSE) {
    return 1;
  } else {
    throw std::invalid_argument("Invalid strand");
  }
}


int Statistics::ReadNumber(bam1_t* b) {
  if (b->core.flag & BAM_FREAD1) {
    return 1;
  } else if (b->core.flag & BAM_FREAD2) {
    return 2;
  } else {
    throw std::invalid_argument("Invalid read number");
  }
}


std::vector<std::string> Statistics::Tokenize(std::string str,
  char delimiter) {
  std::istringstream iss(str);
  std::vector<std::string> tokens;
  std::string token;
  while (std::getline(iss, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}


bool Statistics::BamIsPreProcessed(bam_hdr_t *head) {
  int ops = 0;
  std::vector<std::string> tokens = Tokenize(head->text, '\n');
  for (int j = 0; j < tokens.size(); j++) {
    if (tokens[j].rfind("@PG", 0) == 0) {
      std::vector<std::string> subtokens = Tokenize(tokens[j], '\t');
      for (int k = 0; k < subtokens.size(); k++) {
        if (subtokens[k].rfind("ID", 0) == 0) {
          // account for suffixes on program identifiers
          if (subtokens[k].rfind("ID:bamsormadup", 0) == 0) {
            ops += 1;
          } else if (subtokens[k].rfind("ID:bammarkduplicatesopt", 0) == 0) {
            ops += 1;
          } else if (subtokens[k].rfind("ID:bamaddreadbundles", 0) == 0) {
            ops += 1;
          }
        }
      }
    }
  }
  return(1);
  //return (ops == 3);
}


void Statistics::LoadFiles() {
  this->in1 = hts_open(this->input_bam, "r");
  if (this->in1 == 0) {
    std::stringstream er;
    er << "Fail to open input BAM/CRAM file ";
    er << this->input_bam;
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
  this->in2 = hts_open(this->filtered_bam, "r");
  if (this->in2 == 0) {
    std::stringstream er;
    er << "Fail to open input BAM/CRAM file ";
    er << this->filtered_bam;
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
  this->idx1 = sam_index_load(this->in1, this->input_bam);
  if (this->idx1 == NULL) {
    std::stringstream er;
    er << "Fail to open BAM index for ";
    er << this->input_bam;
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
  this->idx2 = sam_index_load(this->in2, this->filtered_bam);
  if (this->idx2 == NULL) {
    std::stringstream er;
    er << "Fail to open BAM index for ";
    er << this->filtered_bam;
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
  this->head1 = sam_hdr_read(this->in1);
  if ((this->head1 == NULL) || (this->head1->n_targets == 0)) {
    std::stringstream er;
    er << "Error: BAM/CRAM file ";
    er << this->input_bam;
    er << " does not have header.";
    er << std::endl;
    throw std::runtime_error(er.str());
  }
  if (BamIsPreProcessed(this->head1) == false) {
    std::stringstream er;
    er << this->input_bam;
    er << " is not correctly preprocessed. Input is a BAM file that has been ";
    er << "processed using biobambam2 bamsormadup using rcsupport=1, ";
    er << "bammarkduplicatesopt and bamaddreadbundles.";
  }
  this->head2 = sam_hdr_read(this->in2);
  if ((this->head2 == NULL) || (this->head2->n_targets == 0)) {
    std::stringstream er;
    er << "Error: BAM/CRAM file ";
    er << this->filtered_bam;
    er << " does not have header.";
    er << std::endl;
    throw std::runtime_error(er.str());
  }
  if (BamIsPreProcessed(this->head2) == false) {
    std::stringstream er;
    er << this->filtered_bam;
    er << " is not correctly preprocessed. Input is a BAM file that has been ";
    er << "processed using biobambam2 bamsormadup using rcsupport=1, ";
    er << "bammarkduplicatesopt and bamaddreadbundles.";
    er << std::endl;
    throw std::runtime_error(er.str());
  }
}


float Statistics::MeanOfOneVector(std::vector<int> v1) {
  float total = std::accumulate(v1.begin(), v1.end(), 0);
  return std::round(total/static_cast<float>(v1.size()));
}


void Statistics::CollectMetrics() {
  this->reads = 0;
  for (int i = 0; i < RNAMES.size(); i++) {
    char rname[RNAMES[i].length() + 1];
    strcpy(rname, RNAMES[i].c_str());
    char region[65536];
    sprintf(region, "%s:%d-%d", rname, 0, INT_MAX);
    std::cout << RNAMES[i] << std::endl;
    // count reads in original bam
    hts_itr_t *iter1 = sam_itr_querys(this->idx1, this->head1, region);
    bam1_t *b1       = bam_init1();
    while (sam_itr_next(this->in1, iter1, b1) >= 0) {
      this->reads++;
    }
    // count bundles in filtered bam
    hts_itr_t *iter2 = sam_itr_querys(this->idx2, this->head2, region);
    bam1_t *b2       = bam_init1();
    std::map<std::string, std::map<int, std::map<int, int>>> counts;
    std::map<std::string, std::vector<int>> tlens;
    while (sam_itr_next(this->in2, iter2, b2) >= 0) {
      std::string rb = std::string(bam_aux2Z(bam_aux_get(b2, "RB")));
      int strand     = Statistics::ReadStrand(b2);
      int read       = Statistics::ReadNumber(b2);
      counts[rb][strand][read]++;
      tlens[rb].push_back(std::abs(b2->core.isize));
    }
    for (auto count : counts) {
      int f1 = count.second[0][1];
      int f2 = count.second[0][2];
      int r1 = count.second[1][1];
      int r2 = count.second[1][2];
      this->reads_per_bundle[(f1+f2+r1+r2)]++;
      if (((f1 >= this->min_dplx_depth) && (f2 >= this->min_dplx_depth)) ||
          ((r1 >= this->min_dplx_depth) && (r2 >= this->min_dplx_depth))) {
        this->reads_per_usable_bundle[(f1+f2+r1+r2)]++;
        this->tlens_per_usable_bundle[Statistics::MeanOfOneVector(tlens[count.first])]++;
      }
    }
    counts.clear();
    tlens.clear();
    bam_destroy1(b1);
    bam_destroy1(b2);
    hts_itr_destroy(iter1);
    hts_itr_destroy(iter2);
  }
  bam_hdr_destroy(this->head1);
  bam_hdr_destroy(this->head2);
  hts_close(this->in1);
  hts_close(this->in2);
}


void Statistics::WriteOut() {
  this->bundle_metrics.open(this->dir + "ReadBundleMetrics.bin");
  this->bundle_metrics << "paired-end reads per read-bundle\tcount\n";
  int bundles = 0;
  int reads_in_bundles = 0;
  for (auto rpb : this->reads_per_bundle) {
    bundles += rpb.second;
    reads_in_bundles += (rpb.first * rpb.second);
  }
  int usable_bundles = 0;
  float two = 2;
  float reads_in_usable_bundles = 0;
  for (auto rpb : this->reads_per_usable_bundle) {
    usable_bundles += rpb.second;
    reads_in_usable_bundles += (rpb.first * rpb.second);
    this->bundle_metrics << rpb.first/two;
    this->bundle_metrics << "\t";
    this->bundle_metrics << rpb.second;
    this->bundle_metrics << "\n";
  }
  this->tlen_metrics.open(this->dir + "TemplateLengthMetrics.bin");
  this->tlen_metrics << "template-length\tusable read-bundles\n";
  for (auto tlen : this->tlens_per_usable_bundle) {
    this->tlen_metrics << tlen.first;
    this->tlen_metrics << "\t";
    this->tlen_metrics << tlen.second;
    this->tlen_metrics << "\n";
  }
  this->filter_metrics.open(this->dir + "FilteringMetrics.bin");
  this->filter_metrics << "pre-filtered paired-end reads";
  this->filter_metrics << "\t";
  this->filter_metrics << reads/two;
  this->filter_metrics << "\n";
  this->filter_metrics << "read-bundles";
  this->filter_metrics << "\t";
  this->filter_metrics << bundles;
  this->filter_metrics << "\n";
  this->filter_metrics << "paired-end reads in read-bundles";
  this->filter_metrics << "\t";
  this->filter_metrics << reads_in_bundles/two;
  this->filter_metrics << "\n";
  this->filter_metrics << "usable read-bundles";
  this->filter_metrics << "\t";
  this->filter_metrics << usable_bundles;
  this->filter_metrics << "\n";
  this->filter_metrics << "paired-end reads in usable read-bundles";
  this->filter_metrics << "\t";
  this->filter_metrics << reads_in_usable_bundles/two;
  this->filter_metrics << "\n";
  this->filter_metrics << "duplicate rate";
  this->filter_metrics << "\t";
  this->filter_metrics << 1 - (bundles/(reads_in_bundles/two));
  this->filter_metrics << "\n";
  this->filter_metrics.close();
  this->bundle_metrics.close();
}


void Usage() {
  fprintf(stderr, "\nUsage:\n");
  fprintf(stderr, "\t-I\tInput BAM file name\n");
  fprintf(stderr, "\t-F\tFiltered BAM file name\n");
  fprintf(stderr, "\t-O\tOutput directory name\n");
  fprintf(stderr, "\t-M\tMinimum duplex depth (default = 2)\n");
  fprintf(stderr, "\t-h\tHelp\n");
}


int main(int argc, char **argv) {
  Statistics stats;
  stats.min_dplx_depth = 2;
  int opt = 0;
  while ((opt = getopt(argc, argv, "I:F:M:O:h")) >= 0) {
    switch (opt) {
      case 'I':
        stats.input_bam = optarg;
        break;
      case 'F':
        stats.filtered_bam = optarg;
        break;
      case 'M':
        stats.min_dplx_depth = std::stoi(optarg);
        break;
      case 'O':
        stats.dir = optarg;
        break;
      case 'h':
        Usage();
        exit(0);
      default:
        break;
    }
  }
  stats.LoadFiles();
  stats.CollectMetrics();
  stats.WriteOut();
  return 0;
}
