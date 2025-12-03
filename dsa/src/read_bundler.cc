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


#include "read_bundler.h"


const int ALPH_LEN     = 4;
const char ALPH[4]     = {'A', 'C', 'G', 'T'};
const double ALT_BASES = static_cast<double>(ALPH_LEN - 1);
const double POWER     = static_cast<double>(10);

std::map<int, std::map<int, int>> RTYPES =
  { {0, {{1, 0}, {2, 1}}}, {1, {{2, 0}, {1, 1}}} };


int ReadBundler::AuxTagIsPresent(bam1_t* b, const char* tag) {
  return (bam_aux_get(b, tag))? 1 : 0;
}


char* ReadBundler::AuxTagToChar(bam1_t* b, const char* tag) {
  uint8_t *ptr = bam_aux_get(b, tag);
  if (!ptr) {
    std::stringstream er;
    er << "Error: read ";
    er << bam_get_qname(b);
    er << " does not have a ";
    er << tag;
    er << " tag";
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
  return bam_aux2Z(ptr);
}


int ReadBundler::AuxTagToInt(bam1_t* b, const char* tag) {
  uint8_t *ptr = bam_aux_get(b, tag);
  //fa8:
  //std::cerr << tag << " : " << ptr << " : " << bam_aux2i(ptr) << " : " << bam_get_qname(b) << std::endl;
  if (!ptr) {
    std::stringstream er;
    er << "Error: read ";
    er << bam_get_qname(b);
    er << " does not have a ";
    er << tag;
    er << " tag";
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
  return bam_aux2i(ptr);
}


int ReadBundler::ASMinusXS(bam1_t* b) {
  //fa8:
  //std::cerr << ReadBundler::AuxTagToInt(b, "AS") << "/" << ReadBundler::AuxTagToInt(b, "XS") << " : " << bam_get_qname(b) << std::endl;
  return ReadBundler::AuxTagToInt(b, "AS") - ReadBundler::AuxTagToInt(b, "XS");
}


int ReadBundler::IsFivePrimeClipped(bam1_t* b, int strand) {
  uint32_t *cigar = bam_get_cigar(b);
  if ((strand == 0) && (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP)) {
    return 1;
  } else if ((strand == 1) &&
      (bam_cigar_op(cigar[b->core.n_cigar-1]) == BAM_CSOFT_CLIP)) {
    return 1;
  } else {
    return 0;
  }
}


// Returns false if position is not within template; for example, if the read
// extends past its mate into adaptor
bool ReadBundler::IsTemplate(int beg, int end) {
  if ((this->pos >= (beg - this->offset)) &&
      (this->pos <= (end - this->offset))) {
    return true;
  } else {
    return false;
  }
}


int ReadBundler::ReadStrand(bam1_t* b) {
  if (b->core.flag & BAM_FMREVERSE) {
    return 0;
  } else if (b->core.flag & BAM_FREVERSE) {
    return 1;
  } else {
    // read is not mapped to a strand, if proper pairs are required
    // this can be set to throw std::invalid_argument("Invalid strand");
    return -1;
  }
}


int ReadBundler::ReadNumber(bam1_t* b) {
  if (b->core.flag & BAM_FREAD1) {
    return 1;
  } else if (b->core.flag & BAM_FREAD2) {
    return 2;
  } else {
    throw std::invalid_argument("Invalid read number");
  }
}


std::pair<char, int> ReadBundler::BaseAndQual(const bam_pileup1_t* p) {
  // We don't make use of indel quality scores, so gives these value -1
  if ((p->is_del) || (p->indel != 0)) {
    return std::make_pair('i', -1);
  } else {
    uint8_t *seq = bam_get_seq(p->b);
    char base = static_cast<char>(seq_nt16_str[bam_seqi(seq, p->qpos)]);
    int qual  = bam_get_qual(p->b)[p->qpos];
    return std::make_pair(base, qual);
  }
}


// No pre-processing is performed on bulk bam
// Keep non-properly-paired reads
bool ReadBundler::BulkIsUsable(bam1_t *b) {
  int unmapped  = (b->core.flag & BAM_FUNMAP)? 1 : 0;
  int secondary = (b->core.flag & BAM_FSECONDARY)? 1 : 0;
  int qcfail    = (b->core.flag & BAM_FQCFAIL)? 1 : 0;
  int suppl     = (b->core.flag & BAM_FSUPPLEMENTARY)? 1 : 0;
  int duplicate = (b->core.flag & BAM_FDUP)? 1 : 0;
  if ((unmapped + secondary + qcfail + suppl + duplicate) == 0) {
    return true;
  }
  return false;
}


int ReadBundler::IsProperPair(bam1_t *b) {
  return (b->core.flag & BAM_FPROPER_PAIR)? 1 : 0;
}


identifier ReadBundler::DplxIdentifier(const bam_pileup1_t* p) {
  std::string idf1 = ReadBundler::AuxTagToChar(p->b, "RB");
  std::istringstream iss(idf1);
  std::vector<std::string> tokens;
  std::string token;
  while (std::getline(iss, token, ',')) {
    tokens.push_back(token);
  }
  identifier idf;
  idf.id     = idf1;
  idf.beg    = std::stoi(tokens[1]);
  idf.end    = std::stoi(tokens[2]);
  idf.fwd_bc = tokens[3];
  idf.rev_bc = tokens[4];
  return idf;
}


void ReadBundler::UpdateDplxBundle(identifier idf, bundle* bndl,
  const bam_pileup1_t* p) {
  bndl->idf = idf;
  int strand = ReadBundler::ReadStrand(p->b);
  int read   = ReadBundler::ReadNumber(p->b);
  int rtype  = RTYPES[strand][read];
  std::pair<char, int> bq = ReadBundler::BaseAndQual(p);
  bndl->dplx_depth[strand][read]++;
  bndl->counts[rtype][bq.first]++;
  bndl->call[rtype].push_back(bq);
  //fa8
  //std::cerr << "bundle" << std::endl;
  bndl->asxs[rtype].push_back(ReadBundler::ASMinusXS(p->b));
  bndl->clip[rtype].push_back(ReadBundler::IsFivePrimeClipped(p->b, strand));
  bndl->nmms[rtype].push_back(ReadBundler::AuxTagToInt(p->b, "NM"));
  bndl->ppair[rtype].push_back(ReadBundler::IsProperPair(p->b));
}


void ReadBundler::UpdateBulkBundle(bundle* bndl, const bam_pileup1_t* p,
  int min_base_quality) {
  int strand = ReadBundler::ReadStrand(p->b);
  std::pair<char, int> bq = ReadBundler::BaseAndQual(p);
  // only use bulk bundles where base quality is >= threshold
  // fa8:
  //std::cerr << bq.first << ":" << bq.second << "(" << min_base_quality << ")" << std::endl;
  // rob's
  //if (bq.second >= min_base_quality) {
  // fa8 (disabling the filter for indels because they may have -1)
  if (bq.first == 'i' || bq.second >= min_base_quality) {
    // fa8:
  	//std::cerr << "   passed" << std::endl;
    bndl->counts[strand][bq.first]++;
    bndl->call[strand].push_back(bq);
    //fa8
    //std::cerr << "bulk" << std::endl;
    bndl->asxs[strand].push_back(ReadBundler::ASMinusXS(p->b));
    bndl->nmms[strand].push_back(ReadBundler::AuxTagToInt(p->b, "NM"));
    bndl->ppair[strand].push_back(ReadBundler::IsProperPair(p->b));
  }
}


void ReadBundler::DplxConsensus(bundle* bndl) {
  for (int i = 0; i < 2; i++) {
    // sum log10 probability of error
    std::vector<double> probs(ALPH_LEN, static_cast<double>(0));
    for (int j = 0; j < bndl->call[i].size(); j++) {
      char base = bndl->call[i][j].first;
      int qual  = bndl->call[i][j].second;
      // base is canonical
      if (memchr(ALPH, base, sizeof(ALPH))) {
        double perror   = std::pow(POWER, (-qual/POWER));
        double pcorrect = (static_cast<double>(1) - perror)/ALT_BASES;
        for (int k = 0; k < ALPH_LEN; k++) {
          if (base == ALPH[k]) {
            probs[k] += std::log10(perror);
          } else {
            probs[k] += std::log10(pcorrect);
          }
        }
      }
    }
    // log10 sum exp
    double maxp   = *std::max_element(probs.cbegin(), probs.cend());
    double sumexp = 0;
    for (int j = 0; j < probs.size(); j++) {
      sumexp += std::pow(POWER, probs[j] - maxp);
    }
    double logsumexp = std::log10(sumexp) + maxp;
    // normalize sum log10 probs and convert to Phred based quality score
    for (int j = 0; j < probs.size(); j++) {
      bndl->consensus[i].push_back((probs[j] - logsumexp) * -POWER);
    }
  }
}


bundles ReadBundler::DplxBundles(int pos, int offset, int min_dplx_depth,
  pileups plps) {
  this->pos    = pos;
  this->offset = offset;
  bundles bouts;
  for (int i = 0; i < plps.size(); i++) {
    const bam_pileup1_t* p = plps[i];
    identifier idf = ReadBundler::DplxIdentifier(p);
    if (ReadBundler::IsTemplate(idf.beg, idf.end) == 1) {
      ReadBundler::UpdateDplxBundle(idf, &bouts[idf.id], p);
    }
  }
  for (auto it = bouts.begin(); it != bouts.end(); ) {
    // define duplex type: 0 = no duplex, 1 = fwd, 2 = rev, 3 = fwd and rev
    int bundle_type = 0;
    if ((it->second.dplx_depth[0][1] >= min_dplx_depth) &&
        (it->second.dplx_depth[0][2] >= min_dplx_depth)) {
      bundle_type += 1;
    }
    if ((it->second.dplx_depth[1][1] >= min_dplx_depth) &&
        (it->second.dplx_depth[1][2] >= min_dplx_depth)) {
      bundle_type += 2;
    }
    // remove low depth bundles
    if (bundle_type == 0) {
      bouts.erase(it++);
    } else {
      // calculate consensus base qualities
      ReadBundler::DplxConsensus(&bouts[it->second.idf.id]);
      bouts[it->second.idf.id].bundle_type = bundle_type;
      it++;
    }
  }
  return bouts;
}


bundle ReadBundler::BulkBundle(pileups plps, int min_base_quality) {
  bundle bndl;
  for (int i =0; i < plps.size(); i++) {
    const bam_pileup1_t* p = plps[i];
    if (ReadBundler::BulkIsUsable(p->b)) {
      ReadBundler::UpdateBulkBundle(&bndl, p, min_base_quality);
    }
  }
  return bndl;
}
