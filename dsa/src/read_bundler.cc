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

#include "read_bundler.h"


/*
std::map<int, std::map<int, int>> RTYPES =
  { {0, {{1, 0}, {2, 1}}}, {1, {{2, 0}, {1, 1}}} };
*/

static const int RTYPES[STRAND_COUNT][READ_TYPE_COUNT] = {
  {0, 1},
  {1, 0}
};

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
bool ReadBundler::IsTemplate(const int beg, const int end) {
  return
    (this->pos >= (beg - this->offset)) &&
    (this->pos <= (end - this->offset));
}

static inline int get_strand_index(const bam1_t *b) {
  // ASSUMPTION: proper pair and strand have already been verified
  static_assert(STRAND_INDEX_FORWARD == 0);
  static_assert(STRAND_INDEX_REVERSE == 1);

  // TODO: optimise
  if (read_has_flag(b, BAM_FMREVERSE)) {
    return STRAND_INDEX_FORWARD;
  } else if (read_has_flag(b, BAM_FREVERSE)) {
    return STRAND_INDEX_REVERSE;
  } else {
    return STRAND_INDEX_IGNORE;
  }
}

static inline int get_read_type_index(const bam1_t *b) {
  // ASSUMPTION: that one and only one BAM_FREAD* flag is set to be verified earlier
  //  (this would be an unpredictable branch, while flag consistency verification is predictable)
  static_assert(READ_TYPE_INDEX_READ_1 == 0);
  static_assert(READ_TYPE_INDEX_READ_2 == 1);
  return read_has_flag(b, BAM_FREAD2);
}

/*
int ReadBundler::ReadStrand(bam1_t* b) {
  // TODO: consider the implications...
  if (b->core.flag & BAM_FMREVERSE) {
    return STRAND_INDEX_FORWARD;
  } else if (b->core.flag & BAM_FREVERSE) {
    return STRAND_INDEX_REVERSE;
  } else {
    // read is not mapped to a strand, if proper pairs are required
    // this can be set to throw std::invalid_argument("Invalid strand");
    return -1;
  }
}

int ReadBundler::ReadNumber(bam1_t* b) {
  // TODO: optimise!
  if (b->core.flag & BAM_FREAD1) {
    return READ_TYPE_INDEX_READ_1;
  } else if (b->core.flag & BAM_FREAD2) {
    return READ_TYPE_INDEX_READ_2;
  } else {
    throw std::invalid_argument("Invalid read number");
  }
}
*/

const int nt16_allele[16] = {
  [0] = ALLELE_DISCARDED,
  [BAM_NT_A] = ALLELE_A,
  [BAM_NT_C] = ALLELE_C,
  [3] = ALLELE_DISCARDED,
  [BAM_NT_G] = ALLELE_G,
  [5] = ALLELE_DISCARDED,
  [6] = ALLELE_DISCARDED,
  [7] = ALLELE_DISCARDED,
  [BAM_NT_T] = ALLELE_T,
  [9] = ALLELE_DISCARDED,
  [10] = ALLELE_DISCARDED,
  [11] = ALLELE_DISCARDED,
  [12] = ALLELE_DISCARDED,
  [13] = ALLELE_DISCARDED,
  [14] = ALLELE_DISCARDED,
  [15] = ALLELE_DISCARDED,
};

std::pair<int, int> ReadBundler::BaseAndQual(const bam_pileup1_t *p) {
  // We don't make use of indel quality scores, so gives these value -1
  if ((p->is_del) || (p->indel != 0)) {
    return std::make_pair(ALLELE_DEL, -1);
  } else {
    uint8_t *seq = bam_get_seq(p->b);
    int base = nt16_allele[bam_seqi(seq, p->qpos)];
    int qual  = bam_get_qual(p->b)[p->qpos];
    return std::make_pair(base, qual);
  }
}

#define BULK_UNUSABLE (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FSUPPLEMENTARY | BAM_FDUP)

// No pre-processing is performed on bulk bam
// Keep non-properly-paired reads
bool ReadBundler::BulkIsUsable(bam1_t *b) {
  return !read_has_flag(b, BULK_UNUSABLE);
}

static inline duplex_tag_info parse_identifier(std::string idf1) {
  std::istringstream iss(idf1);
  std::vector<std::string> tokens;
  std::string token;

  // TODO: optimise!
  while (std::getline(iss, token, ',')) {
    tokens.push_back(token);
  }
  duplex_tag_info idf = {
    .beg = std::stoi(tokens[1]),
    .end = std::stoi(tokens[2]),
    .fwd_bc = tokens[3],
    .rev_bc = tokens[4]
  };
  return idf;
}

std::string ReadBundler::DplxIdentifier(const bam_pileup1_t* p) {
  return ReadBundler::AuxTagToChar(p->b, "RB");
}

void ReadBundler::UpdateDplxBundle(bundle *bndl, const bam_pileup1_t *p) {
  const int strand = get_strand_index(p->b);
  const int read = get_read_type_index(p->b);
  int rtype;
  std::pair<int, int> bq = ReadBundler::BaseAndQual(p);
  if (strand != STRAND_INDEX_IGNORE) {
    bndl->dplx_depth[strand][read]++;
    rtype = RTYPES[strand][read];
  } else {
    // NOTE: should replicate missing key (strand) in the original implementation
    rtype = 0;
  }
  bndl->counts[rtype][bq.first]++;
  bndl->call[rtype].push_back(bq);
  bndl->asxs[rtype].push_back(ReadBundler::ASMinusXS(p->b));

  // TODO: verify STRAND_INDEX_IGNORE is supported correctly
  bndl->clip[rtype].push_back(ReadBundler::IsFivePrimeClipped(p->b, strand));

  bndl->nmms[rtype].push_back(ReadBundler::AuxTagToInt(p->b, "NM"));
  // bndl->ppair[rtype].push_back(read_is_in_proper_pair(p->b));
  bndl->rtype_ppair_counts[rtype] += read_is_in_proper_pair(p->b);
  bndl->rtype_read_counts[rtype]++;
}

void ReadBundler::UpdateBulkBundle(bundle *bndl, const bam_pileup1_t *p, int min_base_quality) {
  // TODO: verify behaviour on invalid strand!
  const int strand = get_strand_index(p->b);
  std::pair<int, int> bq = ReadBundler::BaseAndQual(p);
  // only use bulk bundles where base quality is >= threshold
  // fa8:
  //std::cerr << bq.first << ":" << bq.second << "(" << min_base_quality << ")" << std::endl;
  // rob's
  //if (bq.second >= min_base_quality) {
  // fa8 (disabling the filter for indels because they may have -1)
  if (bq.first == ALLELE_DEL || bq.second >= min_base_quality) {
    bndl->counts[strand][bq.first]++;
    bndl->call[strand].push_back(bq);
    bndl->asxs[strand].push_back(ReadBundler::ASMinusXS(p->b));
    bndl->nmms[strand].push_back(ReadBundler::AuxTagToInt(p->b, "NM"));
    // bndl->ppair[strand].push_back(read_is_in_proper_pair(p->b));
    bndl->rtype_ppair_counts[strand] += read_is_in_proper_pair(p->b);
    bndl->rtype_read_counts[strand]++;
  }
}

// TODO: optimise calculation!
void ReadBundler::DplxConsensus(bundle *bndl) {
  for (int i = 0; i < 2; i++) {
    std::vector<double> probs(ALPH_LEN, static_cast<double>(0));
    assert(probs.size() == ALPH_LEN);
    // sum log10 probability of error
    for (int j = 0; j < bndl->call[i].size(); j++) {
      int base = bndl->call[i][j].first;
      int qual  = bndl->call[i][j].second;
      // base is canonical
      if (base != ALLELE_DISCARDED && base != ALLELE_DEL) {
        int base_index = base - 1;
        double perror   = std::pow(POWER, (-qual/POWER));
        double pcorrect = (static_cast<double>(1) - perror) / ALT_BASES;
        for (int k = 0; k < ALPH_LEN; k++) {
          if (base_index == k) {
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

static inline int high_duplex_depth(const bundle *b, const int strand, const int min_dplx_depth) {
  return static_cast<int>(
    (b->dplx_depth[strand][READ_TYPE_INDEX_READ_1] >= min_dplx_depth) &&
    (b->dplx_depth[strand][READ_TYPE_INDEX_READ_2] >= min_dplx_depth));
}

static inline int get_duplex_bundle_type(const bundle *b, const int min_dplx_depth) {
  return
    (high_duplex_depth(b, STRAND_INDEX_REVERSE, min_dplx_depth) << 1) |
    (high_duplex_depth(b, STRAND_INDEX_FORWARD, min_dplx_depth) << 0);
}

bundles ReadBundler::DplxBundles(int pos, int offset, int min_dplx_depth, pileups plps) {
  this->pos    = pos;
  this->offset = offset;
  bundle *b;
  bundles bouts;
  for (size_t i = 0; i < plps.size(); i++) {
    const bam_pileup1_t *p = plps[i];
    std::string id = ReadBundler::DplxIdentifier(p);
    duplex_tag_info idf = parse_identifier(id);
    if (ReadBundler::IsTemplate(idf.beg, idf.end) == 1) {
      b = &bouts[id];
      b->duplex_tag_info = idf;
      ReadBundler::UpdateDplxBundle(b, p);
    }
  }
  for (auto it = bouts.begin(); it != bouts.end(); ) {
    // define duplex type: 0 = no duplex, 1 = fwd, 2 = rev, 3 = fwd and rev
    int bundle_type = 0;
    if ((it->second.dplx_depth[0][0] >= min_dplx_depth) &&
        (it->second.dplx_depth[0][1] >= min_dplx_depth)) {
      bundle_type += 1;
    }
    if ((it->second.dplx_depth[1][0] >= min_dplx_depth) &&
        (it->second.dplx_depth[1][1] >= min_dplx_depth)) {
      bundle_type += 2;
    }
    // remove low depth bundles
    if (bundle_type == 0) {
      bouts.erase(it++);
    } else {
      // calculate consensus base qualities
      b = &it->second;
      ReadBundler::DplxConsensus(b);
      b->bundle_type = bundle_type;
      it++;
    }
  }
  return bouts;
}

/*
bundles ReadBundler::DplxBundles(int pos, int offset, int min_dplx_depth, pileups plps) {
  this->pos    = pos;
  this->offset = offset;
  bundles bouts;
  bundle *b;

  {
    std::string idf;
    const bam_pileup1_t* p;
    duplex_tag_info info;
    bundles::iterator it;
    for (int i = 0; i < plps.size(); i++) {
      p = plps[i];
      idf = ReadBundler::DplxIdentifier(p);
      it = bouts.find(idf);
      if (it == bouts.end()) {
        info = parse_identifier(idf);
        if (!ReadBundler::IsTemplate(info.beg, info.end)) {
          continue;
        }

        // Initialise
        b = &bouts[idf];
        b->duplex_tag_info = info;
        // b->counts[0][0] = 0;
        // memset(b->counts, 0, 2 * 6 * sizeof(uint64_t));
        ReadBundler::UpdateDplxBundle(b, p);

      } else if (ReadBundler::IsTemplate(bouts[idf].duplex_tag_info.beg, bouts[idf].duplex_tag_info.end)) {

        // Update
        b = &it->second;
        ReadBundler::UpdateDplxBundle(b, p);

      }
    }
  }

  for (auto it = bouts.begin(); it != bouts.end(); ) {
    b = &it->second;
    // define duplex type: 0 = no duplex, 1 = fwd, 2 = rev, 3 = fwd and rev
    const int bundle_type = get_duplex_bundle_type(b, min_dplx_depth);

    // remove low depth bundles
    if (bundle_type == 0) {
      bouts.erase(it++);
    } else {
      // calculate consensus base qualities
      ReadBundler::DplxConsensus(b);
      b->bundle_type = bundle_type;
      it++;
    }
  }
  return bouts;
}
*/

bundle ReadBundler::BulkBundle(pileups plps, int min_base_quality) {
  bundle bndl = {};
  for (size_t i = 0; i < plps.size(); i++) {
    const bam_pileup1_t *p = plps[i];
    if (ReadBundler::BulkIsUsable(p->b)) {
      ReadBundler::UpdateBulkBundle(&bndl, p, min_base_quality);
    }
  }
  return bndl;
}
