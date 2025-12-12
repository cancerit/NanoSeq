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

#include <format>
#include "writeout.h"
#include "utils.h"
#include "compressor.h"
#include "static_string_builder.hpp"

WriteOut::WriteOut(Options *opts) :
  opts(opts),
  compressor(
    opts->oname,
    "dsa.bed.gz",
    opts->compression_level
  ) {}

static inline int32_t get_asxs(const bundle *bin) {
  const float asxs0 = vector_mean(bin->asxs[0]);
  const float asxs1 = vector_mean(bin->asxs[1]);
  const int asxs    = custom_round(std::min(asxs0, asxs1));
  return asxs;
}

static inline int32_t get_clip(const bundle *bin) {
  const float clip1 = vector_pair_mean(bin->clip[0], bin->clip[1]);
  const int clip    = custom_round(clip1);
  return clip;
}

static float get_nmms(const bundle *bin) {
  const float nmms0 = vector_mean(bin->nmms[0]);
  const float nmms1 = vector_mean(bin->nmms[1]);
  const float nmms  = 0.1 * custom_round(std::max(nmms0, nmms1) * 10.0);
  return nmms;
}

static inline int32_t get_consensus(const bundle *bin, const int32_t rtype, const int32_t allele_index) {
  return custom_round(bin->consensus[rtype][allele_index]);
}

static inline std::string dsa_counts_string(const bundle *bin) {
  return std::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
    bin->counts[RTYPE_A][ALLELE_A],
    bin->counts[RTYPE_A][ALLELE_C],
    bin->counts[RTYPE_A][ALLELE_G],
    bin->counts[RTYPE_A][ALLELE_T],
    bin->counts[RTYPE_A][ALLELE_DEL],

    bin->counts[RTYPE_B][ALLELE_A],
    bin->counts[RTYPE_B][ALLELE_C],
    bin->counts[RTYPE_B][ALLELE_G],
    bin->counts[RTYPE_B][ALLELE_T],
    bin->counts[RTYPE_B][ALLELE_DEL]);
}

static inline std::string dsa_base_quals_string(const bundle *bin) {
  return std::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
    custom_round(bin->consensus[RTYPE_A][0]),
    custom_round(bin->consensus[RTYPE_A][1]),
    custom_round(bin->consensus[RTYPE_A][2]),
    custom_round(bin->consensus[RTYPE_A][3]),
    custom_round(bin->consensus[RTYPE_B][0]),
    custom_round(bin->consensus[RTYPE_B][1]),
    custom_round(bin->consensus[RTYPE_B][2]),
    custom_round(bin->consensus[RTYPE_B][3]));
}

static inline float get_ppair_mean(const bundle *bin, const int rtype) {
  // ASSUMPTION: if the proper pair count is not zero, neither will be the read count
  return
    bin->rtype_ppair_counts[rtype] == 0 ? 0.0f :
    ((float)bin->rtype_ppair_counts[rtype] / (float)bin->rtype_read_counts[rtype]);
}

// TODO: pass the same stream to each of these functions?
static inline const float dsa_proper_pair(const bundle *bin) {
  return custom_round(
    std::min(
      get_ppair_mean(bin, RTYPE_A),
      get_ppair_mean(bin, RTYPE_B)));
}

static inline void upper(std::string &str) {
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

static std::string dsa_identifier_string(const bundle *bin) {
  std::string fwd_bc = bin->duplex_tag_info.fwd_bc;
  std::string rev_bc = bin->duplex_tag_info.rev_bc;
  upper(fwd_bc);
  upper(rev_bc);

  return std::format("{}\t{}\t{}|{}\t{}\t",
    bin->duplex_tag_info.beg,
    bin->duplex_tag_info.end,
    fwd_bc, rev_bc,
    bin->bundle_type);
}

void WriteOut::WriteRows(bundle bulk, bundles dplx, std::string posn) {
  // TODO: retry implementing formatting on static buffer (format_to et sim.)

  const std::string prefix = std::format("{}\t{}\t{}\t{}",
    posn,
    get_asxs(&bulk),
    get_nmms(&bulk),
    dsa_counts_string(&bulk));

  std::stringstream b;
  bundle *d;
  for (auto bndls : dplx) {
    d = &bndls.second;

    b
    << prefix
    << dsa_identifier_string(d)
    << std::format("{}\t{}\t{}\t",
      get_asxs(d),
      get_clip(d),
      get_nmms(d))
    << dsa_counts_string(d)
    << dsa_base_quals_string(d)
    << dsa_proper_pair(&bulk) << '\t'
    << dsa_proper_pair(d)     << '\n';

  }

  this->compressor.compress(b.str());
  this->compressor.write();
}

void WriteOut::Finalise() {
  this->compressor.finalise();
}
