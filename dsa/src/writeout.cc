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

#include "writeout.h"
#include "utils.h"

std::ostream &get_output_stream(ogzstream &gzout, const bool out2stdout) {
  return static_cast<std::ostream&>(out2stdout ? gzout : std::cout);
}

WriteOut::WriteOut(Options *opts, ogzstream &gzout, const bool out2stdout) :
  opts(opts),
  out(get_output_stream(gzout, out2stdout)) {}

static inline void push_counts_string(const bundle *bin, std::ostream &ss) {
  ss

  << bin->counts[RTYPE_A][ALLELE_A]   << '\t'
  << bin->counts[RTYPE_A][ALLELE_C]   << '\t'
  << bin->counts[RTYPE_A][ALLELE_G]   << '\t'
  << bin->counts[RTYPE_A][ALLELE_T]   << '\t'
  << bin->counts[RTYPE_A][ALLELE_DEL] << '\t'

  << bin->counts[RTYPE_B][ALLELE_A]   << '\t'
  << bin->counts[RTYPE_B][ALLELE_C]   << '\t'
  << bin->counts[RTYPE_B][ALLELE_G]   << '\t'
  << bin->counts[RTYPE_B][ALLELE_T]   << '\t'
  << bin->counts[RTYPE_B][ALLELE_DEL];
}

static inline void push_base_quals_string(const bundle *bin, std::ostream &ss) {
  ss

  << custom_round(bin->consensus[RTYPE_A][0]) << '\t'
  << custom_round(bin->consensus[RTYPE_A][1]) << '\t'
  << custom_round(bin->consensus[RTYPE_A][2]) << '\t'
  << custom_round(bin->consensus[RTYPE_A][3]) << '\t'

  << custom_round(bin->consensus[RTYPE_B][0]) << '\t'
  << custom_round(bin->consensus[RTYPE_B][1]) << '\t'
  << custom_round(bin->consensus[RTYPE_B][2]) << '\t'
  << custom_round(bin->consensus[RTYPE_B][3]);
}

static void push_duplex_covariates_string(const bundle *bin, std::ostream &ss) {
  const float asxs0 = vector_mean(bin->asxs[0]);
  const float asxs1 = vector_mean(bin->asxs[1]);
  const int asxs    = custom_round(std::min(asxs0, asxs1));

  const float clip1 = vector_pair_mean(bin->clip[0], bin->clip[1]);
  const int clip    = custom_round(clip1);

  const float nmms0 = vector_mean(bin->nmms[0]);
  const float nmms1 = vector_mean(bin->nmms[1]);
  const float nmms  = 0.1 * custom_round(std::max(nmms0, nmms1) * 10.0);

  ss << asxs;
  ss << "\t";
  ss << clip;
  ss << "\t";
  ss << nmms;
}

static void push_bulk_covariates_string(const bundle *bin, std::ostream &ss) {
  const float asxs0 = vector_mean(bin->asxs[0]);
  const float asxs1 = vector_mean(bin->asxs[1]);

  // fa8 modified:
  // int asxs    = custom_round(std::min(asxs0, asxs1));
  int asxs;
  if(bin->asxs[0].size() == 0 && bin->asxs[1].size() != 0) {
  	asxs    = custom_round(asxs1);
  } else if(bin->asxs[0].size() != 0 && bin->asxs[1].size() == 0) {
  	asxs    = custom_round(asxs0);
  } else {
    asxs    = custom_round(std::min(asxs0, asxs1));
  }
  //end

  const float nmms0 = vector_mean(bin->nmms[0]);
  const float nmms1 = vector_mean(bin->nmms[1]);
  // fa8 modified:
  // float nmms  = 0.1 * custom_round(std::max(nmms0, nmms1) * 10.0);
  float nmms;
  if(bin->nmms[0].size() == 0 && bin->nmms[1].size() != 0) {
  	nmms    = custom_round(nmms1);
  } else if(bin->nmms[0].size() != 0 && bin->nmms[1].size() == 0) {
  	nmms    = custom_round(nmms0);
  } else {
    nmms    = custom_round(std::min(nmms0, nmms0));
  }
  //end

  ss
  << asxs
  << "\t"
  << nmms;
}

static inline float get_ppair_mean(const bundle *bin, const int rtype) {
  // ASSUMPTION: if the proper pair count is not zero, neither will be the read count
  return
    bin->rtype_ppair_counts[rtype] == 0 ? 0.0f :
    ((float)bin->rtype_ppair_counts[rtype] / (float)bin->rtype_read_counts[rtype]);
}

// TODO: pass the same stream to each of these functions?
static void push_proper_pair_string(const bundle *bin, std::ostream &ss) {
  const float ppair = custom_round(
    std::min(
      get_ppair_mean(bin, RTYPE_A),
      get_ppair_mean(bin, RTYPE_B)));

  /*
  float ppair0 = MeanOfOneVector(bin->ppair[0]);
  float ppair1 = MeanOfOneVector(bin->ppair[1]);
  //fa8
  // float ppair  = 0.1 * custom_round(std::min(ppair0, ppair1) * 10.0);
  float ppair;
  if(bin->ppair[0].size() == 0 && bin->ppair[1].size() != 0) {
  	ppair    = custom_round(ppair1);
  } else if(bin->ppair[0].size() != 0 && bin->ppair[1].size() == 0) {
  	ppair    = custom_round(ppair0);
  } else {
    ppair    = custom_round(std::min(ppair0, ppair1));
  }
  // end
  */

  ss << ppair;
}

static void push_identifier_string(const bundle *bin, std::ostream &ss) {
  ss << bin->duplex_tag_info.beg;
  ss << "\t";
  ss << bin->duplex_tag_info.end;
  ss << "\t";
  std::string str = bin->duplex_tag_info.fwd_bc;
  std::transform(str.begin(), str.end(),str.begin(), ::toupper);  // uppercase
  ss << str;
  ss << "|";
  str = bin->duplex_tag_info.rev_bc;
  std::transform(str.begin(), str.end(),str.begin(), ::toupper);
  ss << str;
  ss << "\t";
  ss << bin->bundle_type;
}

void WriteOut::WriteRows(bundle bulk, bundles dplx, std::string posn) {
  // Bulk row
  this->out << posn;
  push_bulk_covariates_string(&bulk, this->out);
  this->out << '\t';
  push_counts_string(&bulk, this->out);

  for (auto bndls : dplx) {
    push_identifier_string(&bndls.second, this->out);        this->out << '\t';
    push_duplex_covariates_string(&bndls.second, this->out); this->out << '\t';
    push_counts_string(&bndls.second, this->out);            this->out << '\t';
    push_base_quals_string(&bndls.second, this->out);        this->out << '\t';
    push_proper_pair_string(&bulk, this->out);               this->out << '\t';
    push_proper_pair_string(&bndls.second, this->out);       this->out << std::endl;
  }
}
