/*########## LICENCE ##########
# Copyright (c) 2022, 2025, 2026 Genome Research Ltd
#
# Authors: Luca Barbon <lb29@sanger.ac.uk>, Alex Byrne <ab63@sanger.ac.uk>
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
#include <string.h>  // memchr
#include <cstring>  // memchr
#include "ref.h"

void Ref::Init(const char *fai_fp) {
  this->fai = fai_load(fai_fp);
  if (this->fai == NULL) {
    throw std::runtime_error(
      std::format("Failed to open '{}'!", fai_fp));
  }
}

void Ref::Fetch(const char *contig, const range_t range_) {
  int32_t seq_length;
  // TODO: verify whether the partitioning step respects the BED conventions...

  const int32_t max_pos = faidx_seq_len(this->fai, contig) - 1;
  assert(max_pos > 0);
  range_t range = range_;
  if (range.start < 0) {
    range.start = 0;
  }
  if (range.end > max_pos) {
    range.end = max_pos;
  }

  this->seq.Set(range, faidx_fetch_seq(
    this->fai, contig, range.start, range.end, &seq_length));

  // Validate reference sequence
  {
    const char *ref_seq = this->seq.Data();
    if (memchr(ref_seq, '\n', static_cast<size_t>(seq_length)) != nullptr) {
      throw std::runtime_error(
        "New line characters in reference sequence! "
        "Check for FASTA vs. FAI mismatch!");
    }
  }

  // TODO: ensure this externally by checking the upper bound as well
  /*
  if (range_length(&range) != seq_length) {
    throw std::runtime_error(std::format(
      "Out of bound reference sequence in {}:{}-{}!",
      contig, range.start, range.end));
  }
  */

  switch (seq_length) {
  case 0:
    throw std::runtime_error(std::format(
      "Empty sequence in {}:{}-{}!",
      contig, range.start, range.end));
  case -1:
    throw std::runtime_error(
      std::format("Contig '{}' not found in FAI!", contig));
  case -2:
    throw std::runtime_error(std::format(
      "Failed to fetch sequence in {}:{}-{}!",
      contig, range.start, range.end));
  default:
    break;
  }

  if (this->seq.IsNull()) {
    throw std::runtime_error(std::format(
      "Failed to fetch sequence in {}:{}-{}!",
      contig, range.start, range.end));
  }
}

char *Ref::From(const int32_t pos) {
  return this->seq.From(pos);
}

const std::string Ref::ToString() {
  return seq.ToString();
}

const std::string_view Ref::GetTripletAround(const int32_t pos) {
  return std::string_view(this->From(pos - 1), 3);
}
