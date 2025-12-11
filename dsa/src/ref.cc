#include <format>
#include "ref.h"

void Ref::Init(const char *fai_fp) {
  this->fai = fai_load(fai_fp);
  if (this->fai == NULL) {
    throw std::runtime_error(
      std::format("Failed to open '{}'!", fai_fp));
  }
}

void Ref::Fetch(const char *contig, const range_t range) {
  int32_t seq_length;
  // TODO: verify whether the partitioning step respects the BED conventions...
  this->seq.Set(range, faidx_fetch_seq(
    this->fai, contig, range.start, range.end - 1, &seq_length));

  // TODO: ensure this externally by checking the upper bound as well
  if (range_length(&range) != seq_length) {
    throw std::runtime_error(std::format(
      "Out of bound reference sequence in {}:{}-{}!",
      contig, range.start, range.end));
  }

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
