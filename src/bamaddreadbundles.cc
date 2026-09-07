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

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cerrno>
#include <cstdint>
#include <cstdlib>

#define PROG_NAME "bamaddreadbundles"
#define PROG_VERSION "1.1"

// NOTE: tag case is important
#define TAG_OPT_DUP "od"
#define TAG_BARCODE_BUNDLE "RB"
#define TAG_READ_BARCODE "rb"
#define TAG_MATE_BARCODE "mb"
#define TAG_READ_COORD "rc"
#define TAG_MATE_COORD "mc"
#define TAG_MATE_MAPQ "MQ"
#define TAG_MATE_SCORE "ms"
#define TAG_MATE_CIGAR "MC"


static bool ReadHasAux(const bam1_t* b, const char* tag)
{
  // does not check if absent, or error.
  return bam_aux_get(b, tag) != NULL;
}
static bool ReadHasAux(const bam1_t* b, const std::vector<const char*>& tags)
{
  for (const auto& t : tags) {
    if (bam_aux_get(b, t) == NULL) {
      return false;
    }
  }
  return true;
}

static bool ReadIsUsable(bam1_t* b)
{
  constexpr auto fail_bits = BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FUNMAP | BAM_FSECONDARY;

  bool out = false;

  const auto has_requisite_tags =
      ReadHasAux(b, {TAG_READ_COORD, TAG_MATE_COORD, TAG_READ_BARCODE, TAG_MATE_BARCODE});

  if (b->core.flag & fail_bits || ReadHasAux(b, TAG_OPT_DUP)) {
    out = false;
  }
  else if ((b->core.flag & BAM_FPROPER_PAIR) && has_requisite_tags) {
    out = true;
  }
  return out;
}

static void DelSupersededTags(bam1_t* b)
{
  uint8_t* tag;
  if ((tag = bam_aux_get(b, TAG_MATE_COORD)) != NULL) {
    bam_aux_del(b, tag);
  }
  if ((tag = bam_aux_get(b, TAG_READ_COORD)) != NULL) {
    bam_aux_del(b, tag);
  }
  if ((tag = bam_aux_get(b, TAG_MATE_BARCODE)) != NULL) {
    bam_aux_del(b, tag);
  }
  if ((tag = bam_aux_get(b, TAG_READ_BARCODE)) != NULL) {
    bam_aux_del(b, tag);
  }
}

static bool CheckWriteFilters(bam1_t* b)
{
  bool out = false;
  if (ReadHasAux(b, TAG_BARCODE_BUNDLE)) {
    out = true;
  }
  else if (b->core.flag & BAM_FQCFAIL || ReadHasAux(b, TAG_OPT_DUP)) {
    out = false;
  }
  else if (ReadHasAux(b, {TAG_READ_BARCODE, TAG_MATE_BARCODE})) {
    out = true;
  }
  return out;
}


void BamAddReadBundles::AddBarcodeBundleAuxTag(bam1_t* b)
{
  int rc = bam_aux2i(bam_aux_get(b, TAG_READ_COORD));
  int mc = bam_aux2i(bam_aux_get(b, TAG_MATE_COORD));
  char* rb = bam_aux2Z(bam_aux_get(b, TAG_READ_BARCODE));
  char* mb = bam_aux2Z(bam_aux_get(b, TAG_MATE_BARCODE));
  int strand = (b->core.flag & BAM_FREVERSE) ? 1 : 0;
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
  }
  else {
    ss << mb;
    ss << ",";
    ss << rb;
  }
  std::string str = ss.str();
  const char* cstr = str.c_str();
  int len = str.length() + 1;
  uint8_t* data = const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(cstr));
  int rco = bam_aux_append(b, TAG_BARCODE_BUNDLE, 'Z', len, data);
  if (rco < 0) {
    exit(1);
  }
}


void BamAddReadBundles::WriteOut(bam1_t* b)
{
  if ((sam_write1(this->out, this->head, b) < 0)) {
    throw std::runtime_error("Error: failed to write record.");
  }
}


void BamAddReadBundles::FilterAndTagReads()
{
  bam1_t* b = bam_init1();
  int ret;
  while (1) {
    ret = sam_read1(this->in, this->head, b);
    if (ret == -1) {
      break;
    }
    if (ret < -1) {
      std::stringstream er;
      er << "Error: failure while reading input BAM";
      er << std::endl;
      throw std::runtime_error(er.str());
    }
    if (ReadIsUsable(b)) {
      BamAddReadBundles::AddBarcodeBundleAuxTag(b);
      DelSupersededTags(b);
    }
    if (write_all_reads || CheckWriteFilters(b)) {
      BamAddReadBundles::WriteOut(b);
    }
  }
  bam_destroy1(b);
}


static void Usage()
{
  fprintf(stderr, "\nUsage:\n");
  fprintf(stderr, "\t-I\tInput BAM/CRAM file name\n");
  fprintf(stderr, "\t-O\tOutput BAM/CRAM file name\n");
  fprintf(
      stderr,
      "\t-n\tDo not filter any reads from output."
      "\t\t  Defaults False, excluding reads marked"
      "\t\t  QC fail or with od Optical duplicate tag."
  );
  fprintf(stderr, "\t-h\tHelp\n");
}


bool has_crbam_ext(const char* filepath)
{
  static const char* allowed[] = {".bam", ".cram"};

  const uint16_t flen = std::strlen(filepath);
  for (const char* ext : allowed) {
    const uint16_t elen = std::strlen(ext);
    if (flen >= elen && std::strcmp(filepath + flen - elen, ext) == 0) {
      return true;
    }
  }
  return false;
}

enum class FilterMode : uint8_t { off, on };
struct CLIArgs {
  const char* alnInPath = nullptr;
  const char* bamOutPath = nullptr;
  FilterMode inputFilter = FilterMode::on;
};

struct MemArenas {
  constexpr static uint16_t arenaMaxSz = 10000;
  static bam1_t recArena[arenaMaxSz];
  static int64_t readCoordArena[arenaMaxSz];
  static int64_t mateCoordArena[arenaMaxSz];
  static const char* readBarcdArena[arenaMaxSz];
  static const char* mateBarcdArena[arenaMaxSz];
};

struct AlnFile {
  htsFile* fh_o = nullptr;
  sam_hdr_t* hdr_o = nullptr;

  ~AlnFile () {
    if (hdr_o != nullptr) {
      sam_hdr_destroy(hdr_o);
    }
    if (fh_o != nullptr) {
      // BUG/TODO/NOTE: destructor can't return error...
      // Deferred cleanup would be better
      sam_close(fh_o);
    }
  }

  AlnFile() = default;
  AlnFile(const AlnFile&) = delete;
  AlnFile& operator=(const AlnFile&) = delete;
  AlnFile(AlnFile&&) = delete;
  AlnFile& operator=(AlnFile&&) = delete;
};

static int load_arenas (AlnFile& aln, FilterMode mode) {
  uint16_t arenaI = 0;
  for (; arenaI < MemArenas::arenaMaxSz;) {
    auto& rec = MemArenas::recArena[arenaI];
    const auto read1Rc = sam_read1(aln.fh_o, aln.hdr_o, &rec);
    if (read1Rc == -1) {
      break;  // EOF
    }
    if (read1Rc < -1) {
      return read1Rc;
    }

    MemArenas::readCoordArena[arenaI] = bam_aux2i(bam_aux_get(&rec, TAG_READ_COORD));
    MemArenas::mateCoordArena[arenaI] = bam_aux2i(bam_aux_get(&rec, TAG_MATE_COORD));
    MemArenas::readBarcdArena[arenaI] = bam_aux2Z(bam_aux_get(&rec, TAG_READ_BARCODE));
    MemArenas::mateBarcdArena[arenaI] = bam_aux2Z(bam_aux_get(&rec, TAG_MATE_BARCODE));
    if (mode == FilterMode::off) {
      ++arenaI;
      continue;
    }
    else {
      if (rec.core.flag & BAM_FQCFAIL || bam_aux_get(&rec, TAG_OPT_DUP) != nullptr ||
          (MemArenas::readBarcdArena[arenaI] == nullptr && MemArenas::mateBarcdArena[arenaI] == nullptr)) {
        continue;  // overwrite
      }
    }
  }
  return arenaI;
}

int main(int argc, char** argv)
{
  CLIArgs args;

  int opt = 0;
  while ((opt = getopt(argc, argv, "I:O:n:h")) >= 0) {
    switch (opt) {
      case 'I':
        args.alnInPath = optarg;
        break;
      case 'O':
        args.bamOutPath = optarg;
        break;
      case 'n':
        args.inputFilter = FilterMode::off;
      case 'h':
        Usage();
        exit(EXIT_SUCCESS);
      default:
        break;
    }
  }

  if (args.alnInPath == nullptr) {
    std::cerr << "Error: no input file specified";
    return EXIT_FAILURE;
  }
  if (args.bamOutPath == nullptr) {
    std::cerr << "Error: no output file specified";
    return EXIT_FAILURE;
  }
  if (!has_crbam_ext(args.bamOutPath)) {
    std::cerr << "Error: output extension must be .bam or .cram";
    return EXIT_FAILURE;
  }

  // open input handle
  AlnFile alnIn;
  alnIn.fh_o = hts_open(args.alnInPath, "r");
  if (alnIn.fh_o == nullptr) {
    std::cerr << "Error: failed to open input alignment file at" << args.alnInPath;
    return EXIT_FAILURE;
  }
  alnIn.hdr_o = sam_hdr_read(alnIn.fh_o);
  if (alnIn.hdr_o) {
    std::cerr << "Error: failed to read header from input alignment";
    return EXIT_FAILURE;
  }
  if (alnIn.hdr_o ->n_targets) {
    std::cerr << "Error: input alignment header contains no contigs";
    return EXIT_FAILURE;
  }

  // open output handle
  AlnFile alnOut;
  alnOut.fh_o = hts_open(args.bamOutPath, "w");
  if (alnOut.fh_o == nullptr) {
    std::cerr << "Error: failed to open output alignment file at " << args.bamOutPath;
    return EXIT_FAILURE;
  }

  // create output header
  std::stringstream cmdLineCall;
  for (uint16_t i = 0; i < argc; ++i) {
    cmdLineCall << argv[i];
    cmdLineCall << " ";
  }

  alnOut.hdr_o = sam_hdr_dup(alnIn.hdr_o);
  if (sam_hdr_add_pg(alnOut.hdr_o, PROG_NAME, "VN", PROG_VERSION, "CL", cmdLineCall.str().c_str(), NULL) <
      0) {
    std::cerr << "Error: failed to update header";
    return EXIT_FAILURE;
  }
  if (sam_hdr_rebuild(alnOut.hdr_o) < 0) {
    std::cerr << "Error: failed to rebuild header";
    return EXIT_FAILURE;
  }
  if (sam_hdr_write(alnOut.fh_o, alnOut.hdr_o) < 0) {
    std::cerr << "Error: failed to write header to output alignment file";
    return EXIT_FAILURE;
  }

  while (true) {
    errno = 0;
    const auto arenaN = load_arenas(alnIn, args.inputFilter);
    if (arenaN < 0) {
      std::cerr << "Error: failure while reading input alignment file";
      return EXIT_FAILURE;
    }
    if (arenaN < MemArenas::arenaMaxSz) {
      // EOF, last batch
      break;
    }
    if (errno == EINVAL) {
      // aux_get or _aux2* encountered a corrupt tag
      std::cerr << "Error: input alignment contains corrupt tag data";
      return EXIT_FAILURE;
    }

    // modify tags
    
  }


  return EXIT_FAILURE;
}
