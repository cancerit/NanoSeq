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
#include <getopt.h>
#include <header.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#define PROG_NAME "bamaddreadbundles"
#define PROG_VERSION "1.1"

// NOTE: tag case is important
#define TAG_OPT_DUP "od"
#define TAG_BARCODE_BUNDLE "RB"
#define TAG_READ_BARCODE "rb"
#define TAG_MATE_BARCODE "mb"
#define TAG_READ_COORD "rc"
#define TAG_MATE_COORD "mc"


const char* crbam_open_mode(const char* filepath)
{
  static const char* allowedExt[] = {".bam", ".cram"};
  static const char* allowedMode[] = {"wb", "wc"};

  const size_t flen = std::strlen(filepath);
  for (size_t i = 0; i < 2; ++i) {
    const size_t elen = std::strlen(allowedExt[i]);
    if (flen >= elen && std::strcmp(filepath + flen - elen, allowedExt[i]) == 0) {
      return allowedMode[i];
    }
  }
  return nullptr;
}

struct AlnFile {
  htsFile* fh_o = nullptr;
  sam_hdr_t* hdr_o = nullptr;

  ~AlnFile()
  {
    if (hdr_o != nullptr) {
      sam_hdr_destroy(hdr_o);
    }
    if (fh_o != nullptr) {
      // may fail, but destructor should
      // only be relied upon for error path
      // anyway.
      sam_close(fh_o);
    }
  }

  AlnFile() = default;
  AlnFile(const AlnFile&) = delete;
  AlnFile& operator=(const AlnFile&) = delete;
  AlnFile(AlnFile&&) = delete;
  AlnFile& operator=(AlnFile&&) = delete;

  [[nodiscard]] int close() noexcept
  {
    if (hdr_o != nullptr) {
      sam_hdr_destroy(hdr_o);
      hdr_o = nullptr;
    }
    int closeRc = 0;
    if (fh_o != nullptr) {
      closeRc = sam_close(fh_o);
      fh_o = nullptr;
    }

    return closeRc;
  }
};

// NOTE: compiler hints may or may not have any impact.
struct MemArenas {
  constexpr static uint16_t arenaMaxSz = 10000;
  static bam1_t recArena[arenaMaxSz];
  static uint8_t* readCoordArena[arenaMaxSz];
  static uint8_t* mateCoordArena[arenaMaxSz];
  static uint8_t* readBarcdArena[arenaMaxSz];
  static uint8_t* mateBarcdArena[arenaMaxSz];
  static bool tagFilterHitArena[arenaMaxSz];

  static uint16_t arenaI;
  static uint16_t arenaN;
  static bool everSeenTags;
  static bool applyInputFilter;
  constexpr static size_t maxFilterTags = 64;
  static std::set<std::array<char, 2>> filterTags;

  static bool matches_any_filter_tag(const bam1_t& rec) noexcept
  {
    for (const auto& tag : filterTags) {
      if (bam_aux_get(&rec, tag.data()) != nullptr) {
        return true;
      }
    }
    return false;
  }

  enum class LoadRetCode : uint8_t { fullLoad, eofLoad, readFail, corruptTag, existingRB };
  static LoadRetCode load_batch(AlnFile& aln) noexcept
  {
    errno = 0;
    arenaI = 0;
    arenaN = 0;
    for (; arenaI < arenaMaxSz;) {
      auto& rec = recArena[arenaI];
      const auto read1Rc = sam_read1(aln.fh_o, aln.hdr_o, &rec);
      if (read1Rc == -1) [[unlikely]] {
        arenaN = arenaI;
        return LoadRetCode::eofLoad;
      }
      if (read1Rc < -1) [[unlikely]] {
        return LoadRetCode::readFail;
      }

      readCoordArena[arenaI] = bam_aux_get(&rec, TAG_READ_COORD);
      mateCoordArena[arenaI] = bam_aux_get(&rec, TAG_MATE_COORD);
      readBarcdArena[arenaI] = bam_aux_get(&rec, TAG_READ_BARCODE);
      mateBarcdArena[arenaI] = bam_aux_get(&rec, TAG_MATE_BARCODE);
      if (readCoordArena[arenaI] || mateCoordArena[arenaI] || readBarcdArena[arenaI] ||
          mateBarcdArena[arenaI]) {
        everSeenTags = true;
      }
      if (bam_aux_get(&rec, TAG_BARCODE_BUNDLE) != nullptr) {
        return LoadRetCode::existingRB;
      }
      if (errno == EINVAL) [[unlikely]] {
        return LoadRetCode::corruptTag;
      }
      const bool qcFail = rec.core.flag & BAM_FQCFAIL;
      tagFilterHitArena[arenaI] = !qcFail && matches_any_filter_tag(rec);
      if (applyInputFilter) {
        if (qcFail || tagFilterHitArena[arenaI] || readBarcdArena[arenaI] == nullptr ||
            mateBarcdArena[arenaI] == nullptr) {
          continue;  // overwrite
        }
      }
      ++arenaI;
    }
    arenaN = arenaI;
    return LoadRetCode::fullLoad;
  }

  /* Add barcode bundle RB tag */
  enum class ProcessRetCode : int8_t { success = 0, tagWriteErr = -1, noTid = -2, tagBadType = -3 };
  static ProcessRetCode modify_tags_batch(sam_hdr_t* hdr_br)
  {
    constexpr static auto flagFailBits =
        BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FSECONDARY | BAM_FUNMAP;

    errno = 0;
    ProcessRetCode rc = ProcessRetCode::success;
    std::string tagBuf;
    std::array<int64_t, 2> coBuf;
    std::array<uint8_t*, 4> tagDelBuf;
    arenaI = 0;
    for (; arenaI < arenaN; ++arenaI) {
      auto& rec = recArena[arenaI];
      auto& readCoordTag = readCoordArena[arenaI];
      auto& mateCoordTag = mateCoordArena[arenaI];
      auto& readBarcdTag = readBarcdArena[arenaI];
      auto& mateBarcdTag = mateBarcdArena[arenaI];

      if (!(rec.core.flag & BAM_FPROPER_PAIR)) {
        continue;
      }
      if (rec.core.flag & flagFailBits || tagFilterHitArena[arenaI]) {
        continue;
      }
      if (readCoordTag == nullptr || mateCoordTag == nullptr || readBarcdTag == nullptr ||
          mateBarcdTag == nullptr) {
        continue;
      }

      tagBuf.clear();

      const auto* tidName = sam_hdr_tid2name(hdr_br, rec.core.tid);
      if (tidName == nullptr) [[unlikely]] {
        rc = ProcessRetCode::noTid;
        break;
      }
      tagBuf += tidName;
      tagBuf += ",";

      coBuf[0] = bam_aux2i(readCoordTag);
      coBuf[1] = bam_aux2i(mateCoordTag);
      std::sort(begin(coBuf), end(coBuf));
      tagBuf += std::to_string(coBuf[0]);
      tagBuf += ",";
      tagBuf += std::to_string(coBuf[1]);
      tagBuf += ",";

      if (rec.core.flag & BAM_FREVERSE) {
        tagBuf += bam_aux2Z(mateBarcdTag);
        tagBuf += ",";
        tagBuf += bam_aux2Z(readBarcdTag);
      }
      else {
        tagBuf += bam_aux2Z(readBarcdTag);
        tagBuf += ",";
        tagBuf += bam_aux2Z(mateBarcdTag);
      };

      if (errno == EINVAL) [[unlikely]] {
        // bad tag type
        rc = ProcessRetCode::tagBadType;
        break;
      }

      // Deleting a tag shifts every *subsequent* byte in rec.data, which
      // invalidates any other cached pointer into rec aux blob.
      // Deleting in descending address order means each delete only ever
      // shifts bytes above any unprocessed pointers.
      tagDelBuf[0] = readCoordTag;
      tagDelBuf[1] = mateCoordTag;
      tagDelBuf[2] = readBarcdTag;
      tagDelBuf[3] = mateBarcdTag;
      std::sort(tagDelBuf.begin(), tagDelBuf.end(), std::greater<uint8_t*>());
      for (auto* tag : tagDelBuf) {
        if (bam_aux_del(&rec, tag) < 0) [[unlikely]] {
          rc = ProcessRetCode::tagWriteErr;
          break;
        }
      }
      if (rc != ProcessRetCode::success) [[unlikely]] {
        break;
      }

      rc = static_cast<ProcessRetCode>(bam_aux_append(
          &rec, TAG_BARCODE_BUNDLE, 'Z', tagBuf.length() + 1,
          reinterpret_cast<const uint8_t*>(tagBuf.c_str())
      ));
      if (rc != ProcessRetCode::success) [[unlikely]] {
        break;
      }
    }
    return rc;
  }

  static bool write_batch(AlnFile& out)
  {
    arenaI = 0;
    for (; arenaI < arenaN; ++arenaI) {
      if (sam_write1(out.fh_o, out.hdr_o, &recArena[arenaI]) < 0) {
        return false;
      }
    }
    arenaI = 0;
    arenaN = 0; // reset
    return true;
  }
};
// C++11 mandates that static
// member init is out-of-line...
uint16_t MemArenas::arenaI = 0;
uint16_t MemArenas::arenaN = 0;
bool MemArenas::everSeenTags = false;
bool MemArenas::applyInputFilter = true;
std::set<std::array<char, 2>> MemArenas::filterTags;
bam1_t MemArenas::recArena[MemArenas::arenaMaxSz];
uint8_t* MemArenas::readCoordArena[MemArenas::arenaMaxSz];
uint8_t* MemArenas::mateCoordArena[MemArenas::arenaMaxSz];
uint8_t* MemArenas::readBarcdArena[MemArenas::arenaMaxSz];
uint8_t* MemArenas::mateBarcdArena[MemArenas::arenaMaxSz];
bool MemArenas::tagFilterHitArena[MemArenas::arenaMaxSz];

struct CLIArgs {
  const char* alnInPath = nullptr;
  const char* bamOutPath = nullptr;
  bool uncompressed = false;
};
static constexpr const char* usage{
    "\nUsage:\n"
    "\t-I, --input\tInput BAM/CRAM file name\n"
    "\t-O, --output\tOutput BAM/CRAM file name\n"
    "\t-n, --no-filter\tDisable output filter.\n"
    "\t\tInclude reads flagged QC fail or with\n"
    "\t\tany -t/--filter-tag tag in output.\n"
    "\t-t, --filter-tag\tAux tag whose presence marks a read for\n"
    "\t\texclusion from output (e.g. optical duplicates).\n"
    "\t\tMay be given multiple times. Defaults to \"od\".\n"
    "\t-u, --uncompressed\tWrite uncompressed output.\n"
    "\t-h, --help\tHelp\n"
};
static const struct option longOpts[] = {
    {"input", required_argument, nullptr, 'I'}, {"output", required_argument, nullptr, 'O'},
    {"no-filter", no_argument, nullptr, 'n'},   {"filter-tag", required_argument, nullptr, 't'},
    {"uncompressed", no_argument, nullptr, 'u'}, {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0},
};
// NOTE: chunk/phase-separated approach
// should make multi-threading very easy.
int main(int argc, char** argv)
{
  CLIArgs args;

  int opt = 0;
  while ((opt = getopt_long(argc, argv, "I:O:nt:uh", longOpts, nullptr)) >= 0) {
    switch (opt) {
      case 'I':
        args.alnInPath = optarg;
        break;
      case 'O':
        args.bamOutPath = optarg;
        break;
      case 'n':
        MemArenas::applyInputFilter = false;
        break;
      case 't':
        if (std::strlen(optarg) != 2) {
          std::cerr << "Usage error: -t/--filter-tag must be exactly 2 characters, got \"" << optarg
                    << "\"" << std::endl;
          return EXIT_FAILURE;
        }
        if (std::strncmp(optarg, TAG_READ_COORD, 2) == 0 ||
            std::strncmp(optarg, TAG_MATE_COORD, 2) == 0 ||
            std::strncmp(optarg, TAG_READ_BARCODE, 2) == 0 ||
            std::strncmp(optarg, TAG_MATE_BARCODE, 2) == 0 ||
            std::strncmp(optarg, TAG_BARCODE_BUNDLE, 2) == 0) {
          std::cerr << "Usage error: -t/--filter-tag \"" << optarg
                    << "\" is a tag reserved for internal use by " << PROG_NAME << std::endl;
          return EXIT_FAILURE;
        }
        if (MemArenas::filterTags.size() >= MemArenas::maxFilterTags) {
          std::cerr << "Usage error: too many -t/--filter-tag options (max "
                    << MemArenas::maxFilterTags << ")" << std::endl;
          return EXIT_FAILURE;
        }
        MemArenas::filterTags.insert({optarg[0], optarg[1]});
        break;
      case 'u':
        args.uncompressed = true;
        break;
      case 'h':
        std::cerr << usage << std::endl;
        return EXIT_SUCCESS;
      default:
        std::cerr << "Unknown option" << std::endl;
        return EXIT_FAILURE;
    }
  }
  if (MemArenas::filterTags.empty()) {
    MemArenas::filterTags.insert({TAG_OPT_DUP[0], TAG_OPT_DUP[1]});
  }
  if (args.alnInPath == nullptr) {
    std::cerr << "Usage error: no input file specified" << std::endl;
    return EXIT_FAILURE;
  }
  if (args.bamOutPath == nullptr) {
    std::cerr << "Usage error: no output file specified" << std::endl;
    return EXIT_FAILURE;
  }
  const char* outMode = crbam_open_mode(args.bamOutPath);
  if (outMode == nullptr) {
    std::cerr << "Usage error: output extension must be .bam or .cram" << std::endl;
    return EXIT_FAILURE;
  }
  std::string outModeStr = outMode;
  if (args.uncompressed) {
    if (outModeStr == "wc") {
      std::cerr << "Warning: -u/--uncompressed has no effect on CRAM output" << std::endl;
    }
    else {
      outModeStr += 'u';
    }
  }

  // open input handle
  AlnFile alnIn;
  alnIn.fh_o = hts_open(args.alnInPath, "r");
  if (alnIn.fh_o == nullptr) {
    std::cerr << "Error: failed to open input alignment file at" << args.alnInPath << std::endl;
    return EXIT_FAILURE;
  }
  alnIn.hdr_o = sam_hdr_read(alnIn.fh_o);
  if (alnIn.hdr_o == nullptr) {
    std::cerr << "Error: failed to read header from input alignment" << std::endl;
    return EXIT_FAILURE;
  }
  if (alnIn.hdr_o->n_targets < 1) {
    std::cerr << "Error: input alignment header contains no contigs" << std::endl;
    return EXIT_FAILURE;
  }

  // open output handle
  AlnFile alnOut;
  alnOut.fh_o = hts_open(args.bamOutPath, outModeStr.c_str());
  if (alnOut.fh_o == nullptr) {
    std::cerr << "Error: failed to open output alignment file at " << args.bamOutPath << std::endl;
    return EXIT_FAILURE;
  }

  // create output header
  alnOut.hdr_o = sam_hdr_dup(alnIn.hdr_o);
  if (alnOut.hdr_o == nullptr) {
    std::cerr << "Error: failed to duplicate input header" << std::endl;
    return EXIT_FAILURE;
  }
  std::stringstream cmdLineCall;
  for (int i = 0; i < argc; ++i) {
    cmdLineCall << argv[i];
    cmdLineCall << " ";
  }
  if (sam_hdr_add_pg(
          alnOut.hdr_o, PROG_NAME, "VN", PROG_VERSION, "CL", cmdLineCall.str().c_str(), NULL
      ) < 0) {
    std::cerr << "Error: failed to update header" << std::endl;
    return EXIT_FAILURE;
  }
  if (sam_hdr_rebuild(alnOut.hdr_o) < 0) {
    std::cerr << "Error: failed to rebuild header" << std::endl;
    return EXIT_FAILURE;
  }
  if (sam_hdr_write(alnOut.fh_o, alnOut.hdr_o) < 0) {
    std::cerr << "Error: failed to write header to output alignment file" << std::endl;
    return EXIT_FAILURE;
  }

  std::cerr << PROG_NAME << ": begin processing" << std::endl;
  errno = 0;
  while (true) {
    const auto loadRc = MemArenas::load_batch(alnIn);
    switch (loadRc) {
      case MemArenas::LoadRetCode::fullLoad:
      case MemArenas::LoadRetCode::eofLoad:
        break;
      case MemArenas::LoadRetCode::readFail:
        std::cerr << "Error: failure while reading input alignment file" << std::endl;
        return EXIT_FAILURE;
      case MemArenas::LoadRetCode::corruptTag:
        std::cerr << "Error: input alignment contains corrupt tag data at read "
                  << bam_get_qname(&MemArenas::recArena[MemArenas::arenaI]) << std::endl;
        return EXIT_FAILURE;
      case MemArenas::LoadRetCode::existingRB:
        std::cerr << "Error: input alignment already carries an RB tag at read "
                  << bam_get_qname(&MemArenas::recArena[MemArenas::arenaI])
                  << "; has bamaddreadbundles already been run on this input?" << std::endl;
        return EXIT_FAILURE;
    }

    // modify tags
    switch (MemArenas::modify_tags_batch(alnIn.hdr_o)) {
      case MemArenas::ProcessRetCode::success:
        break;
      case MemArenas::ProcessRetCode::tagWriteErr:
        std::cerr << "Error: failed to write tag to record." << std::endl;
        return EXIT_FAILURE;
      case MemArenas::ProcessRetCode::noTid:
        std::cerr << "Error: read " << bam_get_qname(&MemArenas::recArena[MemArenas::arenaI])
                  << "has invalid tid" << std::endl;
        return EXIT_FAILURE;
      case MemArenas::ProcessRetCode::tagBadType:
        std::cerr << "Error: read " << bam_get_qname(&MemArenas::recArena[MemArenas::arenaI])
                  << "has corrupt or incorrectly typed tag" << std::endl;
        return EXIT_FAILURE;
    }

    if (!MemArenas::write_batch(alnOut)) {
      std::cerr << "Error: failure while writing record" << std::endl;
      return EXIT_FAILURE;
    }

    if (loadRc == MemArenas::LoadRetCode::eofLoad) {
      // was last batch
      break;
    }
  }

  if (alnOut.close() != 0) {
    std::cerr
        << "Error: error while closing output alignment file; check file data integrity before use"
        << std::endl;
    return EXIT_FAILURE;
  };

  if (!MemArenas::everSeenTags) {
    std::cerr << "Warning: no reads in the input carried rc/mc/rb/mb tags; check the input has "
                 "been through the upstream barcode-tagging step."
              << std::endl;
  }

  std::cerr << PROG_NAME << ": complete" << std::endl;
  return EXIT_SUCCESS;
}
