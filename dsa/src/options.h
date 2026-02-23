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

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include "constants.h"
#include <filesystem>
#include <format>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "utils.h"

struct Options {
  const char *bams[BAM_COUNT] = {NULL, NULL};
  const char *beds[MASK_COUNT] = {NULL, NULL};
  const char *ranges_bed = NULL;
  const char *fasta = NULL;
  const char *oname = NULL;
  int min_base_quality;
  int min_mapQ;
  int max_plp_depth;
  int min_dplx_depth;
  int max_dplx_depth;
  int offset;
  bool doTests;
  bool debug_mode;
  int compression_level;
};

static int options_validate(const Options *opt) {

    if (opt->oname == NULL) {
        std::cerr << std::format("Output directory not set!\n");
        return 1;
    }

    if (!path_exists(opt->oname)) {
        std::cerr << std::format(
            "Output directory not found at '{}'!\n", opt->oname);
        return 1;
    }

    {
        const char *bed_fp = NULL;
        for (int i = 0; i < MASK_COUNT; ++i) {
            bed_fp = opt->beds[i];

            if (bed_fp == NULL) {
                std::cerr << std::format(
                    "Mask file path not set!\n");
                return 1;
            }

            if (!path_exists(bed_fp)) {
                std::cerr << std::format(
                    "Mask file not found at '{}'!\n", bed_fp);
                return 1;
            }
        }
    }

    {
        const char *bam_fp = NULL;
        for (int i = 0; i < BAM_COUNT; ++i) {
            bam_fp = opt->bams[i];

            if (bam_fp == NULL) {
                std::cerr << std::format(
                    "BAM/CRAM file path not set!\n");
                return 1;
            }

            if (!path_exists(bam_fp)) {
                std::cerr << std::format(
                    "BAM/CRAM file not found at '{}'!\n", bam_fp);
                return 1;
            }
        }
    }

    // TODO: is level zero valid?
    if (opt->compression_level < 0 || opt->compression_level > 12) {
        std::cerr << std::format(
            "Invalid compression level %d (should be in [1, 12])!\n",
            opt->compression_level);
        return 1;
    }

    if (opt->debug_mode) {
        std::cerr << "Verbose mode (debug)\n";
    }
    return 0;
}

static FILE *options_open_output_file(const Options *opt, const char *fn) {
    const std::filesystem::path fp = std::filesystem::path(opt->oname).append(fn);
    FILE *f = fopen(fp.c_str(), "w");
    if (f == NULL) {
        throw std::runtime_error(std::format(
            "Failed to open output file '{}'!", fp.c_str()));
    }
    return f;
}

static FILE *options_open_output_debug_file(const Options *opt, const char *fn) {
    return opt->debug_mode ? options_open_output_file(opt, fn) : NULL;
}

#endif  // OPTIONS_H_
