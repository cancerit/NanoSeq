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

#ifndef DSA_SERIALISER_H_
#define DSA_SERIALISER_H_

#include <sstream>
#include <string>
#include <format>
#include "bulk_bundle.h"
#include "duplex_bundle.h"
#include "duplex_tag_info.h"
#include "pos_stats.h"
#include "utils.h"

inline std::string dsa_base_counts(const uint64_t counts[RTYPE_COUNT][ALLELE_COUNT]) {
    return std::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
        counts[RTYPE_A][ALLELE_A],
        counts[RTYPE_A][ALLELE_C],
        counts[RTYPE_A][ALLELE_G],
        counts[RTYPE_A][ALLELE_T],
        counts[RTYPE_A][ALLELE_DEL],
        counts[RTYPE_B][ALLELE_A],
        counts[RTYPE_B][ALLELE_C],
        counts[RTYPE_B][ALLELE_G],
        counts[RTYPE_B][ALLELE_T],
        counts[RTYPE_B][ALLELE_DEL]);
}

inline std::string dsa_get_bulk_prefix(const bulk_bundle_t *bulk, const pos_final_stats_t *stats, const std::string pos_prefix) {
    return std::format("{}\t{}\t{}\t{}",
        pos_prefix,
        stats->asxs,
        custom_round(stats->nm),
        dsa_base_counts(bulk->allele_counts));
}

inline std::string dsa_get_duplex_base_consensus_qualities(const duplex_bundle_t *b) {
    // ASSUMPTION: the object has been finalised, and therefore the accumulator is actually the final quality score
    return std::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
        custom_round(b->duplex_consensus_quality_accum[RTYPE_A][0]),
        custom_round(b->duplex_consensus_quality_accum[RTYPE_A][1]),
        custom_round(b->duplex_consensus_quality_accum[RTYPE_A][2]),
        custom_round(b->duplex_consensus_quality_accum[RTYPE_A][3]),
        custom_round(b->duplex_consensus_quality_accum[RTYPE_B][0]),
        custom_round(b->duplex_consensus_quality_accum[RTYPE_B][1]),
        custom_round(b->duplex_consensus_quality_accum[RTYPE_B][2]),
        custom_round(b->duplex_consensus_quality_accum[RTYPE_B][3]));
}

inline std::string dsa_get_duplex_id(const duplex_tag_info_t *info, const uint8_t bundle_type) {
    std::string fwd_bc = info->fwd_bc;
    std::string rev_bc = info->rev_bc;
    upper(fwd_bc);
    upper(rev_bc);

    return std::format("{}\t{}\t{}|{}\t{}\t",
        info->beg,
        info->end,
        fwd_bc, rev_bc,
        bundle_type);
}

inline void dsa_push_row(
    std::stringstream &s,
    const duplex_tag_info_t *tag_info,
    const duplex_bundle_t *duplex_bundle,
    const pos_final_stats_t *duplex_stats,
    const pos_final_stats_t *bulk_stats,
    const std::string pos_prefix,
    const uint8_t bundle_type
) {
    s
    << pos_prefix
    << dsa_get_duplex_id(tag_info, bundle_type)
    << std::format("{}\t{}\t{}\t",
        duplex_stats->asxs,
        duplex_stats->clip,
        static_cast<float>(0.1 * custom_round(duplex_stats->nm * 10.0)))
    << dsa_base_counts(duplex_bundle->bundle.allele_counts)
    << dsa_get_duplex_base_consensus_qualities(duplex_bundle)
    << bulk_stats->proper_pairs << '\t'
    << duplex_stats->proper_pairs << '\n';
}

#endif
