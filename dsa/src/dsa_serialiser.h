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

static inline std::string dsa_base_counts(const uint64_t counts[RTYPE_COUNT][ALLELE_COUNT]) {
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

static inline std::string dsa_get_bulk_prefix(const bulk_bundle_t *bulk, const pos_final_stats_t *stats, const std::string pos_prefix) {
    return std::format("{}\t{}\t{}\t{}",
        pos_prefix,
        stats->asxs,
        custom_round(stats->nm),
        dsa_base_counts(bulk->allele_counts));
}

static inline std::string dsa_get_duplex_base_consensus_qualities(const duplex_bundle_t *b) {
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

static inline std::string dsa_get_duplex_id(const duplex_tag_info_t *info, const uint8_t bundle_type) {
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

static inline void dsa_push_row(
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
