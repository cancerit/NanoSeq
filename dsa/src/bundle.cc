#include "bundle.h"
#include "utils.h"
#include <format>

// #define POWER 10.0

static inline double bundle_mean(const bundle_open_t *b, const uint64_t r_type, const int64_t x[2]) {
    return static_cast<double>(x[r_type]) / static_cast<double>(b->rtype_read_counts[r_type]);
}

/*
static inline double mean(const double a, const double b) {
    return b != 0.0 ? a / b : 0.0;
}

static inline double bundle_mean_or_zero(const bundle_open_t *b, const uint64_t r_type, const int64_t x[2]) {
    return b->rtype_read_counts[r_type] == 0 ? 0.0 : bundle_mean(b, r_type, x);
}
*/

static inline void bundle_closed_set_common(bundle_closed_t *s, const bundle_open_t *b) {
    const bool has_a = b->rtype_ppair_counts[RTYPE_A] != 0;
    const bool has_b = b->rtype_ppair_counts[RTYPE_B] != 0;

    // Mean AS - XS and proper pair counts
    if (has_a && has_b) {

        s->asxs = std::min(
            bundle_mean(b, RTYPE_A, b->asxs_accum),
            bundle_mean(b, RTYPE_B, b->asxs_accum));

        s->proper_pairs = std::min(
            bundle_mean(b, RTYPE_A, b->rtype_ppair_counts),
            bundle_mean(b, RTYPE_B, b->rtype_ppair_counts));

    } else if (has_a) {
        s->asxs = bundle_mean(b, RTYPE_A, b->asxs_accum);
        s->proper_pairs = bundle_mean(b, RTYPE_A, b->rtype_ppair_counts);
    } else if (has_b) {
        s->asxs = bundle_mean(b, RTYPE_B, b->asxs_accum);
        s->proper_pairs = bundle_mean(b, RTYPE_B, b->rtype_ppair_counts);
    } else {
        s->asxs = 0.0;
        s->proper_pairs = 0.0;
    }
}

void bundle_closed_bulk_init(bundle_closed_t *s, const bundle_open_t *b) {
    const bool has_a = b->rtype_ppair_counts[RTYPE_A] != 0;
    const bool has_b = b->rtype_ppair_counts[RTYPE_B] != 0;
    bundle_closed_set_common(s, b);

    // NM
    if (has_a && has_b) {
        // TODO: should this be max(a, b) instead (same as duplex NM)?
        s->nm = bundle_mean(b, RTYPE_A, b->nmms_accum);
    } else if (has_a) {
        s->nm = bundle_mean(b, RTYPE_A, b->nmms_accum);
    } else if (has_b) {
        s->nm = bundle_mean(b, RTYPE_B, b->nmms_accum);
    } else {
        s->nm = 0.0f;
    }
}

void bundle_closed_duplex_init(bundle_closed_t *s, const bundle_open_t *b) {
    const bool has_a = b->rtype_ppair_counts[RTYPE_A] != 0;
    const bool has_b = b->rtype_ppair_counts[RTYPE_B] != 0;
    bundle_closed_set_common(s, b);

    // NM
    if (has_a && has_b) {
        s->nm = std::max(
            bundle_mean(b, RTYPE_A, b->nmms_accum),
            bundle_mean(b, RTYPE_B, b->nmms_accum));
    } else if (has_a) {
        s->nm = bundle_mean(b, RTYPE_A, b->nmms_accum);
    } else if (has_b) {
        s->nm = bundle_mean(b, RTYPE_B, b->nmms_accum);
    } else {
        s->nm = 0.0f;
    }

    const uint64_t total = b->rtype_read_counts[RTYPE_A] + b->rtype_read_counts[RTYPE_B];

    // 5'-clipping
    if (total != 0) {
        s->clip = static_cast<double>(b->clip_accum[RTYPE_A] + b->clip_accum[RTYPE_B]) / static_cast<double>(total);
    } else {
        s->clip = 0.0;
    }

}

/* SERIALISATION */

static inline std::string dsa_duplex_base_counts(const duplex_base_t *b, const int bundle_type_index) {
    return std::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
        b->counts[bundle_type_index][RTYPE_A][ALLELE_A],
        b->counts[bundle_type_index][RTYPE_A][ALLELE_C],
        b->counts[bundle_type_index][RTYPE_A][ALLELE_G],
        b->counts[bundle_type_index][RTYPE_A][ALLELE_T],
        b->counts[bundle_type_index][RTYPE_A][ALLELE_DEL],
        b->counts[bundle_type_index][RTYPE_B][ALLELE_A],
        b->counts[bundle_type_index][RTYPE_B][ALLELE_C],
        b->counts[bundle_type_index][RTYPE_B][ALLELE_G],
        b->counts[bundle_type_index][RTYPE_B][ALLELE_T],
        b->counts[bundle_type_index][RTYPE_B][ALLELE_DEL]);
}

static inline std::string dsa_duplex_base_consensus_qualities(const duplex_base_t *b) {
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

static inline std::string bundle_closed_dsa_id_chunk(const duplex_tag_info *info, const uint8_t bundle_type) {
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

std::string bundle_closed_pair_to_dsa_row(const bundle_closed_pair_t *bp, const std::string pos_prefix, const duplex_base_t *base, const uint8_t bundle_type) {

    // TODO: check if bundle is set or are the defaults all right?
    const duplex_tag_info *info = &bp->duplex_tag_info;
    const bundle_closed_t *bulk = &bp->bundles[BULK_INDEX];
    const bundle_closed_t *duplex = &bp->bundles[DUPLEX_INDEX];

    // TODO: check NM and others for float vs. double rounding differences

    const std::string prefix = std::format("{}\t{}\t{}\t{}",
        pos_prefix,
        custom_round(bulk->asxs),
        custom_round(bulk->nm),
        dsa_duplex_base_counts(base, BULK_INDEX));

    std::stringstream s;
    s
    << prefix
    << bundle_closed_dsa_id_chunk(info, bundle_type)
    << std::format("{}\t{}\t{}\t",
        custom_round(duplex->asxs),
        custom_round(duplex->clip),
        0.1 * custom_round(duplex->nm * 10.0))
    << dsa_duplex_base_counts(base, DUPLEX_INDEX)
    << dsa_duplex_base_consensus_qualities(base)
    << custom_round(bulk->proper_pairs) << '\t'
    << custom_round(duplex->proper_pairs) << '\n';

    return s.str();
}
