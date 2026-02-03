#include "bundle.h"
#include "utils.h"
#include <format>

static inline double bundle_mean(const bundle_open_t *b, const uint64_t r_type, const int64_t x[2]) {
    return static_cast<double>(x[r_type]) / static_cast<double>(b->rtype_read_counts[r_type]);
}

static inline double bulk_base_mean(const bulk_base_open_t *b, const uint64_t r_type, const int64_t x[2]) {
    return static_cast<double>(x[r_type]) / static_cast<double>(b->read_counts[r_type]);
}

static inline void bundle_closed_set_common(bundle_closed_t *s, const bundle_open_t *b) {
    const bool has_a = b->rtype_read_counts[RTYPE_A] != 0;
    const bool has_b = b->rtype_read_counts[RTYPE_B] != 0;

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

void bundle_closed_duplex_init(bundle_closed_t *s, const bundle_open_t *b) {
    const bool has_a = b->rtype_read_counts[RTYPE_A] != 0;
    const bool has_b = b->rtype_read_counts[RTYPE_B] != 0;
    bundle_closed_set_common(s, b);
    s->duplex_tag_info = b->duplex_tag_info;

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

void bulk_base_closed_init(bulk_base_closed_t *s, const bulk_base_open_t *b) {
    std::memcpy(s->counts, b->counts, sizeof(s->counts));

    const bool has_a = b->ppair_accum[RTYPE_A] != 0;
    const bool has_b = b->ppair_accum[RTYPE_B] != 0;

    // Mean AS - XS and proper pair counts
    if (has_a && has_b) {

        s->asxs = std::min(
            bulk_base_mean(b, RTYPE_A, b->asxs_accum),
            bulk_base_mean(b, RTYPE_B, b->asxs_accum));

        s->proper_pairs = std::min(
            bulk_base_mean(b, RTYPE_A, b->ppair_accum),
            bulk_base_mean(b, RTYPE_B, b->ppair_accum));

    } else if (has_a) {
        s->asxs = bulk_base_mean(b, RTYPE_A, b->asxs_accum);
        s->proper_pairs = bulk_base_mean(b, RTYPE_A, b->ppair_accum);
    } else if (has_b) {
        s->asxs = bulk_base_mean(b, RTYPE_B, b->asxs_accum);
        s->proper_pairs = bulk_base_mean(b, RTYPE_B, b->ppair_accum);
    } else {
        s->asxs = 0.0;
        s->proper_pairs = 0.0;
    }

    // NM
    if (has_a && has_b) {
        // TODO: should this be max(a, b) instead (same as duplex NM)?
        s->nm = bulk_base_mean(b, RTYPE_A, b->nm_accum);
    } else if (has_a) {
        s->nm = bulk_base_mean(b, RTYPE_A, b->nm_accum);
    } else if (has_b) {
        s->nm = bulk_base_mean(b, RTYPE_B, b->nm_accum);
    } else {
        s->nm = 0.0f;
    }
}

/* SERIALISATION */

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

std::string bulk_base_get_dsa_chunk(bulk_base_closed_t *bulk, const std::string pos_prefix) {
    return std::format("{}\t{}\t{}\t{}",
        pos_prefix,
        custom_round(bulk->asxs),
        custom_round(bulk->nm),
        dsa_base_counts(bulk->counts));
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

void dsa_push_row(
    std::stringstream &s,
    const bundle_closed_t *duplex_bundle,
    const duplex_base_t *duplex_base,
    const bulk_base_closed_t *bulk_base,
    const std::string pos_prefix,
    const uint8_t bundle_type
) {
    s
    << pos_prefix
    << bundle_closed_dsa_id_chunk(&duplex_bundle->duplex_tag_info, bundle_type)
    << std::format("{}\t{}\t{}\t",
        custom_round(duplex_bundle->asxs),
        custom_round(duplex_bundle->clip),
        static_cast<float>(0.1 * custom_round(duplex_bundle->nm * 10.0)))
    << dsa_base_counts(duplex_base->counts)
    << dsa_duplex_base_consensus_qualities(duplex_base)
    << custom_round(bulk_base->proper_pairs) << '\t'
    << custom_round(duplex_bundle->proper_pairs) << '\n';
}
