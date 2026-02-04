#include <algorithm>  // fill
#include "pileup_batch.h"
#include "utils.h"
#include <iostream>
#include <unordered_map>
#include <assert.h>
#include <set>
#include "pileup_custom.h"
#include "options.h"
#include "bundle.h"

#define BULK_UNUSABLE (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FSUPPLEMENTARY | BAM_FDUP)

static inline bool is_valid_read(const bam1_t *b) {
    return (read_has_flag(b, BAM_FREAD1) ^ read_has_flag(b, BAM_FREAD2));
}

static inline bool is_usable_read(const aux_t *aux, const bam1_t *b) {
    return ((int)b->core.qual >= aux->min_mapQ);
}

static inline bool is_usable_bulk_read(const aux_t *aux, const bam1_t *b) {
    return is_usable_read(aux, b) && !read_has_flag(b, BULK_UNUSABLE);
}

static inline bool is_usable_duplex_read(const aux_t *aux, const bam1_t *b) {
    return is_usable_read(aux, b) && read_has_tag(b, "RB");
}

void PileupBatch::Update(const char *contig, const range_tid_t range, MaskLoader mls[2], Ref *ref) {
    this->range = range_tid_to_range(&range);
    assert(range_length(&this->range) > 0);

    this->contig = contig;
    this->tid = range.tid;

    // Load masks
    this->mask.Reset(this->range);
    assert(mask.CountBytesSet() == 0);

    uint64_t max_masked_positions = 0;
    for (int i = 0; i < 2; ++i) {
        max_masked_positions += mls[i].LoadMask(this->contig, this->range.start, this->range.end, this->mask);
    }
    std::cerr << std::format("Masked positions: {}\n", mask.CountBytesSet());
    assert(mask.CountBytesSet() <= max_masked_positions);

    // Load reference sequence
    const range_t ref_range = range_grow(&this->range);
    std::cerr << std::format("SLICE: {}:{}-{}\n", contig, range.start, range.end);
    std::cerr << std::format("REF: {}:{}-{}\n", contig, ref_range.start, ref_range.end);
    ref->Fetch(contig, ref_range);
}

const std::string PileupBatch::PositionString(const char *contig, const int pos, Ref *ref, const uint8_t mask_values[MASK_COUNT]) {
    const std::string_view ctx = ref->GetTripletAround(pos);

    std::stringstream ss;
    ss << contig;
    ss << "\t";
    ss << pos;
    ss << "\t";
    ss << pos + 1;
    ss << "\t";
    ss << ctx;
    ss << "\t";
    ss << static_cast<int>(mask_values[MASK_INDEX_SNP]);
    ss << "\t";
    ss << static_cast<int>(mask_values[MASK_INDEX_NOISE]);

    return ss.str();
}

void PileupBatch::Pileup(aux_t **data, Ref *ref, const Options *opts, GzipCompressor *compressor) {
    uint8_t mask_flag = 0;
    uint8_t mask_values[MASK_COUNT] = {0, 0};
    std::string posn;

    base_array_t base_buffer = {};
    if (base_info_array_reset(&base_buffer, 256)) {
        throw std::runtime_error("Failed to allocate base info array!");
    }

    range_tid_t r = {this->range.start, this->range.end, this->tid};
    aux_t *bulk_aux = data[BULK_INDEX];
    aux_t *duplex_aux = data[DUPLEX_INDEX];

    aux_set_iterator(bulk_aux, r);
    aux_set_iterator(duplex_aux, r);

    // Bulk pileup
    bam1_t *read = bam_init1();

    // Genomic position -> bundle indices
    std::map<int32_t, pos_stats_t> pos_bundles = {};
    std::unordered_map<uint64_t, bundle_closed_t> duplex_bundles = {};
    std::vector<std::string> bundle_id_decoder = {};

    std::string bundle_id;

    const uint64_t min_dplx_depth = static_cast<uint64_t>(opts->min_dplx_depth);

    // TODO: make global (or part of the Pileup object)
    probs_t probs = {};
    probs_init(&probs);

    // TODO: add processed read counters

    {
        std::unordered_map<uint64_t, bundle_open_t> open_bundles = {};

        {
            std::unordered_map<std::string, uint64_t> bundle_id_encoder = {};

            // A. Aggregate bulk
            std::cerr << "Aggregating bulk reads..." << std::endl;
            bulk_read_info_t bulk_read_info = {};
            while (1) {
                const int rc = aux_iter(bulk_aux, read);
                if (rc < 0) {
                    if (rc == -1) {
                        break;
                    } else {
                        throw std::runtime_error("Failed to read from bulk file!");
                    }
                }

                if (!is_valid_read(read)) {
                    throw std::runtime_error("Invalid read in bulk file!");
                }

                if (is_usable_bulk_read(bulk_aux, read)) {
                    int strand = get_strand_index(read);
                    if (strand == STRAND_INDEX_IGNORE) {
                        // TODO: verify this recapitulates the original behaviour!
                        continue;
                    }

                    if (base_array_update(&base_buffer, read, opts->min_base_quality)) {
                        throw std::runtime_error("Failed to perform pileup on bulk read!");
                    }

                    // Update allele counts
                    base_t *bi;
                    bulk_base_open_t *bbx;
                    bulk_read_info_init(&bulk_read_info, read);
                    for (uint64_t i = 0; i < base_buffer.count; ++i) {
                        bi = &base_buffer.bases[i];
                        bbx = &pos_bundles[bi->aln_pos].bulk_base;
                        if (range_tid_contains(&r, bi->aln_pos)) {
                            bulk_base_update(bbx, &bulk_read_info, bi, strand);
                        }
                    }
                }
            }

            // B. Aggregate duplex
            std::cerr << "Aggregating duplex reads..." << std::endl;
            bundle_open_t *duplex_bundle;
            while (1) {
                const int rc = aux_iter(duplex_aux, read);
                if (rc < 0) {
                    if (rc == -1) {
                        break;
                    } else {
                        throw std::runtime_error("Failed to read from duplex file!");
                    }
                }

                if (!is_valid_read(read)) {
                    throw std::runtime_error("Invalid read in duplex file!");
                }

                if (is_usable_duplex_read(duplex_aux, read)) {
                    uint64_t bundle_index;
                    bundle_id = get_duplex_id(read);
                    if (bundle_id_encoder.contains(bundle_id)) {
                        bundle_index = bundle_id_encoder[bundle_id];
                    } else {
                        // Initialise and assign a the next duplex index
                        bundle_index = bundle_id_encoder.size();
                        bundle_id_encoder[bundle_id] = bundle_index;
                        open_bundles[bundle_index].duplex_tag_info = duplex_tag_info_parse(bundle_id);
                    }

                    // Update bundle
                    duplex_bundle = &open_bundles[bundle_index];
                    const int r_type = bundle_open_duplex_update(duplex_bundle, read);

                    // NOTE: do not filter by quality for duplex bundles (?)!
                    if (base_array_update(&base_buffer, read, 0)) {
                        throw std::runtime_error("Failed to perform pileup on bulk read!");
                    }

                    const int strand = get_strand_index(read);
                    const int read_index = get_read_type_index(read);

                    // Update allele counts
                    base_t *bi;
                    duplex_base_t *dbx;
                    for (uint64_t i = 0; i < base_buffer.count; ++i) {
                        bi = &base_buffer.bases[i];
                        dbx = &pos_bundles[bi->aln_pos].duplex_bases[bundle_index];
                        if (
                            range_tid_contains(&r, bi->aln_pos) &&
                            duplex_tag_info_is_pos_in_template(
                                &duplex_bundle->duplex_tag_info,
                                bi->aln_pos + opts->offset)
                        ) {
                            dbx->counts[r_type][bi->base]++;

                            // Duplex-specific
                            if (strand != STRAND_INDEX_IGNORE) {
                                dbx->duplex_depth[strand][read_index]++;
                            }
                            probs_add_p_error(&probs, bi->qual, bi->base, dbx->duplex_consensus_quality_accum[r_type]);
                        }
                    }
                }
            }

            // C. Generate duplex index decoder
            bundle_id_decoder.resize(bundle_id_encoder.size());
            for (auto kvp : bundle_id_encoder) {
                // From str -> int to int -> str
                bundle_id_decoder[kvp.second] = kvp.first;
            }

        }

        // D. Close bundles
        std::cerr << "Finalising bundle stats..." << std::endl;
        bundle_open_t *duplex_bundle_open;
        bundle_closed_t *duplex_bundle;
        for (auto kvp : open_bundles) {
            const uint64_t bundle_index = kvp.first;
            duplex_bundle_open = &kvp.second;
            duplex_bundle = &duplex_bundles[bundle_index];
            bundle_id = bundle_id_decoder[bundle_index];

            bundle_closed_duplex_init(duplex_bundle, duplex_bundle_open);
        }
    }

    // E. Process bundle stats by position

    std::cerr << "Generating DSA table..." << std::endl;
    int32_t pos;
    // bundle_closed_t *bulk_bundle;
    // bundle_closed_t *duplex_bundle;
    uint8_t bundle_type;
    std::string pos_prefix;
    bundle_closed_t *duplex_bundle;
    duplex_base_t *duplex_base;
    bulk_base_closed_t bulk_base = {};
    pos_stats_t *pos_stats;

    uint64_t dsa_row_count = 0;
    for (auto pos_bundles_kvp : pos_bundles) {
        pos = pos_bundles_kvp.first;
        pos_stats = &pos_bundles_kvp.second;

        // Generate DSA table row prefix
        mask_flag = this->mask.GetFlag(pos);
        mask_values[MASK_INDEX_SNP] = flag_is_set(mask_flag, MASK_FLAG_SNP);
        mask_values[MASK_INDEX_NOISE] = flag_is_set(mask_flag, MASK_FLAG_NOISE);

        std::stringstream s;

        bulk_base_closed_init(&bulk_base, &pos_stats->bulk_base);
        pos_prefix = bulk_base_get_dsa_chunk(&bulk_base, PositionString(this->contig, pos, ref, mask_values));

        for (auto bundle_index_probs_kvp : pos_stats->duplex_bases) {
            const uint64_t bundle_index = bundle_index_probs_kvp.first;

            duplex_base = &bundle_index_probs_kvp.second;
            // BEWARE: the argument gets modified!
            duplex_base_finalise(duplex_base);

            duplex_bundle = &duplex_bundles[bundle_index];

            bundle_type = duplex_base_get_bundle_type(duplex_base, min_dplx_depth);
            if (bundle_type != 0) {
                dsa_push_row(s, duplex_bundle, duplex_base, &bulk_base, pos_prefix, bundle_type);
            }

            dsa_row_count++;
        }

        // Dump every position
        compressor->compress(s.str());
        compressor->write();

    }

    std::cerr << std::format("Generated {} rows", dsa_row_count) << std::endl;

}
