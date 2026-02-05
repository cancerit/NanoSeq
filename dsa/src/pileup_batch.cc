#include "pileup_batch.h"
#include "bulk_bundle.h"
#include "dsa_serialiser.h"
#include "duplex_bundle.h"
#include "pos_stats.h"
#include "read_info.h"
#include "utils.h"
#include <iostream>
#include <unordered_map>
#include <assert.h>
#include <map>
#include "pileup_custom.h"
#include "options.h"
#include <sstream>
#include "duplex_tag_info.h"

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

typedef struct pos_stats_t {
    bulk_bundle_t bulk = {};
    std::map<uint64_t, duplex_bundle_t> duplexes = {};
} pos_stats_t;

typedef struct duplex_info_t {
	uint64_t index;
	duplex_tag_info_t tag;
} duplex_info_t;

void PileupBatch::Pileup(pileup_state_t *state) {
    const Options *opts = state->opts;

    uint8_t mask_flag = 0;
    uint8_t mask_values[MASK_COUNT] = {0, 0};
    std::string posn;

    base_array_t base_buffer = {};
    if (base_info_array_reset(&base_buffer, 256)) {
        throw std::runtime_error("Failed to allocate base info array!");
    }

    const range_tid_t r = {
        .start = this->range.start,
        .end = this->range.end,
        .tid = this->tid
    };

    if (aux_set_iterator(state->bulk_aux, r)) {
        throw std::runtime_error("Failed to create bulk iterator!");
    }
    if (aux_set_iterator(state->duplex_aux, r)) {
        throw std::runtime_error("Failed to create duplex iterator!");
    }

    // Bulk pileup
    bam1_t *read = bam_init1();

    // Genomic position -> bundle indices
    std::map<int32_t, pos_stats_t> pos_bundles = {};
    // std::unordered_map<uint64_t, bundle_closed_t> duplex_bundles = {};
    std::vector<duplex_tag_info_t> bundle_id_decoder = {};

    std::string bundle_id;

    const uint64_t min_dplx_depth = static_cast<uint64_t>(opts->min_dplx_depth);

    // TODO: make global (or part of the Pileup object)
    probs_t probs = {};
    probs_init(&probs);

    // TODO: add processed read counters

    {
        {
            std::unordered_map<std::string, duplex_info_t> bundle_id_encoder = {};

            // A. Aggregate duplex
            std::cerr << "Aggregating duplex reads... ";
            uint64_t duplex_total_reads = 0;
            uint64_t duplex_usable_reads = 0;
            uint64_t duplex_positions = 0;
            uint64_t duplex_positions_in_range = 0;
            duplex_bundle_t *duplex_bundle;
            read_info_t read_info = {};
            duplex_info_t *duplex_info;
            uint64_t duplex_index;
            while (1) {
                const int rc = aux_iter(state->duplex_aux, read);
                if (rc < 0) {
                    if (rc == -1) {
                        break;
                    } else {
                        throw std::runtime_error("Failed to read from duplex file!");
                    }
                }
                duplex_total_reads++;

                if (!is_valid_read(read)) {
                    throw std::runtime_error("Invalid read in duplex file!");
                }

                if (is_usable_duplex_read(state->duplex_aux, read)) {
                    duplex_usable_reads++;


                    bundle_id = get_duplex_id(read);
                    if (bundle_id_encoder.contains(bundle_id)) {
                    	duplex_info = &bundle_id_encoder[bundle_id];
                    } else {
                        // Initialise and assign a the next duplex index
                        duplex_index = bundle_id_encoder.size();
                        duplex_info = &bundle_id_encoder[bundle_id];
                        duplex_info->index = duplex_index;
                        duplex_info->tag = duplex_tag_info_parse(bundle_id);
                    }

                    // Update bundle
                    read_info_duplex_init(&read_info, read);

                    // NOTE: do not filter by quality for duplex bundles (?)!
                    if (base_array_update(&base_buffer, read, 0)) {
                        throw std::runtime_error("Failed to perform pileup on bulk read!");
                    }

                    // Update allele counts
                    base_t *bi;
                    for (uint64_t i = 0; i < base_buffer.count; ++i) {
                        duplex_positions++;

                        bi = &base_buffer.bases[i];
                        if (
                            range_tid_contains(&r, bi->aln_pos) &&
                            duplex_tag_info_is_pos_in_template(
                                &duplex_info->tag,
                                bi->aln_pos + opts->offset)
                        ) {
                            duplex_positions_in_range++;
                            duplex_bundle = &pos_bundles[bi->aln_pos].duplexes[duplex_info->index];
                            duplex_bundle_update(duplex_bundle, &read_info, &probs, bi);
                        }
                    }
                }
            }
            std::cerr << std::format(
                "{}/{} ({:.0f}%) usable reads, {}/{} ({:.0f}%) positions in range\n",
                duplex_usable_reads,
                duplex_total_reads,
                static_cast<double>(duplex_usable_reads) / static_cast<double>(duplex_total_reads) * 100.0,
                duplex_positions_in_range,
                duplex_positions,
                static_cast<double>(duplex_positions_in_range) / static_cast<double>(duplex_positions) * 100.0);

            // B. Aggregate bulk
            std::cerr << "Aggregating bulk reads... ";
            uint64_t bulk_total_reads = 0;
            uint64_t bulk_usable_reads = 0;
            uint64_t bulk_positions = 0;
            uint64_t bulk_positions_in_range = 0;
            while (1) {
                const int rc = aux_iter(state->bulk_aux, read);
                if (rc < 0) {
                    if (rc == -1) {
                        break;
                    } else {
                        throw std::runtime_error("Failed to read from bulk file!");
                    }
                }
                bulk_total_reads++;

                if (!is_valid_read(read)) {
                    throw std::runtime_error("Invalid read in bulk file!");
                }

                if (is_usable_bulk_read(state->bulk_aux, read)) {
                    bulk_usable_reads++;

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
                    bulk_bundle_t *bbx;
                    read_info_init(&read_info, read);
                    for (uint64_t i = 0; i < base_buffer.count; ++i) {
                        bulk_positions++;
                        bi = &base_buffer.bases[i];

                        if (pos_bundles.contains(bi->aln_pos)) {
                            bulk_positions_in_range++;
                        	bbx = &pos_bundles[bi->aln_pos].bulk;
                            bulk_bundle_update(bbx, &read_info, bi);
                        }
                    }
                }
            }
            std::cerr << std::format(
                "{}/{} ({:.0f}%) usable reads, {}/{} ({:.0f}%) positions in range\n",
                bulk_usable_reads,
                bulk_total_reads,
                static_cast<double>(bulk_usable_reads) / static_cast<double>(bulk_total_reads) * 100.0,
                bulk_positions_in_range,
                bulk_positions,
                static_cast<double>(bulk_positions_in_range) / static_cast<double>(bulk_positions) * 100.0);

            // C. Generate duplex index decoder
            bundle_id_decoder.resize(bundle_id_encoder.size());
            for (auto kvp : bundle_id_encoder) {
                // From str -> (int, tag) to int -> tag
                bundle_id_decoder[kvp.second.index] = kvp.second.tag;
            }

        }

        // D. Close bundles
        /*
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
        */
    }

    // E. Process bundle stats by position

    std::cerr << "Generating DSA table..." << std::endl;
    int32_t pos;
    // bundle_closed_t *bulk_bundle;
    // bundle_closed_t *duplex_bundle;
    uint8_t bundle_type;
    std::string pos_prefix;
    duplex_bundle_t *duplex_bundle;
    bulk_bundle_t *bulk_bundle;
    pos_stats_t *pos_stats;
    pos_final_stats_t duplex_stats = {};
    pos_final_stats_t bulk_stats = {};
    duplex_tag_info_t *duplex_tag_info;

    uint64_t dsa_row_count = 0;
    for (auto pos_bundles_kvp : pos_bundles) {
        pos = pos_bundles_kvp.first;
        pos_stats = &pos_bundles_kvp.second;

        // Generate DSA table row prefix
        mask_flag = this->mask.GetFlag(pos);
        mask_values[MASK_INDEX_SNP] = flag_is_set(mask_flag, MASK_FLAG_SNP);
        mask_values[MASK_INDEX_NOISE] = flag_is_set(mask_flag, MASK_FLAG_NOISE);

        std::stringstream s;

        bulk_bundle = &pos_stats->bulk;

        // TODO: check all attributes get overridden!
        bulk_bundle_finalise(bulk_bundle, &bulk_stats);
        pos_prefix = dsa_get_bulk_prefix(bulk_bundle, &bulk_stats, PositionString(this->contig, pos, state->ref, mask_values));

        for (auto bundle_index_probs_kvp : pos_stats->duplexes) {
            const uint64_t bundle_index = bundle_index_probs_kvp.first;
            duplex_bundle = &bundle_index_probs_kvp.second;

            // BEWARE: the argument gets modified!
            // TODO: check all attributes get overridden!
            duplex_bundle_finalise(duplex_bundle, &duplex_stats);

            // TODO: prune the pileup by bundle type before bulk is processed?
            bundle_type = duplex_base_get_bundle_type(duplex_bundle, min_dplx_depth);
            if (bundle_type != 0) {
            	duplex_tag_info = &bundle_id_decoder[bundle_index];
                dsa_push_row(s, duplex_tag_info, duplex_bundle, &duplex_stats, &bulk_stats, pos_prefix, bundle_type);
            }

            dsa_row_count++;
        }

        // Dump every position
        state->compressor->compress(s.str());
        state->compressor->write();

    }

    std::cerr << std::format("Generated {} rows", dsa_row_count) << std::endl;

    bam_destroy1(read);

}
