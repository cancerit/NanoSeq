#include <algorithm>  // fill
#include "pileup_batch.h"
#include "utils.h"
#include <iostream>
#include <unordered_map>
#include <assert.h>
#include <set>
#include "pileup_custom.h"
#include "read_info.h"
#include "options.h"
#include "bundle.h"

#define BULK_UNUSABLE (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FSUPPLEMENTARY | BAM_FDUP)

static inline bool is_valid_read(const bam1_t *b) {
    return (read_has_flag(b, BAM_FREAD1) ^ read_has_flag(b, BAM_FREAD2));
}

static inline bool is_valid_bulk_read(const aux_t *aux, const bam1_t *b) {
    return ((int)b->core.qual >= aux->min_mapQ) && !read_has_flag(b, BULK_UNUSABLE);
}

static int RetrieveAlignments(void *data, bam1_t *b) {
    aux_t *aux = (aux_t *)data;

    if (!aux->iter) {
        return -1;
    }

    int ret;
    while (1) {
        // ret = sam_itr_next(aux->fp, aux->iter, b);
        // ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->head, b);
        ret = aux_iter(aux, b);

        if (ret < 0) {
            if (ret == -1) {
                fprintf(stderr, "TYPE=%d (%d:%d-%d) EOF after %llu iterations!\n", aux->duplex, aux->range.tid, aux->range.start, aux->range.end, aux->iterations);
                break;
            } else {
                std::stringstream er;
                er << "Error: failure while reading input BAM";
                er << std::endl;
                throw std::runtime_error(er.str());
            }
        }

        // TODO: put here all read checks that would lead to a critical error (spares branches)
        assert(read_has_flag(b, BAM_FREAD1) ^ read_has_flag(b, BAM_FREAD2));

        if (
            ((aux->duplex == 1) && (bam_aux_get(b, "RB") == NULL)) ||
            ((int)b->core.qual < aux->min_mapQ)
        ) {
            continue;
        }

        break;
    }
    return ret;
}

static int get_duplex_read(void *data, bam1_t *b) {
    aux_t *aux = (aux_t *)data;

    if (!aux->iter) {
        return -1;
    }

    int ret;
    while (1) {
        ret = aux_iter(aux, b);

        if (ret < 0) {
            if (ret == -1) {
                fprintf(stderr, "DUPLEX (%d:%d-%d) EOF after %llu iterations!\n", aux->range.tid, aux->range.start, aux->range.end, aux->iterations);
                break;
            } else {
                std::stringstream er;
                er << "Error: failure while reading input BAM";
                er << std::endl;
                throw std::runtime_error(er.str());
            }
        }

        // TODO: put here all read checks that would lead to a critical error (spares branches)
        assert(read_has_flag(b, BAM_FREAD1) ^ read_has_flag(b, BAM_FREAD2));

        if (
            ((aux->duplex == 1) && (bam_aux_get(b, "RB") == NULL)) ||
            ((int)b->core.qual < aux->min_mapQ)
        ) {
            continue;
        }

        break;
    }
    return ret;
}

static int get_bulk_read(void *data, bam1_t *b) {
    aux_t *aux = (aux_t *)data;

    if (!aux->iter) {
        return -1;
    }

    int ret;
    while (1) {
        ret = aux_iter(aux, b);

        if (ret < 0) {
            if (ret == -1) {
                fprintf(stderr, "BULK (%d:%d-%d) EOF after %llu iterations!\n", aux->range.tid, aux->range.start, aux->range.end, aux->iterations);
                break;
            } else {
                std::stringstream er;
                er << "Error: failure while reading input BAM";
                er << std::endl;
                throw std::runtime_error(er.str());
            }
        }

        // TODO: put here all read checks that would lead to a critical error (spares branches)
        assert(read_has_flag(b, BAM_FREAD1) ^ read_has_flag(b, BAM_FREAD2));

        if ((int)b->core.qual < aux->min_mapQ) {
            continue;
        }

        break;
    }
    return ret;
}

void PileupBatch::Update(const char *contig, const range_tid_t range, MaskLoader mls[2], Ref *ref) {
    this->range = range_tid_to_range(&range);
    assert(range_length(&this->range) > 0);

    this->contig = contig;
    this->tid = range.tid;

    // Load masks
    this->mask.Reset(this->range);
    for (int i = 0; i < 2; ++i) {
        mls[i].LoadMask(this->contig, this->range.start, this->range.end, this->mask);
    }
    std::cerr << std::format("Masked positions: {}\n", mask.CountBytesSet());

    // Load reference sequence
    const range_t ref_range = range_grow(&this->range);
    std::cerr << std::format("SLICE: {}:{}-{}\n", contig, range.start, range.end);
    std::cerr << std::format("REF: {}:{}-{}\n", contig, ref_range.start, ref_range.end);
    ref->Fetch(contig, ref_range);

    // DEBUG ONLY!
    /*
    std::cerr << "[" << ref->ToString() << "]" << std::endl;
    std::cerr << "<" << ref->GetTripletAround(range.start) << ">" << std::endl;
    std::cerr << "<" << ref->GetTripletAround(range.start + 1) << ">" << std::endl;
    std::cerr << "<" << ref->GetTripletAround(range.end) << ">" << std::endl;
    */
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

// void ReadBundler::UpdateBulkBundle(bundle *bndl, const bam_pileup1_t *p, int min_base_quality) {
/*
void bundle_bulk_update(bundle *bndl, const ReadStats *rs, const BaseInfoArray *ba) {
  // TODO: verify behaviour on invalid strand!
  const int strand = rs->strand;

  if (strand == STRAND_INDEX_IGNORE) {
    return;
  }

  BaseInfo *bi = NULL;
  for (uint64_t i = 0; i < ba->count; ++i) {
    bi = &ba->bases[i];

    bndl->counts[strand][bi->base]++;
    bndl->call[strand].push_back({bi->base, bi->qual});
    bndl->asxs[strand].push_back(ReadBundler::ASMinusXS(p->b));
    bndl->nmms[strand].push_back(ReadBundler::AuxTagToInt(p->b, "NM"));
    // bndl->ppair[strand].push_back(read_is_in_proper_pair(p->b));
    bndl->rtype_ppair_counts[strand] += read_is_in_proper_pair(p->b);
    bndl->rtype_read_counts[strand]++;
  }

  std::pair<int, int> bq = ReadBundler::BaseAndQual(p);

  }
}
*/

// void PileupBatch::Pileup(aux_t **data, Ref *ref, WriteOut *out) {
void PileupBatch::Pileup(aux_t **data, Ref *ref, const Options *opts, GzipCompressor *compressor) {
    // const Options *opts = out->opts;
    // int n_plp[BUNDLE_TYPES_COUNT];
    std::vector<const bam_pileup1_t *> plps[BUNDLE_TYPES_COUNT];
    // const bam_pileup1_t *plp[BUNDLE_TYPES_COUNT];
    uint8_t mask_flag = 0;
    uint8_t mask_values[MASK_COUNT] = {0, 0};
    std::string posn;
    // ReadBundler rb = {};
    // rb.Init();

    // std::map<int32_t, std::map<std::string, bundle>> ppp = {};

    base_array_t base_buffer = {};
    if (base_info_array_reset(&base_buffer, 256)) {
        throw std::runtime_error("Failed to allocate base info array!");
    }

    // bam_mplp_t mplp = bam_mplp_init(BAM_COUNT, RetrieveAlignments, reinterpret_cast<void **>(data));
    // bam_mplp_set_maxcnt(mplp, opts->max_plp_depth);

    range_tid_t r = {this->range.start, this->range.end, this->tid};
    aux_t *bulk_aux = data[BULK_INDEX];
    aux_t *duplex_aux = data[DUPLEX_INDEX];

    aux_set_iterator(bulk_aux, r);
    aux_set_iterator(duplex_aux, r);

    // Bulk pileup
    int rc;
    bam1_t *read = bam_init1();
    base_t *bi = NULL;

    // Genomic position -> bundle indices
    std::map<int32_t, std::map<uint64_t, duplex_base_t>> pos_bundles = {};
    std::unordered_map<uint64_t, bundle_closed_pair_t> closed_bundles = {};
    std::vector<std::string> bundle_id_decoder = {};

    std::string bundle_id;
    uint64_t bundle_index;
    bundle_closed_pair_t *bp;
    duplex_base_t *dbx;

    const uint64_t min_dplx_depth = static_cast<uint64_t>(opts->min_dplx_depth);

    {
        std::unordered_map<uint64_t, bundle_open_t[BUNDLE_TYPES_COUNT]> open_bundles = {};

        {
            std::unordered_map<std::string, uint64_t> bundle_id_encoder = {};

            // A. Aggregate bulk
            while (1) {
                rc = aux_iter(bulk_aux, read);
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

                if (is_valid_bulk_read(bulk_aux, read)) {
                    int strand = get_strand_index(read);
                    if (strand == STRAND_INDEX_IGNORE) {
                        // TODO: verify this recapitulates the original behaviour!
                        continue;
                    }

                    bundle_id = get_duplex_id(read);
                    if (bundle_id_encoder.contains(bundle_id)) {
                        bundle_index = bundle_id_encoder[bundle_id];
                    } else {
                        bundle_index = bundle_id_encoder.size();
                        bundle_id_encoder[bundle_id] = bundle_index;
                    }

                    // Update bundle
                    bundle_open_t *bundle = &open_bundles[bundle_index][BUNDLE_TYPE_BULK];
                    bundle_open_bulk_update(bundle, read, strand);

                    // NOTE: do not filter by quality for duplex bundles (?)!
                    if (base_array_update(&base_buffer, read, opts->min_base_quality)) {
                        throw std::runtime_error("Failed to perform pileup on bulk read!");
                    }

                    // Update allele counts
                    dbx = &pos_bundles[bi->aln_pos][bundle_index];
                    for (uint64_t i = 0; i < base_buffer.count; ++i) {
                        bi = &base_buffer.bases[i];
                        if (range_tid_contains(&r, bi->aln_pos)) {
                            dbx->counts[BUNDLE_TYPE_BULK][strand][bi->base]++;
                        }
                    }
                }
            }

            // B. Aggregate duplex
            while (1) {
                rc = aux_iter(bulk_aux, read);
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

                {
                    bundle_id = get_duplex_id(read);
                    if (bundle_id_encoder.contains(bundle_id)) {
                        bundle_index = bundle_id_encoder[bundle_id];
                    } else {
                        bundle_index = bundle_id_encoder.size();
                        bundle_id_encoder[bundle_id] = bundle_index;
                    }

                    // Update bundle
                    bundle_open_t *bundle = &open_bundles[bundle_index][DUPLEX_INDEX];
                    const int r_type = bundle_open_duplex_update(bundle, read);

                    // NOTE: do not filter by quality for duplex bundles (?)!
                    if (base_array_update(&base_buffer, read, 0)) {
                        throw std::runtime_error("Failed to perform pileup on bulk read!");
                    }

                    const int strand = get_strand_index(read);
                    const int read_index = get_read_type_index(read);

                    // Update allele counts
                    dbx = &pos_bundles[bi->aln_pos][bundle_index];
                    for (uint64_t i = 0; i < base_buffer.count; ++i) {
                        bi = &base_buffer.bases[i];
                        if (range_tid_contains(&r, bi->aln_pos)) {
                            dbx->counts[DUPLEX_INDEX][r_type][bi->base]++;

                            // Duplex-specific
                            dbx->duplex_depth[strand][read_index]++;
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
        for (auto kvp : open_bundles) {
            bp = &closed_bundles[bundle_index];

            bundle_id = bundle_id_decoder[bundle_index];
            bp->duplex_tag_info = duplex_tag_info_parse(bundle_id);

            bundle_closed_bulk_init(&bp->bundles[BULK_INDEX], &kvp.second[BULK_INDEX]);
            bundle_closed_duplex_init(&bp->bundles[DUPLEX_INDEX], &kvp.second[DUPLEX_INDEX]);
        }
    }

    // E. Process bundle stats by position
    int32_t pos;
    // bundle_closed_t *bulk_bundle;
    // bundle_closed_t *duplex_bundle;
    uint8_t bundle_type;
    std::string pos_prefix;
    std::string dsa_row;
    for (auto pos_bundles_kvp : pos_bundles) {
        pos = pos_bundles_kvp.first;

        // Generate DSA table row prefix
        mask_flag = this->mask.GetFlag(pos);
        mask_values[MASK_INDEX_SNP] = flag_is_set(mask_flag, MASK_FLAG_SNP);
        mask_values[MASK_INDEX_NOISE] = flag_is_set(mask_flag, MASK_FLAG_NOISE);
        pos_prefix = PositionString(contig, pos, ref, mask_values);

        for (auto bundle_index_probs_kvp : pos_bundles_kvp.second) {
            bundle_index = pos_bundles_kvp.first;

            dbx = &bundle_index_probs_kvp.second;
            // BEWARE: the argument gets modified!
            duplex_base_finalise(dbx);
            // bundle_closed_to_dsa_row(closed_bundles[bundle_index], pos, dbx);

            bp = &closed_bundles[bundle_index];
            // duplex_bundle = &bp->bundles[DUPLEX_INDEX];

            bundle_type = duplex_base_get_bundle_type(dbx, min_dplx_depth);
            if (bundle_type != 0) {
                // bulk_bundle = &bp->bundles[BULK_INDEX];
                dsa_row = bundle_closed_pair_to_dsa_row(bp, pos_prefix, dbx, bundle_type);

            }

        }

        // Dump every position
        compressor->compress(dsa_row);
        compressor->write();

    }

    // std::map<int32_t, std::map<uint64_t, duplex_base_t>> pos_bundles = {};
    /*
    int32_t pos;
    bundle_open_t *bundle_open;
    for (auto pos_bundles_kvp : pos_bundles) {
        pos = pos_bundles_kvp.first;
        for (auto bundle_index_probs_kvp : pos_bundles_kvp.second) {
            bundle_index = pos_bundles_kvp.first;
            bundle_open = &pos_bundles_kvp.second;
        }
    }
    */

    /*
    bam_plp_t bulk_plp = bam_plp_init(get_bulk_read, data[BULK_INDEX]);
    bam_plp_set_maxcnt(bulk_plp, opts->max_plp_depth);
    bam_plp_t duplex_plp = bam_plp_init(get_duplex_read, data[DUPLEX_INDEX]);
    bam_plp_set_maxcnt(duplex_plp, opts->max_plp_depth);

    // while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {
    while (1) {
        plp[DUPLEX_INDEX] = bam_plp_auto(duplex_plp, &tid, &pos, &n_plp[DUPLEX_INDEX]);
        plp[BULK_INDEX] = bam_plp_auto(bulk_plp, &tid, &pos, &n_plp[BULK_INDEX]);

        plps[BUNDLE_TYPE_BULK].clear();
        plps[BUNDLE_TYPE_DUPLEX].clear();

        // TODO: verify end inclusiveness convention!
        //  Originally: ((pos >= opts->beg) && (pos <= opts->end))
        if (pos < this->range.start) {
            continue;
        } else if (pos >= this->range.end) {
            break;
        }

        // Pileup
        for (int i = 0; i < BUNDLE_TYPES_COUNT; i++) {
            for (int j = 0; j < n_plp[i]; ++j) {
                plps[i].push_back(plp[i] + j);
            }
        }

        // Bundle reads
        bundles dplx = rb.DplxBundles(pos, opts->offset, opts->min_dplx_depth, plps[BUNDLE_TYPE_DUPLEX]);
        // fprintf(stderr, "%d duplex bundles.\n", dplx.size());
        if (dplx.size() == 0) {
            continue;
        }
        bundle bulk = rb.BulkBundle(plps[BUNDLE_TYPE_BULK], opts->min_base_quality);
        // fprintf(stderr, "%d bulk calls.\n", bulk.call->size());


        // Generate DSA table row prefix
        mask_flag = this->mask.GetFlag(pos);
        mask_values[MASK_INDEX_SNP] = flag_is_set(mask_flag, MASK_FLAG_SNP);
        mask_values[MASK_INDEX_NOISE] = flag_is_set(mask_flag, MASK_FLAG_NOISE);

        // TODO: avoid string reallocation!
        posn = PositionString(contig, pos, ref, mask_values);

        // Push DSA table rows to compressor
        out->WriteRows(bulk, dplx, posn);
    }

    bam_plp_destroy(bulk_plp);
    bam_plp_destroy(duplex_plp);
    // bam_mplp_destroy(mplp);
    */

}
