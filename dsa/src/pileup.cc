/*########## LICENCE ##########
# Copyright (c) 2022, 2025 Genome Research Ltd
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

#include <format>
#include "pileup.h"
#include "pileup_batch.h"
#include "utils.h"
#include "ref.h"

static int RetrieveAlignments(void *data, bam1_t *b) {
    aux_t *aux = (aux_t *)data;
    int ret;
    while (1) {
        ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->head, b);

        if (ret == -1) {
            break;
        }

        if (ret < -1) {
            std::stringstream er;
            er << "Error: failure while reading input BAM";
            er << std::endl;
            throw std::runtime_error(er.str());
        }

        // TODO: put here all read checks that would lead to a critical error (spares branches)
        assert (read_has_flag(b, BAM_FREAD1) ^ read_has_flag(b, BAM_FREAD2));

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

std::vector<std::string> tokenize(std::string str, char delimiter) {
    std::istringstream iss(str);
    std::vector<std::string> tokens;
    std::string token;
    while (std::getline(iss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

static inline bool find_pg_id(const std::string s, const char *pattern) {
    // TODO: verify...
    return s.rfind(pattern, 0) == 0;
}

bool BamIsCorrectlyPreprocessed(bam_hdr_t *head, const int bundle_type) {
    bool op1 = false;
    bool op2 = false;
    bool op3 = false;
    bool op4 = false;
    std::vector<std::string> tokens = tokenize(head->text, '\n');
    for (int j = 0; j < tokens.size(); j++) {
        if (tokens[j].rfind("@PG", 0) == 0) {
            std::vector<std::string> subtokens = tokenize(tokens[j], '\t');
            for (int k = 0; k < subtokens.size(); k++) {
                if (subtokens[k].rfind("ID", 0) == 0) {
                    // account for suffixes on program identifiers
                    op1 |= find_pg_id(subtokens[k], "ID:bamsormadup");
                    op2 |= find_pg_id(subtokens[k], "ID:bammarkduplicatesopt");
                    op3 |= find_pg_id(subtokens[k], "ID:bamaddreadbundles");
                    op4 |= find_pg_id(subtokens[k], "ID:randomreadinbundle");
                }
            }
        }
    }
    switch (bundle_type) {
    case BUNDLE_TYPE_BULK:
        return !op3 || (op3 && op4);
    case BUNDLE_TYPE_DUPLEX:
        return op1 && op2 && op3;
    default:
        throw std::runtime_error("Invalid bundle type index!");
    }
}

Pileup::Pileup() : opts(nullptr), fai(nullptr) {
}

int Pileup::GetTID(const char *contig) {
    const int bulk_tid = sam_hdr_name2tid(this->data[BULK_INDEX]->head, contig);
    const int duplex_tid = sam_hdr_name2tid(this->data[DUPLEX_INDEX]->head, contig);
    assert(bulk_tid == duplex_tid);
    return bulk_tid;
}

const char *Pileup::GetContig(const int32_t tid) {
    return sam_hdr_tid2name(this->data[BULK_INDEX]->head, tid);
}

void Pileup::Initiate(Options *opts) {
    this->opts = opts;
    // this->out.opts = opts;

    // TODO: reintroduce write permission test?

    // Open masks (BED files)
    // Initialise mask flags
    static_assert(MASK_INDEX_SNP == 0);
    static_assert(MASK_INDEX_NOISE == 1);
    static_assert(MASK_COUNT == 2);

    for (int i = 0; i < MASK_COUNT; ++i) {
        std::cerr << std::format("Loading mask {} from {}...\n", i, this->opts->beds[i]);
        this->masks[i] = MaskLoader(i);
        this->masks[i].Init(this->opts->beds[i]);
    }

    // Load FAI
    // TODO: remove (now in Ref)
    this->fai = fai_load(this->opts->fasta);
    if (this->fai == nullptr) {
        std::stringstream er;
        er << "Error: failed to open index of ";
        er << this->opts->fasta;
        er << std::endl;
        throw std::runtime_error(er.str());
    }

    int tid = -1;
    this->data = reinterpret_cast<aux_t **>(calloc(BUNDLE_TYPES_COUNT, sizeof(aux_t *)));
    if (!this->data) {
        throw std::runtime_error("Failed to allocate!");
    }

    for (int i = 0; i < BUNDLE_TYPES_COUNT; ++i) {
        this->data[i] = reinterpret_cast<aux_t *>(calloc(1, sizeof(aux_t)));
        if (!this->data[i]) {
            throw std::runtime_error("Failed to allocate!");
        }

        this->data[i]->fp = hts_open(this->opts->bams[i], "r");
        if (i == BUNDLE_TYPE_DUPLEX) {
            this->data[i]->min_mapQ = this->opts->min_mapQ;
            this->data[i]->duplex = 1;
        } else {
            this->data[i]->min_mapQ = 0;
            this->data[i]->duplex = 0;
        }

        if (this->data[i]->fp == NULL) {
            std::stringstream er;
            er << "Error: failed to open ";
            er << this->opts->bams[i];
            er << std::endl;
            throw std::runtime_error(er.str());
        }
        this->data[i]->head = NULL;
        this->data[i]->head = sam_hdr_read(this->data[i]->fp);
        if (this->data[i]->head == NULL) {
            std::stringstream er;
            er << "Error: failed to read the header of ";
            er << this->opts->bams[i];
            er << std::endl;
            throw std::runtime_error(er.str());
        }

        // allow to skip tests
        if (this->opts->doTests) {
            if (!BamIsCorrectlyPreprocessed(this->data[i]->head, i)) {
                std::stringstream er;
                er << "Error : bam ";
                er << this->opts->bams[i];
                er << " is not properly preprocessed.";
                er << std::endl;
                throw std::runtime_error(er.str());
            }
        }

        this->indices[i] = sam_index_load(this->data[i]->fp, this->opts->bams[i]);
        if (this->indices[i] == NULL) {
            std::stringstream er;
            er << "Error: failed to load the index of ";
            er << this->opts->bams[i];
            er << std::endl;
            throw std::runtime_error(er.str());
        }
    }

    LoadRanges();

    if (this->opts->doTests) {
        // Check that the headers of both BAMs match each other
        const int n_targets_bulk = sam_hdr_nref(this->data[BULK_INDEX]->head);
        const int n_targets_duplex = sam_hdr_nref(this->data[DUPLEX_INDEX]->head);
        if (n_targets_bulk != n_targets_duplex) {
            std::stringstream er;
            er << "Error : number of chromosomes in bulk and duplex don't match (" << n_targets_bulk << ":" << n_targets_duplex << ")";
            er << std::endl;
            throw std::runtime_error(er.str());
        }
        for (int i = 0; i < n_targets_bulk; i++) {
            if (strcmp(sam_hdr_tid2name(this->data[0]->head, i), sam_hdr_tid2name(this->data[1]->head, i))) {
                std::stringstream er;
                er << "Error : order of chromosomes in bulk and duplex BAMs don't match";
                er << std::endl;
                throw std::runtime_error(er.str());
            }
        }
        // Check BAM contig names against the reference
        for (int i = 0; i < n_targets_bulk; i++) {
            if (!faidx_has_seq(this->fai, sam_hdr_tid2name(this->data[0]->head, i))) {
                std::stringstream er;
                er << "Error: BAM file chromosome " << sam_hdr_tid2name(this->data[0]->head, i) << " doesn't match any reference chromosome";
                er << std::endl;
                throw std::runtime_error(er.str());
            }
        }
        // Check that the BAM chomosome lenghts match the reference
        for (int i = 0; i < n_targets_bulk; i++) {
            if (sam_hdr_tid2len(this->data[0]->head, i) != faidx_seq_len(this->fai, sam_hdr_tid2name(this->data[0]->head, i))) {
                std::stringstream er;
                er << "Error: BAM file chromosome length for " << sam_hdr_tid2name(this->data[0]->head, i) << " doesn't match reference chromosome length";
                er << std::endl;
                throw std::runtime_error(er.str());
            }
        }
    }

    this->mplp = bam_mplp_init(BAM_COUNT, RetrieveAlignments, reinterpret_cast<void **>(this->data));
    bam_mplp_set_maxcnt(this->mplp, this->opts->max_plp_depth);
}

std::string Pileup::Header() {
    std::stringstream ss;
    ss << "#\n";
    ss << "# DESCRIPTION OF FIELDS\n";
    ss << "# 00  chrom                ";
    ss << "name of the chromosome or scaffold\n";
    ss << "# 01  chromBeg           ";
    ss << "beginning position of the feature in the chromosome or scaffold\n";
    ss << "# 02  chromEnd             ";
    ss << "ending position of the feature in the chromosome or scaffold\n";
    ss << "# 03  context              ";
    ss << "trinucleotide context\n";
    ss << "# 04  commonSNP            ";
    ss << "overlaps with dbSNP common 146\n";
    ss << "# 05  shearwater           ";
    ss << "overlaps with shearwater mask\n";
    ss << "# 06  bulkASXS             ";
    ss << "minimum of bulk forward and reverse mean AS-XS values\n";
    ss << "# 07  bulkNM               ";
    ss << "maximum of bulk forward and reverse mean NM values\n";
    ss << "# 08  bulkForwardA         ";
    ss << "count of bulk forward A bases\n";
    ss << "# 09  bulkForwardC         ";
    ss << "count of bulk forward C bases\n";
    ss << "# 10  bulkForwardG         ";
    ss << "count of bulk forward G bases\n";
    ss << "# 11  bulkForwardT         ";
    ss << "count of bulk forward T bases\n";
    ss << "# 12  bulkForwardIndel     ";
    ss << "count of bulk forward indels\n";
    ss << "# 13  bulkReverseA         ";
    ss << "count of bulk reverse A bases\n";
    ss << "# 14  bulkReverseC         ";
    ss << "count of bulk reverse C bases\n";
    ss << "# 15  bulkReverseG         ";
    ss << "count of bulk reverse G bases\n";
    ss << "# 16  bulkReverseT         ";
    ss << "count of bulk reverse T bases\n";
    ss << "# 17  bulkReverseIndel     ";
    ss << "count of bulk reverse indels\n";
    ss << "# 18  dplxBreakpointBeg    ";
    ss << "read-bundle beginning position\n";
    ss << "# 19  dplxBreakpointEnd    ";
    ss << "read-bundle ending position\n";
    ss << "# 20  dplxBarcode          ";
    ss << "read-bundle barcode sorted by (forward, reverse) strand\n";
    ss << "# 21  dplxBundleType       ";
    ss << "read-bundle type (1 = forward duplex, 2 = reverse duplex,";
    ss << " 3 = overlapping forward/reverse duplex)\n";
    ss << "# 22  dplxASXS             ";
    ss << "minimum of read-bundle f1r2 and f2r1 mean AS-XS values\n";
    ss << "# 23  dplxCLIP             ";
    ss << "proportion of read-bundle paired-end reads with a 5' clip\n";
    ss << "# 24  dplxNM               ";
    ss << "maximum of read-bundle f1r2 and f2r1 mean NM values\n";
    ss << "# 25  dplxf1r2A            ";
    ss << "count of read-bundle f1r2 A bases\n";
    ss << "# 26  dplxf1r2C            ";
    ss << "count of read-bundle f1r2 C bases\n";
    ss << "# 27  dplxf1r2G            ";
    ss << "count of read-bundle f1r2 G bases\n";
    ss << "# 28  dplxf1r2T            ";
    ss << "count of read-bundle f1r2 T bases\n";
    ss << "# 29  dplxf1r2Indel        ";
    ss << "count of read-bundle f1r2 indels\n";
    ss << "# 30  dplxf2r1A            ";
    ss << "count of read-bundle f2r1 A bases\n";
    ss << "# 31  dplxf2r1C            ";
    ss << "count of read-bundle f2r1 C bases\n";
    ss << "# 32  dplxf2r1G            ";
    ss << "count of read-bundle f2r1 G bases\n";
    ss << "# 33  dplxf2r1T            ";
    ss << "count of read-bundle f2r1 T bases\n";
    ss << "# 34  dplxf2r1Indel        ";
    ss << "count of read-bundle f2r1 indels\n";
    ss << "# 35  dplxCQf1r2A          ";
    ss << "consensus Phred quality value for read-bundle f1r2 A bases\n";
    ss << "# 36  dplxCQf1r2C          ";
    ss << "consensus Phred quality value for read-bundle f1r2 C bases\n";
    ss << "# 37  dplxCQf1r2G          ";
    ss << "consensus Phred quality value for read-bundle f1r2 G bases\n";
    ss << "# 38  dplxCQf1r2T          ";
    ss << "consensus Phred quality value for read-bundle f1r2 T bases\n";
    ss << "# 39  dplxCQf2r1A          ";
    ss << "consensus Phred quality value for read-bundle f2r1 A bases\n";
    ss << "# 40  dplxCQf2r1C          ";
    ss << "consensus Phred quality value for read-bundle f2r1 C bases\n";
    ss << "# 41  dplxCQf2r1G          ";
    ss << "consensus Phred quality value for read-bundle f2r1 G bases\n";
    ss << "# 42  dplxCQf2r1T          ";
    ss << "consensus Phred quality value for read-bundle f2r1 T bases\n";
    ss << "# 43  bulkProperPair       ";
    ss << "proportion of bulk reads that are properly-paired\n";
    ss << "# 44  dplxProperPair       ";
    ss << "proportion of dplx reads that are properly-paired\n";
    ss << "#";
    ss << "\n# PARAMETER VALUES";
    ss << "\n# Bulk BAM file:                         ";
    ss << this->opts->bams[0];
    ss << "\n# Duplex BAM file:                       ";
    ss << this->opts->bams[1];
    ss << "\n# SNP BED file:                          ";
    ss << this->opts->beds[0];
    ss << "\n# Mask BED file:                         ";
    ss << this->opts->beds[1];
    ss << "\n# Reference FASTA file:                  ";
    ss << this->opts->fasta;
    ss << "\n# Coordinate offset:                     ";
    ss << this->opts->offset;
    ss << "\n# Maximum pileup depth:                  ";
    ss << this->opts->max_plp_depth;
    ss << "\n# Minimum read-bundle depth per strand:  ";
    ss << this->opts->min_dplx_depth;
    ss << "\n#\n# TAB SEPARATED HEADER\n";
    ss << "#chrom\tchromStart\tchromEnd\tcontext\t";
    ss << "commonSNP\tshearwater\t";
    ss << "bulkASXS\tbulkNM\t";
    ss << "bulkForwardA\tbulkForwardC\tbulkForwardG\tbulkForwardT\t";
    ss << "bulkForwardIndel\t";
    ss << "bulkReverseA\tbulkReverseC\tbulkReverseG\tbulkReverseT\t";
    ss << "bulkReverseIndel\t";
    ss << "dplxBreakpointBeg\tdplxBreakpointEnd\tdplxBarcode\tdplxBundleType\t";
    ss << "dplxASXS\tdplxCLIP\tdplxNM\t";
    ss << "dplxf1r2A\tdplxf1r2C\tdplxf1r2G\tdplxf1r2T\tdplxf1r2Indel\t";
    ss << "dplxf2r1A\tdplxf2r1C\tdplxf2r1G\tdplxf2r1T\tdplxf2r1Indel\t";
    ss << "dplxCQf1r2A\tdplxCQf1r2C\tdplxCQf1r2G\tdplxCQf1r2T\t";
    ss << "dplxCQf2r1A\tdplxCQf2r1C\tdplxCQf2r1G\tdplxCQf2r1T\t";
    ss << "bulkProperPair\tdplxProperPair";
    return ss.str();
}

std::string Pileup::PositionString(const char *contig, const int pos, const uint8_t mask_values[MASK_COUNT]) {

    // TODO: pass as argument or store as attribute!
    Ref ref;
    std::string_view ctx = ref.GetTripletAround(pos);

    std::stringstream ss;
    ss << contig;
    ss << "\t";
    ss << pos;
    ss << "\t";
    ss << pos + 1;
    ss << "\t";
    ss << ctx;
    ss << "\t";
    ss << mask_values[MASK_INDEX_SNP];
    ss << "\t";
    ss << mask_values[MASK_INDEX_NOISE];
    ss << "\t";

    return ss.str();
}

static void init_iterator(hts_itr_t **it, const hts_idx_t *idx, const range_tid_t *r) {
    if (*it != NULL) {
        sam_itr_destroy(*it);
    }

    *it = sam_itr_queryi(idx, r->tid, r->start, r->end + 1);
    if (*it == NULL) {
        std::stringstream er;
        er << "Error: failed to parse region";
        er << std::endl;
        throw std::runtime_error(er.str());
    }
}

void Pileup::InitIterators(const range_tid_t *r) {
    for (int i = 0; i < BUNDLE_TYPES_COUNT; ++i) {
        init_iterator(&this->data[i]->iter, this->indices[i], r);
    }
}

void log_skip_contig_range(const char *contig, const range_t range) {
    std::cerr << std::format(
        "Contig '{}' not found, skipping range {}:{}-{}!\n",
        contig, contig, range.start, range.end);
}

void Pileup::MultiplePileupInRange(const char *contig, const range_t range) {
    const int tid = GetTID(contig);
    if (tid < 0) {
        log_skip_contig_range(contig, range);
        return;
    }

    // Create iterators
    // InitIterators(&r);

    /*
    {
        int pos;
        int n_plp[BUNDLE_TYPES_COUNT];
        std::vector<const bam_pileup1_t *> plps[BUNDLE_TYPES_COUNT];
        const bam_pileup1_t **plp;
        while (bam_mplp_auto(this->mplp, &tid, &pos, n_plp, plp) > 0) {
            if ((pos >= opts->beg) && (pos <= opts->end)) {
                // pileup
                for (int i = 0; i < BUNDLE_TYPES_COUNT; i++) {
                    for (int j = 0; j < n_plp[i]; ++j) {
                        plps[i].push_back(this->plp[i] + j);
                    }
                }
                // bundle reads
                std::unique_ptr<ReadBundler> rb(new ReadBundler());
                bundles dplx = rb->DplxBundles(pos, this->opts->offset, this->opts->min_dplx_depth, plps[BUNDLE_TYPE_DUPLEX]);
                if (dplx.size() == 0) {
                    continue;
                }
                bundle bulk = rb->BulkBundle(plps[BUNDLE_TYPE_BULK], this->opts->min_base_quality);
                // output
                std::string posn = Pileup::PositionString(contig, pos);
                out->WriteRows(bulk, dplx, posn, this->gzout, this->opts->out2stdout);
            }
            if (pos > opts->end) {
                break;
            }
            plps[BUNDLE_TYPE_BULK].clear();
            plps[BUNDLE_TYPE_DUPLEX].clear();
        }
        free(plp);
    }
    */
}


void Pileup::MultiplePileup() {
    std::cerr << "Output directory: " << opts->oname << std::endl;
    std::unique_ptr<WriteOut> out(new WriteOut(opts));

    PileupBatch batch;
    for (auto r : this->ranges) {
        std::cerr << std::format("<TID:{}>:{}-{}\n", r.tid, r.start, r.end);
        batch.Update(GetContig(r.tid), {r.start, r.end}, masks);
    }

    // TODO: handle the empty output case better
    out->Finalise();

    return;

    //
    /*

    */

    // TODO: iterate over input ranges
    // TODO: convert from cgranges type


    out->Finalise();

    fai_destroy(this->fai);
    bam_mplp_destroy(this->mplp);
    for (int i = 0; i < BUNDLE_TYPES_COUNT; ++i) {
        sam_close(this->data[i]->fp);
        if (this->data[i]->iter) {
            hts_itr_destroy(this->data[i]->iter);
        }
        sam_hdr_destroy(this->data[i]->head);
        free(this->data[i]);
    }
    free(this->data);
}

void Pileup::LoadRanges() {

    const char *fp = this->opts->ranges_bed;
    htsFile *f = hts_open(fp, "r");
    if (fp == nullptr) {
        std::runtime_error(std::format(
            "Failed to open {}!", fp));
    }

    range_tid_t r = {};
    kstring_t ks = {0, 0, nullptr};
    std:: string contig_s;
    const char *contig;
    while (hts_getline(f, '\n', &ks) >= 0) {
        std::string line(ks.s, ks.l);
        std::stringstream ss(line);
        if (!(ss >> contig_s >> r.start >> r.end)) {
            std::runtime_error(std::format(
                "Failed to parse line in {}!", fp));
        }
        if (r.end <= r.start) {
            std::runtime_error(std::format(
                "Invalid range {}:{}-{} (end must be greater than start) in {}!",
                contig, r.start, r.end, fp));
        }
        contig = contig_s.c_str();
        r.tid = GetTID(contig);
        if (r.tid >= 0) {
            ranges.push_back(r);
        } else {
            log_skip_contig_range(contig, {r.start, r.end});
        }
    }

    free(ks.s);
    hts_close(f);
}
