#include <format>
#include <iostream>
#include "utils.h"
#include "mask_loader.h"

MaskLoader::MaskLoader(const uint8_t index) {
    if (index > 7) {
        throw std::runtime_error(std::format(
            "Invalid mask index {} (must be in [0, 7])!", index));
    }
    this->index = index;
    this->flag = 1 << this->index;
}

void MaskLoader::Init(const char *bed_fp) {
    this->bed_fp = bed_fp;

    if (bed_fp == nullptr || bed_fp[0] == '\0') {
        throw std::runtime_error("No BED file name provided!");
    }

    this->f = hts_open(bed_fp, "r");
    if (this->f == nullptr) {
        throw std::runtime_error(
            std::format("Error: failed to open {}\n", bed_fp));
    }

    this->tbx = tbx_index_load(bed_fp);
    if (this->tbx == nullptr) {
        throw std::runtime_error(
            std::format("Error: failed to open .tbi index of {}\n", bed_fp));
    }
}

void MaskLoader::LoadMask(const char *contig, const int start, const int end, Mask &mask) {
    assert(end - start > 0);
    const range_t t = {start, end};

    // char region[MAX_REGION_STR_LENGTH];
    // get_region(contig, start, end, region);
    // tbx_itr_queryi(tbx, tid, beg, end) hts_itr_query((tbx)->idx, (tid), (beg), (end), tbx_readrec)
    const int tid = tbx_name2id(tbx, contig);
    hts_itr_t *itr = tbx_itr_queryi(tbx, tid, start, end);
    // hts_itr_t *itr = tbx_itr_querys(this->tbx, region);
    range_t r;
    if (itr == nullptr) {
        // TODO: warn
        return;
    }

    kstring_t str = {};
    int32_t field_count = 0;
    int32_t *fields;
    uint64_t variant_count = 0;
    uint64_t position_count = 0;
    while (tbx_itr_next(this->f, this->tbx, itr, &str) >= 0) {
        fields = ksplit(&str, 0, &field_count);
        if (fields == nullptr || field_count < 3) {
            throw std::runtime_error(std::format(
                "Invalid entry in mask file!"));
        }
        r.start = std::stoi(&str.s[fields[1]]);
        r.end = std::stoi(&str.s[fields[2]]);

        range_clamp(&r, &t);

        variant_count++;
        position_count += range_length(&r);

        // std::cerr << contig << ":" << r.start << "-" << r.end << std::endl;
        mask.Update(r, this->flag);
    }
    tbx_itr_destroy(itr);

    std::cerr << std::format(
        "Loaded {} variants ({}/{} positions covered).\n",
        variant_count, position_count, range_length(&t));
}
