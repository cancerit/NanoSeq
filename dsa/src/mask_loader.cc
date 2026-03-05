#include <cstdint>
#include <format>
#include <iostream>
#include "mask_loader.h"
#include "range.h"

MaskLoader::MaskLoader(const uint8_t bit_index) {
    constexpr auto FLAG_WIDTH
        = std::numeric_limits<decltype(this->flag)>::digits;
    // bit_index is an index into a bitwise flag
    if (bit_index >= FLAG_WIDTH) {
        throw std::runtime_error(std::format(
            "Error: Invalid mask index {} (must be in [0, {}])!",
            index,
            FLAG_WIDTH - 1
        ));
    }
    this->index = bit_index;
    this->flag = static_cast<uint8_t>(1 << this->index);  // set single bit
}

// long term, string view may be easier than const char*
void MaskLoader::Init(const char *bed_fp) {
    if (bed_fp == nullptr || bed_fp[0] == '\0') {
        throw std::runtime_error("Error: No BED file name provided!");
    }

    this->mask_fp = bed_fp;

    this->f = hts_open(mask_fp, "r");
    if (this->f == nullptr) {
        throw std::runtime_error(
            std::format("Error: failed to open {}\n", mask_fp));
    }

    this->tbx = tbx_index_load(mask_fp);
    if (this->tbx == nullptr) {
        throw std::runtime_error(
            std::format("Error: failed to open .tbi index of {}\n", mask_fp));
    }
}

uint64_t MaskLoader::LoadMask(const char *contig, const int start, const int end, Mask &mask) {
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
        std::cerr << std::format(
            "Failed to load mask iterator for range {}:{}-{}\n",
            contig, start, end);
        return 0;
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
                "Invalid entry in mask {}!", this->mask_fp));
        }
        try {
            r.start = std::stoi(&str.s[fields[1]]);
            r.end = std::stoi(&str.s[fields[2]]);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                std::format(
                    "Could not convert coordinate in mask {}"
                    " to integer: {}",
                    this->mask_fp,
                    e.what()
                )
            );
        }

        if (!range_is_valid(&r)) {
            std::cerr << std::format(
                    "Warning: skipping invalid range [{}-{}] in mask {} -"
                    " ranges must be valid positive half open coordinates "
                    "(start >= 0 && end > start)",
                    r.start,
                    r.end,
                    this->mask_fp
            );
            continue;
        }

        range_clamp(&r, &t);

        variant_count++;
        position_count += static_cast<uint64_t>(range_length(&r));

        // std::cerr << contig << ":" << r.start << "-" << r.end << std::endl;
        mask.Update(r, this->flag);
    }
    tbx_itr_destroy(itr);

    std::cerr << std::format(
        "Loaded {} variants ({}/{} positions covered).\n",
        variant_count, position_count, range_length(&t));

    return position_count;
}
