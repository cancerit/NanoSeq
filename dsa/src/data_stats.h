/*########## LICENCE ##########
# Copyright (c) 2022, 2025, 2026 Genome Research Ltd
#
# Authors: Luca Barbon <lb29@sanger.ac.uk>
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

#ifndef DATA_STATS_H_
#define DATA_STATS_H_

#include <openssl/evp.h>

#include <cstdint>
#include <format>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

class DataStats {
public:
    DataStats()
        : md5_ctx_(EVP_MD_CTX_new())
    {
        if (!md5_ctx_) {
            throw std::runtime_error("Failed to allocate MD5 context");
        }
        EVP_DigestInit_ex(md5_ctx_, EVP_md5(), nullptr);
    }

    ~DataStats() {
        if (md5_ctx_) {
            EVP_MD_CTX_free(md5_ctx_);
        }
    }

    // Non-copyable
    DataStats(const DataStats&) = delete;
    DataStats& operator=(const DataStats&) = delete;

    // Movable
    DataStats(DataStats&& other) noexcept
        : md5_ctx_(std::exchange(other.md5_ctx_, nullptr))
        , total_size_(other.total_size_)
        , finalized_(other.finalized_)
        , md5_hex_(std::move(other.md5_hex_))
    {}

    DataStats& operator=(DataStats&& other) noexcept {
        if (this != &other) {
            if (md5_ctx_) {
                EVP_MD_CTX_free(md5_ctx_);
            }
            md5_ctx_ = std::exchange(other.md5_ctx_, nullptr);
            total_size_ = other.total_size_;
            finalized_ = other.finalized_;
            md5_hex_ = std::move(other.md5_hex_);
        }
        return *this;
    }

    void update(const void* data, std::size_t size) {
        if (finalized_) {
            throw std::logic_error("Cannot update finalized DataStats");
        }
        if (size == 0) {
            return;
        }
        total_size_ += size;
        EVP_DigestUpdate(md5_ctx_, data, size);
    }

    void update(std::span<const std::uint8_t> data) {
        update(data.data(), data.size());
    }

    void update(const std::string& data) {
        update(data.data(), data.size());
    }

    void finalize() {
        if (finalized_) {
            return;
        }
        finalized_ = true;

        unsigned char digest[EVP_MAX_MD_SIZE];
        unsigned int digest_len = 0;
        EVP_DigestFinal_ex(md5_ctx_, digest, &digest_len);

        md5_hex_.reserve(digest_len * 2);
        for (unsigned int i = 0; i < digest_len; ++i) {
            md5_hex_ += std::format("{:02x}", digest[i]);
        }
    }

    [[nodiscard]] std::uint64_t size() const noexcept {
        return total_size_;
    }

    [[nodiscard]] const std::string& md5() const {
        if (!finalized_) {
            throw std::logic_error("MD5 not available until finalized");
        }
        return md5_hex_;
    }

    [[nodiscard]] bool is_finalized() const noexcept {
        return finalized_;
    }

private:
    EVP_MD_CTX* md5_ctx_ = nullptr;
    std::uint64_t total_size_ = 0;
    bool finalized_ = false;
    std::string md5_hex_;
};

#endif // DATA_STATS_H_
