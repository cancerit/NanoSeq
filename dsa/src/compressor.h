#ifndef GZIP_COMPRESSOR_H_
#define GZIP_COMPRESSOR_H_

#include "data_stats.h"

#include <libdeflate.h>

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

class GzipCompressor {
public:
    explicit GzipCompressor(std::filesystem::path output_path, std::string fn, int compression_level = 6)
        : final_path_(std::move(output_path.append(fn)))
        , temp_path_(final_path_.string() + ".part")
        , output_file_(temp_path_, std::ios::binary | std::ios::trunc)
        , compressor_(libdeflate_alloc_compressor(compression_level))
    {
        if (!output_file_) {
            throw std::runtime_error(
                std::format("Failed to open output file: {}", temp_path_.string()));
        }

        if (!compressor_) {
            throw std::runtime_error("Failed to allocate libdeflate compressor");
        }
    }

    ~GzipCompressor() {
        if (compressor_) {
            libdeflate_free_compressor(compressor_);
        }
    }

    // Non-copyable
    GzipCompressor(const GzipCompressor&) = delete;
    GzipCompressor& operator=(const GzipCompressor&) = delete;

    // Movable
    GzipCompressor(GzipCompressor&& other) noexcept
        : final_path_(std::move(other.final_path_))
        , temp_path_(std::move(other.temp_path_))
        , output_file_(std::move(other.output_file_))
        , compressor_(std::exchange(other.compressor_, nullptr))
        , uncompressed_stats_(std::move(other.uncompressed_stats_))
        , compressed_stats_(std::move(other.compressed_stats_))
        , compression_buffer_(std::move(other.compression_buffer_))
        , pending_buffer_(std::move(other.pending_buffer_))
        , pending_size_(other.pending_size_)
        , finalized_(other.finalized_)
    {}

    GzipCompressor& operator=(GzipCompressor&& other) noexcept {
        if (this != &other) {
            if (compressor_) {
                libdeflate_free_compressor(compressor_);
            }

            final_path_ = std::move(other.final_path_);
            temp_path_ = std::move(other.temp_path_);
            output_file_ = std::move(other.output_file_);
            compressor_ = std::exchange(other.compressor_, nullptr);
            uncompressed_stats_ = std::move(other.uncompressed_stats_);
            compressed_stats_ = std::move(other.compressed_stats_);
            compression_buffer_ = std::move(other.compression_buffer_);
            pending_buffer_ = std::move(other.pending_buffer_);
            pending_size_ = other.pending_size_;
            finalized_ = other.finalized_;
        }
        return *this;
    }

    void compress(const std::string& data) {
        if (finalized_) {
            throw std::logic_error("Cannot compress after finalization");
        }

        if (data.empty()) {
            return;
        }

        // Update uncompressed stats
        uncompressed_stats_.update(data);

        // Ensure compression buffer is large enough
        std::size_t required_size = libdeflate_gzip_compress_bound(compressor_, data.size());
        if (compression_buffer_.size() < required_size) {
            compression_buffer_.resize(required_size);
        }

        // Compress
        std::size_t compressed_size = libdeflate_gzip_compress(
            compressor_,
            data.data(),
            data.size(),
            compression_buffer_.data(),
            compression_buffer_.size()
        );

        if (compressed_size == 0) {
            throw std::runtime_error("Compression failed");
        }

        // Update compressed stats
        compressed_stats_.update(compression_buffer_.data(), compressed_size);

        // Append to pending buffer (grow only)
        std::size_t new_pending_size = pending_size_ + compressed_size;
        if (pending_buffer_.size() < new_pending_size) {
            pending_buffer_.resize(new_pending_size);
        }
        std::memcpy(pending_buffer_.data() + pending_size_,
                    compression_buffer_.data(),
                    compressed_size);
        pending_size_ = new_pending_size;
    }

    void write() {
        if (finalized_) {
            throw std::logic_error("Cannot write after finalization");
        }

        if (pending_size_ == 0) {
            return;
        }

        output_file_.write(
            reinterpret_cast<const char*>(pending_buffer_.data()),
            static_cast<std::streamsize>(pending_size_)
        );

        if (!output_file_) {
            throw std::runtime_error("Failed to write compressed data to file");
        }

        // Reset pending size but keep buffer allocated
        pending_size_ = 0;

        compression_buffer_.clear();
        pending_buffer_.clear();
    }

    void finalise() {
        if (finalized_) {
            return;
        }

        // NOTE: checking everything has been written
        assert(pending_size_ == 0);

        // Flush any pending data
        // write();

        finalized_ = true;

        // Release buffer memory
        compression_buffer_ = std::vector<std::uint8_t>();
        pending_buffer_ = std::vector<std::uint8_t>();

        // Close the output file
        output_file_.close();

        // Finalize stats
        uncompressed_stats_.finalize();
        compressed_stats_.finalize();

        // Rename .part file to final name
        std::error_code ec;
        std::filesystem::rename(temp_path_, final_path_, ec);
        if (ec) {
            throw std::runtime_error(
                std::format("Failed to rename {} to {}: {}",
                    temp_path_.string(), final_path_.string(), ec.message()));
        }

        // Write report.json alongside the output file
        auto report_path = final_path_.parent_path() / "report.json";
        if (final_path_.parent_path().empty()) {
            report_path = "report.json";
        }

        write_report(report_path);
    }

    // Accessors
    [[nodiscard]] const DataStats& uncompressed_stats() const noexcept {
        return uncompressed_stats_;
    }

    [[nodiscard]] const DataStats& compressed_stats() const noexcept {
        return compressed_stats_;
    }

    [[nodiscard]] double compression_ratio() const noexcept {
        if (uncompressed_stats_.size() == 0) return 0.0;
        return static_cast<double>(compressed_stats_.size()) /
               static_cast<double>(uncompressed_stats_.size());
    }

    [[nodiscard]] const std::filesystem::path& output_path() const noexcept {
        return final_path_;
    }

    [[nodiscard]] std::size_t compression_buffer_capacity() const noexcept {
        return compression_buffer_.capacity();
    }

    [[nodiscard]] std::size_t pending_buffer_capacity() const noexcept {
        return pending_buffer_.capacity();
    }

    [[nodiscard]] std::size_t pending_size() const noexcept {
        return pending_size_;
    }

    [[nodiscard]] bool is_finalized() const noexcept {
        return finalized_;
    }

private:
    std::filesystem::path final_path_;
    std::filesystem::path temp_path_;
    std::ofstream output_file_;

    libdeflate_compressor* compressor_ = nullptr;

    DataStats uncompressed_stats_;
    DataStats compressed_stats_;

    std::vector<std::uint8_t> compression_buffer_;  // Reused for each compress() call
    std::vector<std::uint8_t> pending_buffer_;      // Accumulates until write()
    std::size_t pending_size_ = 0;

    bool finalized_ = false;

    void write_report(const std::filesystem::path& report_path) const {
        std::ofstream report(report_path, std::ios::trunc);
        if (!report) {
            throw std::runtime_error(
                std::format("Failed to open report file: {}", report_path.string()));
        }

        report << std::format(R"({{
  "output_file": "{}",
  "uncompressed": {{
    "size_bytes": {},
    "md5": "{}"
  }},
  "compressed": {{
    "size_bytes": {},
    "md5": "{}"
  }},
  "compression_ratio": {:.4f}
}}
)",
            final_path_.filename().string(),
            uncompressed_stats_.size(),
            uncompressed_stats_.md5(),
            compressed_stats_.size(),
            compressed_stats_.md5(),
            compression_ratio()
        );
    }
};

#endif // GZIP_COMPRESSOR_H_
