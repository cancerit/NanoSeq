#pragma once

#include <format>
#include <cstddef>

template<std::size_t MaxLength>
class StaticStringBuilder {
    char buffer_[MaxLength];
    std::size_t length_ = 0;

public:
    const std::size_t& length = length_;

    const char *data() const { return buffer_; }
    void reset(std::size_t length = 0) { length_ = length; }

    template<typename... Args>
    StaticStringBuilder& append(std::format_string<Args...> fmt, Args&&... args) {
        auto result = std::format_to_n(buffer_ + length_, MaxLength - length_, fmt, std::forward<Args>(args)...);
        length_ += result.size;
        return *this;
    }

    StaticStringBuilder& rtrim(const std::size_t x) {
        if (x >= length_) {
            length_ = 0;
        } else {
            length_ -= x;
        }
    }

    std::string to_string() {
        return std::string(buffer_, length_);
    }
};
