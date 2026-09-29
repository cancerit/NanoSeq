/*########## LICENCE ##########
# Copyright (c) 2022, 2025, 2026 Genome Research Ltd
#
# Authors: Luca Barbon <lb29@sanger.ac.uk>, Alex Byrne <ab63@sanger.ac.uk>
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

#ifndef DUPLEX_TAG_INFO_H_
#define DUPLEX_TAG_INFO_H_

#include <charconv>
#include <stdexcept>
#include <string>

typedef struct duplex_tag_info_t {
    int32_t beg;
    int32_t end;
    std::string fwd_bc;
    std::string rev_bc;
} duplex_tag_info_t;

inline bool duplex_tag_info_is_pos_in_template(const duplex_tag_info_t *info, const int32_t pos) {
    // ASSUMPTION: offset correction (based on the convention used for the position)
    //  has been applied to the input position.
    return (pos >= info->beg) && (pos <= info->end);
}

// Format: token0,beg,end,fwd_bc,rev_bc
inline duplex_tag_info_t duplex_tag_info_parse(const std::string& s) {
    // Skip token0
    const size_t p0 = s.find(',');

    // beg
    const size_t p1 = s.find(',', p0 + 1);
    int32_t beg;
    auto [ptr1, ec1] = std::from_chars(s.data() + p0 + 1, s.data() + p1, beg);
    if (ec1 != std::errc{}) [[unlikely]] {
        throw std::runtime_error("Failed to parse beg in duplex tag");
    }

    // end
    const size_t p2 = s.find(',', p1 + 1);
    int32_t end;
    auto [ptr2, ec2] = std::from_chars(s.data() + p1 + 1, s.data() + p2, end);
    if (ec2 != std::errc{}) [[unlikely]] {
        throw std::runtime_error("Failed to parse end in duplex tag");
    }

    // fwd_bc
    const size_t p3 = s.find(',', p2 + 1);

    // rev_bc (remainder)
    return {
        .beg = beg,
        .end = end,
        .fwd_bc = s.substr(p2 + 1, p3 - p2 - 1),
        .rev_bc = s.substr(p3 + 1)
    };
}

#endif
