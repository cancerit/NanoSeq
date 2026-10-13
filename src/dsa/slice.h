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

#ifndef SLICE_H_
#define SLICE_H_

#include <assert.h>
#include <cstring>  // memset
#include <string>
#include <cstdlib>
#include <stdexcept>

#include "range.h"

template<typename T>
class Slice {
  private:
    range_t range;
    size_t capacity;
    T *values;

  public:
    Slice() : range({0, 0}), capacity(0), values(nullptr) {};
    bool IsNull();
    uint64_t CountBytesSet();
    void Update(const range_t range, const T value);
    void Set(const range_t range, T *new_values);
    void Reset(const range_t range, const bool zero);
    const std::string ToString();
    int32_t GetLength();
    int32_t GetIndex(const int32_t pos);
    T Get(const int32_t pos);
    T *From(const int32_t pos);
    T *Data();
};

template <typename T>
inline bool Slice<T>::IsNull() {
  return this->values == nullptr;
}

template <typename T>
inline uint64_t Slice<T>::CountBytesSet() {
  const uint64_t length = GetLength();
  uint64_t count = 0;
  for (uint64_t i = 0; i < length; ++i) {
    count += (values[i] == 0);
    /*
    if (values[i]) {
      std::cout << std::format("{}", values[i]);
    } else {
      std::cout << std::format(" ");
    }
    */
  }
  // std::cout << std::endl;
  return length - count;
}

template <typename T>
inline void Slice<T>::Set(const range_t r, T *new_values) {
  this->range = r;
  if (this->values != nullptr) {
    free(this->values);
  }
  this->values = new_values;
}

template <typename T>
inline int32_t Slice<T>::GetLength() {
  return range_length(&this->range);
}

template <typename T>
const std::string Slice<T>::ToString() {
  return std::string(this->values, this->GetLength());
}

template <typename T>
inline int32_t Slice<T>::GetIndex(const int32_t pos) {
  return pos - this->range.start;
}

template<typename T>
inline T Slice<T>::Get(const int32_t pos) {
  return this->values[this->GetIndex(pos)];
}

template<typename T>
T *Slice<T>::From(const int32_t pos) {
  return &this->values[this->GetIndex(pos)];
}

template<typename T>
inline T *Slice<T>::Data() {
  return this->values;
}

template<typename T>
void Slice<T>::Reset(const range_t r, const bool zero) {
  const int32_t new_range_length = range_length(&r);
  assert(new_range_length > 0);
  this->range = r;

  /*
  if (this->values) {
    free(this->values);
  }
  this->values = (T*)calloc(new_range_length + 1, sizeof(T));
  return;
  */

  const size_t n = this->capacity;
  const size_t m = static_cast<size_t>(new_range_length + 1);

  if (m > n) {
    this->capacity = m;
  }

  const size_t new_size = this->capacity * sizeof(T);

  if (this->values == nullptr) {
    this->values = static_cast<T*>(malloc(new_size));
  } else if (m > n) {
    this->values = static_cast<T*>(realloc(this->values, new_size));
  }

  if (this->values == nullptr) {
    throw std::runtime_error("Failed to allocate slice!");
  }

  if (zero) {
    memset(this->values, 0, new_size);
  }

}

template<typename T>
void Slice<T>::Update(const range_t r, const T value) {
  assert(range_is_valid(&r));
  const int32_t a = r.start - this->range.start;
  assert(a >= 0 && a < range_length(&this->range));
  for (int i = a; i < a + range_length(&r); ++i) {
    this->values[i] |= value;
  }
}

#endif
