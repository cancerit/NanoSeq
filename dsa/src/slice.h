#ifndef SLICE_H_
#define SLICE_H_

#include <vector>
#include <algorithm>  // fill
#include <assert.h>
#include "range.h"
#include <iostream>

template<typename T>
class Slice {
  private:
    range_t range;
    T *values;

  public:
    Slice() : range({0, 0}), values(nullptr) {};
    bool IsNull();
    uint64_t CountBytesSet();
    void Update(const range_t range, const T value);
    void Set(const range_t range, T *new_values);
    void Reset(const range_t range, const bool zero);
    int32_t GetLength();
    int32_t GetIndex(const int32_t pos);
    T Get(const int32_t pos);
    T *From(const int32_t pos);
};

template <typename T>
inline bool Slice<T>::IsNull() {
  return this->values == nullptr;
}

template <typename T>
inline uint64_t Slice<T>::CountBytesSet() {
  const uint64_t length = GetLength();
  uint64_t count = 0;
  for (int i = 0; i < length; ++i) {
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
inline void Slice<T>::Set(const range_t range, T *new_values) {
  this->range = range;
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
void Slice<T>::Reset(const range_t range, const bool zero) {
  const int32_t m = range_length(&range);
  assert(m > 0);
  this->range = range;

  if (this->values == nullptr) {

    this->values = (T*)calloc(m, sizeof(T));

  } else {

    const int32_t n = range_length(&this->range);

    // Expand the mask (if necessary) and reset it to zero
    if (m > n) {
      this->values = (T*)realloc(this->values, m);
      if (this->values == nullptr) {
        throw std::runtime_error("Failed to reallocate slice!");
      }
    }

    if (zero) {
      memset(this->values, 0, m);
    }

  }
}

template<typename T>
void Slice<T>::Update(const range_t range, const T value) {
  const int32_t a = range.start - this->range.start;
  assert(a >= 0 && a < range_length(&this->range));
  // std::cerr << std::format("{}-{}\t[{}] = {}\n", range.start, range.end, a, value);
  if (range.end == range.start) {
    this->values[a] |= value;
  } else {
    for (int i = a; i < a + range_length(&range); ++i) {
      this->values[i] |= value;
    }
  }
  // assert(CountBytesSet() >= range_length(&range));
}

#endif
