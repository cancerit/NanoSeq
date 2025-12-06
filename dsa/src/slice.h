#ifndef SLICE_H_
#define SLICE_H_

#include <vector>
#include <algorithm>  // fill
#include <assert.h>
#include "range.h"

template<typename T>
class Slice {
  private:
    range_t range;
    std::vector<T> values;
  public:
    void Update(const range_t range, const T value);
    void Reset(const range_t range);
    int32_t GetIndex(const int32_t pos);
    T Get(const int32_t pos);
    T *From(const int32_t pos);
};

template <typename T>
inline int32_t Slice<T>::GetIndex(const int32_t pos) {
  return pos - this->range.start;
}

template<typename T>
T Slice<T>::Get(const int32_t pos) {
  return this->values[this->GetIndex(pos)];
}

template<typename T>
T *Slice<T>::From(const int32_t pos) {
  return &this->values[this->GetIndex(pos)];
}

template<typename T>
void Slice<T>::Reset(const range_t range) {
  const int32_t m = range_length(&range);
  const int32_t n = range_length(&this->range);
  assert(this->values.size() == n);

  // Expand the mask (if necessary) and reset it to zero
  if (m > n) {
    this->values.resize(m);
  }
  std::fill(this->values.begin(), this->values.end(), 0);
}

template<typename T>
void Slice<T>::Update(const range_t range, const T value) {
  const int32_t a = range.start - this->range.start;
  if (range.end == range.start) {
    this->values[a] |= value;
  } else {
    for (int i = a; i <= range_length(&range); ++i) {
      this->values[i] |= value;
    }
  }
}

#endif
