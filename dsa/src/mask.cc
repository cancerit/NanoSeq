#include "mask.h"

void Mask::Update(const range_t range, const uint8_t flag) {
  this->set_byte_count++;
  this->mask.Update(range, flag);
}

void Mask::Reset(const range_t range) {
  this->set_byte_count = 0;
  this->mask.Reset(range, true);
}

uint8_t Mask::GetFlag(const int32_t pos) {
  return this->mask.Get(pos);
}
