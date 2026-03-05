#include <catch2/catch_test_macros.hpp>

#include "range.h"

TEST_CASE("regular range properties") {
  range_t r{0, 10}; // 0-indexed half-open

  REQUIRE (r.start == 0);
  REQUIRE (r.end == 10);
  REQUIRE (range_is_valid(&r));
  REQUIRE (range_length(&r) == 10);
  REQUIRE (range_contains(&r, 0));
  REQUIRE (range_contains(&r, 9));
  REQUIRE (!range_contains(&r, 10));  // excludes end
}

TEST_CASE("invalid ranges") {
  range_t r1{1, 1};
  REQUIRE (!range_is_valid(&r1));

  range_t r2{20, 10};
  REQUIRE (!range_is_valid(&r2));

  range_t r3{-10, -5};
  REQUIRE (!range_is_valid(&r3));
}

TEST_CASE("range_clamp") {
  range_t r{0, 30};
  range_t c{10, 20};

  range_clamp(&r, &c);

  REQUIRE(r.start == 10);
  REQUIRE(r.end == 20);
}

TEST_CASE("range_triplet_grow") {
  // TODO
}
