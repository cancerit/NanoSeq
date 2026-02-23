set -e

clang++ \
    -o dsa \
    -g \
    -O3 \
    -std=c++23 \
    ${CXXFLAGS:-} \
    -Wall -Wextra -Wpedantic -Wnull-dereference -Warray-bounds -Wformat=2 \
    -Wno-unused-private-field \
    -Wno-unused-function \
    src/ref.cc \
    src/mask.cc \
    src/mask_loader.cc \
    src/pileup_custom.cc \
    src/pileup.cc \
    src/pileup_batch.cc \
    src/dsa.cc \
    -lhts -ldeflate -lz -lpthread -lcurl -ldl -llzma -lbz2 -lm -lssl -lcrypto \
    ${LDFLAGS:-}

set +e
