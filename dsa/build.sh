g++ \
    -o dsa \
    -std=c++0x -Wall \
    -L"$HOME/.homebrew/lib" \
    -I"$HOME/.homebrew/include" \
    -L"external/gzstream" \
    -I"external/gzstream" \
    src/pileup.cc \
    src/bed_reader.cc \
    src/writeout.cc \
    src/read_bundler.cc \
    src/dsa.cc \
    -lhts -ldeflate -lgzstream -lz -lpthread -lcurl -ldl -llzma -lbz2 -lm -lcrypto
