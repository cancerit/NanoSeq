g++ \
    -o dsa \
    -std=c++23 -Wall \
    -Wno-unused-private-field \
    -Wno-unused-function \
    -L"$HOME/.homebrew/lib" \
    -I"$HOME/.homebrew/include" \
    -L"external/gzstream" \
    -I"external/gzstream" \
    -I"external/bedtk" \
    src/pileup.cc \
    external/bedtk/cgranges.o \
    src/bedtk_lite.cc \
    src/bed_reader.cc \
    src/writeout.cc \
    src/read_bundler.cc \
    src/dsa.cc \
    -lhts -ldeflate -lgzstream -lz -lpthread -lcurl -ldl -llzma -lbz2 -lm -lcrypto
