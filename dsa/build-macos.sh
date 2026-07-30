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
    -L"$HOME/.homebrew/lib" \
    -I"$HOME/.homebrew/include" \
    src/ref.cc \
    src/mask.cc \
    src/mask_loader.cc \
    src/pileup_custom.cc \
    src/pileup.cc \
    src/pileup_batch.cc \
    src/dsa.cc \
    -lhts -ldeflate -lz -lpthread -lcurl -ldl -llzma -lbz2 -lm -lssl -lcrypto \
    ${LDFLAGS:-}

codesign -s - -f --entitlements /dev/stdin dsa <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>com.apple.security.get-task-allow</key>
    <true/>
</dict>
</plist>
EOF

set +e
