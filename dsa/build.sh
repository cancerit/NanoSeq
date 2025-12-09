g++ \
    -o dsa \
    -g \
    -std=c++23 -Wall \
    -Wno-unused-private-field \
    -Wno-unused-function \
    -L"$HOME/.homebrew/lib" \
    -I"$HOME/.homebrew/include" \
    -L"external/gzstream" \
    -I"external/gzstream" \
    -I"external/bedtk" \
    src/ref.cc \
    src/mask.cc \
    src/mask_loader.cc \
    src/pileup.cc \
    src/pileup_batch.cc \
    src/writeout.cc \
    src/read_bundler.cc \
    src/dsa.cc \
    -lhts -ldeflate -lgzstream -lz -lpthread -lcurl -ldl -llzma -lbz2 -lm -lssl -lcrypto

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
