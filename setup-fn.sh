set -ue

detect_cpu () {
  if [ -f /proc/cpuinfo ]; then
    CPU=$(grep -c ^processor /proc/cpuinfo)
    if [ "$CPU" -gt "6" ]; then
      CPU=6
    fi
  else
    CPU=1
  fi
}

require_install_path () {
  # extra args beyond the install path are passed straight through to cmake
  if [ "$#" -lt "1" ] ; then
    echo "Please provide an installation path such as /opt/ICGC"
    exit 0
  fi
}

prep_dirs () {
  # $1 = requested install path; sets INST_PATH (canonicalized) and SETUP_DIR
  # requires INIT_DIR to already be set
  mkdir -p "$1"
  cd "$1"
  INST_PATH=$(pwd)
  mkdir -p "$INST_PATH/bin"
  cd "$INIT_DIR"
  SETUP_DIR="$INIT_DIR/install_tmp"
  mkdir -p "$SETUP_DIR"
}

install_repo_scripts () {
  cp "$INIT_DIR"/python/* "$INST_PATH/bin"
  cp "$INIT_DIR"/R/* "$INST_PATH/bin/"
  cp "$INIT_DIR"/perl/* "$INST_PATH/bin"
  chmod a+x "$INST_PATH"/bin/*
}

build_repo () {
  echo "Compiling code form this repository"
  if [ -e "$SETUP_DIR"/botseq.success ]; then
    echo " previously compiled";
  else
    cd "$INIT_DIR"
    cmake -S . -B build "$@" -DCMAKE_INSTALL_PREFIX="$INST_PATH" -DNANOSEQ_INSTALL=ON
    cmake --build build -j"$CPU"
    cmake --install build
    ctest --test-dir build --output-on-failure
    touch "$SETUP_DIR"/botseq.success
  fi
}
