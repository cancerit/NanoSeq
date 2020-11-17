#!/bin/bash

########## LICENCE ##########
# Copyright (c) 2020 Genome Research Ltd.
# 
# Author: Cancer Genome Project <cgphelp@sanger.ac.uk>
# 
# This file is part of NanoSeq.
# 
# NanoSeq is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
###########################


get_distro () {
  EXT=""
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
  elif [[ $2 == *.zip* ]] ; then
    EXT="zip"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
  elif [[ $2 == *.tgz* ]] ; then
    EXT="tgz"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi
  rm -f $1.$EXT
  if hash curl 2>/dev/null; then
    curl --retry 10 -sS -o $1.$EXT -L $2
  else
    wget --tries=10 -nv -O $1.$EXT $2
  fi
}

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/ICGC"
  exit 0
fi

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

set -e
# cleanup inst_path
mkdir -p $INST_PATH
cd $INST_PATH
INST_PATH=`pwd`
mkdir -p $INST_PATH/bin
cd $INIT_DIR

export PATH="$INST_PATH/bin:$PATH"

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

echo -n "Get libdeflate ..."
if [ -e $SETUP_DIR/libdeflateGet.success ]; then
  echo " already staged ...";
else
  echo
  cd $SETUP_DIR
  get_distro "libdeflate" "https://github.com/ebiggers/libdeflate/archive/$VER_LIBDEFLATE.tar.gz"
  touch $SETUP_DIR/libdeflateGet.success
fi

echo -n "Building libdeflate ..."
if [ -e $SETUP_DIR/libdeflate.success ]; then
  echo " previously installed ...";
else
  echo
  cd $SETUP_DIR
  mkdir -p libdeflate
  tar --strip-components 1 -C libdeflate -zxf libdeflate.tar.gz
  cd libdeflate
  make -j$CPU CFLAGS="-fPIC -O3" libdeflate.a
  PREFIX=$INST_PATH make install
  cd $SETUP_DIR
  rm -r libdeflate.tar.gz
  touch $SETUP_DIR/libdeflate.success
fi

echo -n "Get htslib ..."
if [ -e $SETUP_DIR/htslibGet.success ]; then
  echo " already staged ...";
else
  echo
  cd $SETUP_DIR
  get_distro "htslib" "https://github.com/samtools/htslib/releases/download/$VER_HTSLIB/htslib-$VER_HTSLIB.tar.bz2"
  touch $SETUP_DIR/htslibGet.success
fi

echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo " previously installed ...";
else
  echo
  cd $SETUP_DIR
  mkdir -p htslib
  tar --strip-components 1 -C htslib -jxf htslib.tar.bz2
  cd htslib
  export CFLAGS="-I$INST_PATH/include"
  export LDFLAGS="-L$INST_PATH/lib"
  ./configure --enable-plugins  --enable-libcurl --with-libdeflate --prefix=$INST_PATH
  make -j$CPU
  make install
  mkdir $INST_PATH/include/cram
  cp ./cram/*.h $INST_PATH/include/cram/
  cp header.h $INST_PATH/include
  cd $SETUP_DIR
  rm -r htslib.tar.bz2
  touch $SETUP_DIR/htslib.success
fi

#Non-critical for the BotSeq pipeline execution
#(used to compute purity of the BotSeq BAM files)
echo -n "Downloading VerifyBamID ..."
if [ -e $SETUP_DIR/verifyBAMID.success ]; then
  echo " previously downloaded ...";
else
  echo
  cd $SETUP_DIR
  rm -rf verifyBamID
  get_distro "VerifyBamID" "https://github.com/Griffan/VerifyBamID/archive/$VER_VERIFYBAMID.tar.gz"
  mkdir -p verifyBamID
  tar --strip-components 1 -C verifyBamID -xzf VerifyBamID.tar.gz
  cd verifyBamID
  mkdir build
  cd build
  cmake ..
  make
  make test
  cp ../bin/VerifyBamID $INST_PATH/bin
  cd $SETUP_DIR
  rm -f VerifyBamID.tar.gz
  touch $SETUP_DIR/verifyBAMID.success
fi

echo -n "Downloading samtools ..."
if [ -e $SETUP_DIR/samtools.success ]; then
  echo " previously downloaded ...";
else
  echo
  cd $SETUP_DIR
  rm -rf samtools
  get_distro "samtools" "https://github.com/samtools/samtools/releases/download/$VER_SAMTOOLS/samtools-$VER_SAMTOOLS.tar.bz2"
  mkdir -p samtools
  tar --strip-components 1 -C samtools -xjf samtools.tar.bz2
  cd samtools
  cp *.h $INST_PATH/include
  ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
  make -j$CPU all all-htslib
  make install all all-htslib
  cd $SETUP_DIR
  rm -f samtools.tar.bz2
  touch $SETUP_DIR/samtools.success
fi

#Non-critical R functions (can be easily run manually)
#(create plots and summary stats)
echo -n "Installing R libraries..."
mkdir -p $R_LIBS
if [ -e $SETUP_DIR/Rlib.success ]; then
  echo " previously installed ...";
else
  echo
  Rscript $INIT_DIR/build/libInstall.R
  cd $SETUP_DIR
  touch $SETUP_DIR/Rlib.success
fi

echo -n "Building gzstream ..."
if [ -e $SETUP_DIR/gzstream.success ]; then
  echo " previously downloaded ...";
else
  echo
  cd $SETUP_DIR
  rm -rf gzstream
  get_distro "gzstream" "https://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz"
  mkdir -p gzstream
  tar --strip-components 1 -C gzstream -xzf gzstream.tgz
  cd gzstream
  cp *.h $INST_PATH/include
  make
  cp libgzstream.a $INST_PATH/lib
  cd $SETUP_DIR
  rm -f gzstream.tgz
  touch $SETUP_DIR/gzstream.success
fi

