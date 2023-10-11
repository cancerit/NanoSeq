#!/bin/bash

########## LICENCE ##########
# Copyright (c) 2022 Genome Research Ltd
# 
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
# 
# This file is part of NanoSeq.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
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

echo -n "Building libdeflate ..."
if [ -e $SETUP_DIR/libdeflate.success ]; then
  echo " previously built ...";
else
  echo
  cd $SETUP_DIR
  mkdir -p libdeflate
  get_distro "libdeflate" "https://github.com/ebiggers/libdeflate/archive/$VER_LIBDEFLATE.tar.gz"
  tar --strip-components 1 -C libdeflate -zxf libdeflate.tar.gz
  cd libdeflate
  mkdir build
  cmake -B build
  cmake --build build
  cmake --install build
  cmake --install build --prefix $INST_PATH
  cd $SETUP_DIR
  rm -r libdeflate.tar.gz
  touch $SETUP_DIR/libdeflate.success
fi

echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo " previously built ...";
else
  echo
  cd $SETUP_DIR
  mkdir -p htslib
  get_distro "htslib" "https://github.com/samtools/htslib/releases/download/$VER_HTSLIB/htslib-$VER_HTSLIB.tar.bz2"
  tar --strip-components 1 -C htslib -jxf htslib.tar.bz2
  cd htslib
  export CFLAGS="-I$INST_PATH/include"
  export LDFLAGS="-L$INST_PATH/lib"
  ./configure --disable-plugins  --enable-libcurl --with-libdeflate --prefix=$INST_PATH
  make -j$CPU
  make install
  mkdir $INST_PATH/include/cram
  cp ./cram/*.h $INST_PATH/include/cram/
  cp header.h $INST_PATH/include
  cd $SETUP_DIR
  rm -r htslib.tar.bz2
  unset CFLAGS
  unset LDFLAGS
  unset LIBS
  touch $SETUP_DIR/htslib.success
fi

echo -n "Building bcftools ..."
if [ -e $SETUP_DIR/bcftools.success ]; then
  echo " previously built ...";
else
  echo
  cd $SETUP_DIR
  rm -rf bcftools
  get_distro "bcftools" "https://github.com/samtools/bcftools/releases/download/$VER_BCFTOOLS/bcftools-$VER_BCFTOOLS.tar.bz2"
  mkdir -p bcftools
  tar --strip-components 1 -C bcftools -xjf bcftools.tar.bz2
  cd bcftools
  ./configure --enable-libgsl --enable-perl-filters --prefix=$INST_PATH
  make -j$CPU
  make install
  cd $SETUP_DIR
  rm -f bcftools.tar.bz2
  touch $SETUP_DIR/bcftools.success
fi

echo -n "Building samtools ..."
if [ -e $SETUP_DIR/samtools.success ]; then
  echo " previously built ...";
else
  echo
  cd $SETUP_DIR
  rm -rf samtools
  get_distro "samtools" "https://github.com/samtools/samtools/releases/download/$VER_SAMTOOLS/samtools-$VER_SAMTOOLS.tar.bz2"
  mkdir -p samtools
  tar --strip-components 1 -C samtools -xjf samtools.tar.bz2
  cd samtools
  ./configure --with-htslib=$SETUP_DIR/htslib --enable-plugins --enable-libcurl --prefix=$INST_PATH
  make -j$CPU
  make install
  cd $SETUP_DIR
  rm -f samtools.tar.bz2
  touch $SETUP_DIR/samtools.success
fi

#used for QC of NanoSeq BAM files
echo -n "Building VerifyBamID ..."
if [ -e $SETUP_DIR/verifyBAMID.success ]; then
  echo " previously built ...";
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
  cmake -DCMAKE_PREFIX_PATH=$INST_PATH ..
  make -j$CPU
  make test
  cp ../bin/VerifyBamID $INST_PATH/bin
  cd $SETUP_DIR
  rm -f VerifyBamID.tar.gz
  touch $SETUP_DIR/verifyBAMID.success
fi

echo -n "Installing R libraries..."
mkdir -p $R_LIBS
if [ -e $SETUP_DIR/Rlib.success ]; then
  echo " previously installed ...";
else
  echo
  export R_LIBS=$INST_PATH/R-lib
  export R_LIBS_USER=$INST_PATH/R-lib
  mkdir -p $R_LIBS_USER
  Rscript $INIT_DIR/build/libInstall.R $R_LIBS_USER
  cd $SETUP_DIR
  touch $SETUP_DIR/Rlib.success
fi

echo -n "Building gzstream ..."
if [ -e $SETUP_DIR/gzstream.success ]; then
  echo " previously built ...";
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

