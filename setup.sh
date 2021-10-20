#!/bin/bash

########## LICENCE ##########
# Copyright (c) 2020-2021 Genome Research Ltd
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


SOURCE_LIBDEFLATE="https://github.com/ebiggers/libdeflate/archive/v1.6.tar.gz"
SOURCE_HTSLIB="https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2"
SOURCE_GZSTREAM="https://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz"

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

# make sure that build is self contained
unset PERL5LIB
ARCHNAME=`perl -e 'use Config; print $Config{archname};'`
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"
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
  get_distro "libdeflate" $SOURCE_LIBDEFLATE
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
  get_distro "htslib" $SOURCE_HTSLIB
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
  #make SEQS_PER_SLICE 1000 for improved CRAM handling
  sed -ri 's/#define +SEQS_PER_SLICE .+/#define SEQS_PER_SLICE 1000/' ./cram/cram_structs.h
  export CFLAGS="-I$INST_PATH/include -D HAVE_LIBDEFLATE"
  export LDFLAGS="-L$INST_PATH/lib"
  #do this to force libdeflate to be statically liked, (for tabix and bgzip)
  #so that programs can be executed directly after a setup installation
  export HAVE_LIBDEFLATE=1
  export LIBS="-l:libdeflate.a"
  ./configure --disable-plugins  --enable-libcurl --without-libdeflate --prefix=$INST_PATH
  make -j$CPU
  make install
  mkdir -p $INST_PATH/include/cram
  cp ./cram/*.h $INST_PATH/include/cram/
  cp header.h $INST_PATH/include
  cd $SETUP_DIR
  rm -r htslib.tar.bz2
  unset HAVE_LIBDEFLATE
  unset CFLAGS
  unset LDFLAGS
  unset LIBS
  touch $SETUP_DIR/htslib.success
fi

echo -n "Building gzstream ..."
if [ -e $SETUP_DIR/gzstream.success ]; then
  echo " previously downloaded ...";
else
  echo
  cd $SETUP_DIR
  rm -rf gzstream
  get_distro "gzstream" $SOURCE_GZSTREAM
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


echo "Compiling code form this repository"
if [ -e $SETUP_DIR/botseq.success ]; then
  echo " previously compiled";
else
  export PREFIX="$INST_PATH"
  export LD_LIBRARY_PATH=$INST_PATH/lib:$LD_LIBRARY_PATH
  cd $INIT_DIR
  make
  make install
  make test
  touch $SETUP_DIR/botseq.success
fi

cp $INIT_DIR/python/* $INST_PATH/bin
cp $INIT_DIR/R/* $INST_PATH/bin/
cp $INIT_DIR/perl/* $INST_PATH/bin

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo
