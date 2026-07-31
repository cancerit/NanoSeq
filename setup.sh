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

set -uex

. ./setup-fn.sh  # source functions

require_install_path "$@"

detect_cpu
echo "Max compilation CPUs set to $CPU"

INST_PATH=$1

# get current directory
INIT_DIR=$(pwd)

prep_dirs "$INST_PATH"

# make sure that build is self contained
unset PERL5LIB
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"
export PATH="$INST_PATH/bin:$PATH"

# any args beyond the install path are passed straight through to cmake, e.g.
# ./setup.sh path_to_install -DNANOSEQ_FETCH_DEPS=OFF -DCMAKE_PREFIX_PATH=/existing/prefix
build_repo "${@:2}"

install_repo_scripts

# cleanup all intermediates
rm -rf "$SETUP_DIR"

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo
