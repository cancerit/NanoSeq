#!/usr/bin/env python3

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


import sys
import argparse
import math
import statistics
import numpy
import os
import glob
import gzip

#simple script to compute statistics from cov files. The coverage values in these files are computed
#in windows of 100bp so they differer from single base coverage statistics.

parser = argparse.ArgumentParser(description="Compute various statistics from coverage (cov) files")
#arguments for all subcommands
parser.add_argument('--cov',   action='store', default='./tmpNanoSeq/cov', help='path of the cov directory (./tmpNanoSeq/cov)')

args = parser.parse_args()
tmpDir = args.cov

#check that coverage calculation completed correctly
if (not os.path.isfile(tmpDir+'/nfiles') ):
  sys.exit(tmpDir+'/nfiles not found!\n')
else :
  with open(tmpDir+'/nfiles') as iofile :
    nfiles = int(iofile.readline())

for i in range( nfiles ) :
  if ( len(glob.glob(tmpDir+"/%s.done"%(i+1) )) != 1 ) :
    sys.exit("\ncov job %s did not complete correctly\n"%(i+1))
  if ( len(glob.glob(tmpDir+"/%s.cov.bed.gz"%(i+1))) != 1 ) :
    sys.exit("\ncov job %s did not complete correctly\n"%(i+1))


print("\nParsing coverage files\n")
coverageL = []
for i in range(nfiles) :
  with gzip.open( tmpDir+"/%s.cov.bed.gz"%(i+1),'rt') as iofile :
    for iline in iofile :
      ichr = str(iline.split('\t')[0])
      ib = int(iline.split('\t')[1])
      ie = int(iline.split('\t')[2])
      cc = int(iline.split('\t')[3])
      if ( cc == 0 ) : continue
      coverageL.append( cc )
coverage = numpy.sort(numpy.array(coverageL))
print("\nCompleted parsing coverage files\n")
print("min: %s"%(coverage.min())) 
print("max: %s"%(coverage.max())) 
print("mean: %s"%(numpy.mean(coverage) ))
print("std: %s"%( numpy.std(coverage) ))
print("median: %s"%(numpy.median(coverage) ))
print("quartiles: %s"%(numpy.quantile(coverage,[0.25,0.5,0.75]) ))
print("last 90ile: %s"%(numpy.quantile(coverage,.90) ))
print("last 98ile: %s"%(numpy.quantile(coverage,.98) ))
print("last 99ile: %s"%(numpy.quantile(coverage,.99) ))

