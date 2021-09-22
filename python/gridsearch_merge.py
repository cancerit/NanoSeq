#!/usr/bin/env python
# coding=utf-8

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

from __future__ import division
import argparse
import os
from collections import defaultdict, Counter


class MergeResults(object):

  """Merge results from multiple BED files and output metrics files
  """

  def __init__(self, dirname):
    """Initialize.
    Args:
      dirname: directory containing variantcaller results files
    """
    self.dirname    = dirname
    self.clipped    = { 0 : 'no.filter', 1 : '<=0.1', 2 : '==0' }
    self.asxs       = { 0 : 'no.filter', 1 : '>=0', 2 : '>=50', 3 : '>=100' }
    self.nm         = { 0 : 'no.filter', 1 : '<=5', 2 : '<=3', 3: '<=2' }
    self.matched    = { 0 : 'no.filter', 1 : '>=5', 2 : '>=10' }
    self.dplx       = { 0 : '>=2', 1 : '>=3'}
    self.consensus  = { 0 : 'no.filter', 1 : '>=30', 2 : '>=60', 3 : '>=90'}
    self.indels     = { 0 : 'no.filter', 1 : '<=0.05', 2 : '==0' }
    self.fiveprime  = { 0 : 'no.filter', 1 : '>=5', 2 : '>=10' }
    self.threeprime = { 0 : 'no.filter', 1 : '>=5', 2 : '>=10' }
    self.ppairs     = { 0 : 'no.filter', 1 : '>=0.95', 2 : '==1' }
    self.isvar      = { 0 : 'ref', 1 : 'var' }
    self.main()


  def main(self):
    """Main function.
    """
    bedfiles   = self.bedfiles()
    counts     = self.count_variants(bedfiles)
    cumulative = self.cumulative(counts)
    self.writeout(cumulative)


  def bedfiles(self):
    """Returns list of BED filenames for all .bed files in the directory
    """
    fns = []
    for fn in os.listdir(self.dirname):
      if fn.endswith('.bed'):
        fns.append('%s/%s' % (self.dirname, fn))
    return fns


  def count_variants(self, bedfiles):
    """Read files from individual beds and output into combined metric files.
    Args:
      bedfiles : list of bedfile names
    """
    counts = Counter()
    for bed in bedfiles:
      for row in open(bed, 'rU'):
        row        = row.strip().split(',')
        clipped    = int(row[0])
        asxs       = int(row[1])
        nm         = int(row[2])
        matched    = int(row[3])
        dplx       = int(row[4])
        consensus  = int(row[5])
        indels     = int(row[6])
        fiveprime  = int(row[7])
        threeprime = int(row[8])
        ppairs     = int(row[9])
        isvar      = int(row[10])
        count      = int(row[11])
        counts[(clipped, asxs, nm, matched, dplx, consensus, indels,
                fiveprime, threeprime, ppairs, isvar)] += count
    return counts


  def cumulative(self, counts):
    """Sum calls at or above a given threshold value.
    Args:
      counts : dict of key = (clipped, asxs, nm, matched, consensus, rpbundle,
        indels, fiveprime, ppairs, isvar) and values = counts
    Returns:
      cumulative : nested dict key1 = (clipped, asxs, nm, matched, dplx, consensus,
        rpbundle, indels, fiveprime, ppairs), key2 = 1 if variant otherwise 0, 
        values, values = counts
    """
    cumulative = defaultdict(Counter)
    for i0 in self.clipped.keys():
      for i1 in self.asxs.keys():
        for i2 in self.nm.keys():
          for i3 in self.matched.keys():
            for i4 in self.dplx.keys():
              for i5 in self.consensus.keys():
                for i6 in self.indels.keys():
                  for i7 in self.fiveprime.keys():
                    for i8 in self.threeprime.keys():
                      for i9 in self.ppairs.keys():
                        for (j0, j1, j2, j3, j4, j5, j6, j7, j8, j9, isvar),count in counts.iteritems():
                          if (j0 >= i0 and
                              j1 >= i1 and
                              j2 >= i2 and
                              j3 >= i3 and
                              j4 >= i4 and
                              j5 >= i5 and
                              j6 >= i6 and
                              j7 >= i7 and
                              j8 >= i8 and
                              j9 >= i9):
                            cumulative[(i0,i1,i2,i3,i4,i5,i6,i7,i8,i9)][isvar] += count
    return cumulative


  def writeout(self, cumulative):
    """Write cumulative counts.
    Args:
      cumulative nested dict
    """
    print ('clips,asxs,nm,matched,dplx,consensus,indels,fiveprime,threeprime,'
           'ppairs,variant,reference,burden')
    for k in sorted(cumulative.keys()):
      reference = cumulative[k][0]
      variant   = cumulative[k][1]
      (clipped, asxs, nm, matched, dplx, consensus, indels, fiveprime, threeprime, ppairs) = k
      row = []
      row.append(self.clipped[clipped])
      row.append(self.asxs[asxs])
      row.append(self.nm[nm])
      row.append(self.matched[matched])
      row.append(self.dplx[dplx])
      row.append(self.consensus[consensus])
      row.append(self.indels[indels])
      row.append(self.fiveprime[fiveprime])
      row.append(self.threeprime[threeprime])
      row.append(self.ppairs[ppairs])
      row.append(variant)
      row.append(reference)
      row.append(variant/(variant+reference))
      print ','.join(['%s' % r for r in row])


if __name__ == "__main__":
  msg = 'Merge results from gridsearch result files.'
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument('-d', '--directory', required=True)
  args = parser.parse_args()
  MergeResults(args.directory)
