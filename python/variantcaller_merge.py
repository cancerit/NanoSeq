#!/usr/bin/python

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

import argparse
import os
import subprocess
import glob
import sys
class MergeResults(object):

  """Merge results from multiple BED files and output metrics files
  """

  def __init__(self, dirname):
    """Initialize.
    Args:
      dirname: directory containing variantcaller results files
    """
    self.dirname  = dirname
    self.bedfiles = self.bedfiles()
    self.outfiles = self.outfiles()
    self.headers()
    self.read_and_split()
    self.close()

  def bedfiles(self):
    """Returns list of BED filenames for all .bed files in the directory
    """
    #cmp function to sort iles on natural order of chromosomes and coordinates
    def cmp_files(a,b):
      filea = a.split('/')[-1]
      fileb = b.split('/')[-1]
      try :
        chra = int(filea.split('.')[0])
      except :
        chra = filea.split('.')[0]
      try :
        chrb = int(fileb.split('.')[0])
      except :
        chrb = fileb.split('.')[0]
      ia = int(filea.split('.')[1])
      ib = int(fileb.split('.')[1])
      if ( chra > chrb ) : return 1
      elif ( chra < chrb) : return -1
      elif ( ia > ib ) : return 1
      elif ( ia < ib ) : return -1
      else : return 0

    fns = []
    listfiles = glob.glob("%s/*.bed"%self.dirname)
    listfiles += glob.glob("%s/*.var"%self.dirname)
    listfiles.sort(cmp_files) 
    for fn in listfiles :
      if (not fn.endswith('.csv')):
        fns.append(fn)
    #assert len(fns) == 298
    return fns


  def outfiles(self):
    """Returns dict with key = filename, value = output file
    """
    return {
      'Params'     : open('%s/%s' % (self.dirname, 'params.csv'), 'w'),
      'Coverage'   : open('%s/%s' % (self.dirname, 'coverage.csv'), 'w'),
      'CallVsQpos' : open('%s/%s' % (self.dirname, 'callvsqpos.csv'), 'w'),
      'PyrVsMask'  : open('%s/%s' % (self.dirname, 'pyrvsmask.csv'), 'w'),
      'ReadBundles': open('%s/%s' % (self.dirname, 'readbundles.csv'), 'w'),
      'Burdens'    : open('%s/%s' % (self.dirname, 'burdens.csv'), 'w'),
      'Variants'   : open('%s/%s' % (self.dirname, 'variants.csv'), 'w'),
      'Mismatches' : open('%s/%s' % (self.dirname, 'mismatches.csv'), 'w')
      }


  def headers(self):
    """Write headers to output files
    """
    self.outfiles['Coverage'].write('count\n')
    self.outfiles['CallVsQpos'].write('base,qpos,ismasked,count\n')
    self.outfiles['PyrVsMask'].write('pyrcontext,ismasked,count\n')
    self.outfiles['ReadBundles'].write('fwd,rev,ismasked,isvariant,count\n')
    self.outfiles['Burdens'].write('ismasked,isvariant,count\n')
    self.outfiles['Variants'].write('chrom,chromStart,context,commonSNP,'
      'shearwater,bulkASXS,bulkNM,bulkForwardA,bulkForwardC,bulkForwardG,'
      'bulkForwardT,bulkForwardIndel,bulkReverseA,bulkReverseC,bulkReverseG,'
      'bulkReverseT,bulkReverseIndel,dplxBreakpointBeg,dplxBreakpointEnd,'
      'bundleType,dplxASXS,dplxCLIP,dplxNM,dplxfwdA,dplxfwdC,dplxfwdG,dplxfwdT,'
      'dplxfwdIndel,dplxrevA,dplxrevC,dplxrevG,dplxrevT,dplxrevIndel,'
      'dplxCQfwdA,dplxCQfwdC,dplxCQfwdG,dplxCQfwdT,dplxCQrevA,'
      'dplxCQrevC,dplxCQrevG,dplxCQrevT,bulkForwardTotal,bulkReverseTotal,'
      'dplxfwdTotal,dplxrevTotal,left,right,qpos,call,isvariant,pyrcontext,'
      'stdcontext,pyrsub,stdsub,ismasked\n')
    self.outfiles['Mismatches'].write('chrom,chromStart,context,commonSNP,'
      'shearwater,bulkASXS,bulkNM,bulkForwardA,bulkForwardC,bulkForwardG,'
      'bulkForwardT,bulkForwardIndel,bulkReverseA,bulkReverseC,bulkReverseG,'
      'bulkReverseT,bulkReverseIndel,dplxBreakpointBeg,dplxBreakpointEnd,'
      'dplxASXS,dplxCLIP,dplxNM,dplxfwdA,dplxfwdC,dplxfwdG,dplxfwdT,'
      'dplxfwdIndel,dplxrevA,dplxrevC,dplxrevG,dplxrevT,dplxrevIndel,'
      'dplxCQfwdA,dplxCQfwdC,dplxCQfwdG,dplxCQfwdT,dplxCQrevA,'
      'dplxCQrevC,dplxCQrevG,dplxCQrevT,bulkForwardTotal,bulkReverseTotal,'
      'dplxfwdTotal,dplxrevTotal,left,right,qpos,mismatch,ismasked\n')


  def read_and_split(self):
    """Read files from individual beds and output into combined metric files
    """
    for bed in self.bedfiles:
      for row in open(bed, 'rU'):
        if ( row[0] == '#') :
          self.outfiles['Params'].write(row[2:])
          continue
        row = row.strip().split('\t')
        if self.outfiles.get(row[0], None):
          self.outfiles[row[0]].write('%s\n' % ','.join(row[1:]))


  def close(self):
    """Close all outfiles
    """
    for k,v in self.outfiles.iteritems():
      v.close()


if __name__ == "__main__":
  msg = 'Merge results from multiple duplex sequencing result files.'
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument('-d', '--directory', required=True)
  args = parser.parse_args()
  MergeResults(args.directory)
