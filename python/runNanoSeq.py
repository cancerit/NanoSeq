#!/usr/bin/python3

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

import sys
from multiprocessing import Pool
import subprocess
import json
import pickle
import os
import argparse
import math
import glob
import gzip
import shutil
import time
import re
import tempfile

version='2.0.0'

parser = argparse.ArgumentParser()
#arguments for all subcommands
parser.add_argument('--out',   action='store', default='.', help='path of the output files and scratch directory (.)')
parser.add_argument('-j','--index', type=int, action='store', help='index of the LSF job array. One based')
parser.add_argument('-k','--max_index', type=int, action='store', help='maximum index of the LSF job array')
parser.add_argument('-t','--threads', type=int, action='store', default= 1, help='number of threads (1)')
parser.add_argument('-v','--version', action='version', version=version)

#subcommands and their arguments
#coverage
subparsers = parser.add_subparsers(dest='subcommand', help='subcommands')
subparsers.required = True #work around for older python versions
parser_cov = subparsers.add_parser('cov', help='coverage calculation')
parser_cov.add_argument('-R','--ref', action='store', required=True, help="referene sequence")
parser_cov.add_argument('-B','--tumour', action='store', required=True, help="tumour (duplex) BAM")
parser_cov.add_argument('--exclude', action='store', default='MT,GL%%,NC_%,hs37d5', help='List of contigs to exclude. Comma separated, %% acts as a wild card. (MT,GL%%,NC_%%,hs37d5)')
parser_cov.add_argument('-w','--win', type=int, action='store', default=100, help='bin size for coverage distribution (100)')
parser_cov.add_argument('-Q', type=int, action='store', default=0, help="minimum mapQ to include a tumour read (0)")
parser_cov.add_argument('--excludeBED', action='store', help='BED file with regions to exclude from analysis')

#dsa
parser_dsa = subparsers.add_parser('dsa', help='compute tables')
parser_dsa.add_argument('-A','--normal', action='store', required=True, help="normal BAM")
parser_dsa.add_argument('-B','--tumour', action='store', required=True, help="tumour (duplex) BAM")
parser_dsa.add_argument('-C','--snp', action='store', required=True, help="SNP BED file")
parser_dsa.add_argument('-D','--mask', action='store', required=True, help="mask BED file")
parser_dsa.add_argument('-R','--ref', action='store', required=True, help="referene sequence")
parser_dsa.add_argument('-d', type=int, action='store', default=2, help="minimum duplex depth (2)")
parser_dsa.add_argument('-q', type=int, action='store', default=30, help="minimum base quality for normal (30)")

#variant caller
parser_var = subparsers.add_parser('var', help='variant caller')
parser_var.add_argument('-U', '--coverage', action='store', required=False, help="file with output coverage")
parser_var.add_argument('-a', type=int, action='store', default=2, help="minimum AS-XS (2)")
parser_var.add_argument('-b', type=int, action='store', default=5, help="minimum duplex reads per strand (5)")
parser_var.add_argument('-c', type=int, action='store', default=0, help="maximum number of clips (0)")
parser_var.add_argument('-f', type=float, action='store', default=0.9, help="minimum fraction of reads for consensus (0.9)")
parser_var.add_argument('-i', type=int, action='store', default=1, help="maximum fraction of reads with an indel (1)")
parser_var.add_argument('-m', type=int, action='store', default=8, help="minimum cycle number (8)")
parser_var.add_argument('-n', type=int, action='store', default=3, help="maximun number of mismatches (3)")
parser_var.add_argument('-p', type=int, action='store', default=0, help="minimum fraction of reads that are proper-pairs (0)")
parser_var.add_argument('-q', type=int, action='store', default=60, help="minimum consensus base quality (60)")
parser_var.add_argument('-r', type=int, action='store', default=144, help="read length (after 5' trimming) (144)")
parser_var.add_argument('-v', type=float, action='store', default=0.01, help="maximum bulk VAF (0.01)")
parser_var.add_argument('-x', type=int, action='store', default=8, help="maximum cycle number (8)")
parser_var.add_argument('-z', type=int, action='store', default=12, help="minimum total number of normal reads (12)")

#compute indels
parser_indel = subparsers.add_parser('indel', help='indel caller')
parser_indel.add_argument('-R','--ref', action='store', required=True, help="referene sequence")
parser_indel.add_argument('-A','--normal', action='store', required=True, help="normal BAM")
parser_indel.add_argument('-B','--tumour', action='store', required=True, help="tumour (duplex) BAM")
parser_indel.add_argument('-s','--sample', action='store', default='sample_1', help="sample name in output vcf")
parser_indel.add_argument('--rb', type=int, action='store', default=2, help="minimum reads in a bundle. (2)")
parser_indel.add_argument('--t3', type=int, action='store', default=135, help="excess bases above this value are trimmed from 3' (135)")
parser_indel.add_argument('--t5', type=int, action='store', default=10, help="bases to trim from 5' reads (10)")
parser_indel.add_argument('--mc', type=int, action='store', default=16, help="minimum bulk coverage (16)")

#carry out gather operations and compute summaries
parser_post = subparsers.add_parser('post', help='gather final results produce summaries')

args = parser.parse_args()

#check job partition arguments
if ( (args.index == None and args.max_index != None) or 
     (args.index != None and args.max_index == None) ) :
    raise ValueError( "Must specify index and max_index for array execution!")

if ( args.index == None and args.max_index == None ) :
  jobArray = False
else :
  jobArray = True

if ( jobArray and args.threads > 1 ):
  print("Warning: execution with a job array is single threaded so thread argument is ignored")

#check the location of the tmpNanoSeq dir
if not os.path.isdir(args.out) :
  parser.error("Specified out directory %s is not accessible!"%args.out)

try :
  testfile = tempfile.TemporaryFile(dir = args.out)
  testfile.close()
except OSError :
  sys.exit("\nCan't write to out directory %s\n"% args.out )


#check required files
if ( hasattr(args,'tumour') ) :
  if not os.path.isfile(args.tumour) :
    parser.error("BAM file %s was not found!" % args.tumour)
  if not os.path.isfile(args.tumour + '.bai') :
    parser.error("Index file %s was not found!" % (args.tumour + '.bai'))

if ( hasattr(args,'normal') ) :
  if not os.path.isfile(args.normal) :
    parser.error("BAM file %s was not found!" % args.normal)
  if not os.path.isfile(args.normal + '.bai') :
    parser.error("Index file %s was not found!" % (args.normal + '.bai'))

if ( hasattr(args,'snp') ) :
  if not os.path.isfile(args.snp) :
    parser.error("SNP file %s was not found!" % args.snp)
  if not os.path.isfile(args.snp + '.tbi') :
    parser.error("Index file %s was not found!" % (args.snp + '.tbi'))

if ( hasattr(args,'mask') ) :
  if not os.path.isfile(args.mask) :
    parser.error("Mask file %s was not found!" % args.mask)
  if not os.path.isfile(args.mask + '.tbi') :
    parser.error("Index file %s was not found!" % (args.mask + '.tbi'))

if ( hasattr(args,'ref') ) :
  if not os.path.isfile(args.ref) :
    parser.error("Reference file %s was not found!" % args.ref)
  if not os.path.isfile(args.ref + '.fai') :
    parser.error("Index file %s was not found!" % (args.ref + '.fai'))

#check that all the code for the analysis is in PATH
scripts = [ "bamcov",
            "dsa",
            "variantcaller",
            "variantcaller_merge.py",
            "variantcaller.R",
            "indelCaller_step1.pl",
            "indelCaller_step2.pl",
            "indelCaller_step3.R",
            "efficiency_nanoseq.R" ]

for icode in scripts :
  if ( shutil.which( icode ) == None ) :
    raise ValueError( "%s was not found in path!"%icode )

#handle genomic intervals with a class
class GInterval :
  def __init__(self, chrr, beg, end ) : 
    self.chr = chrr
    self.beg = beg
    self.end = end
    self.l = end - beg + 1

  def convert2DSAInput(self) :
    #zero based and inclusive of end
    return( GInterval(self.chr,self.beg-1,self.end-1) )

  def __repr__(self) :
    return repr((self.chr,self.beg,self.end))

  def __str__(self):
    return "%s:%s-%s"%(self.chr,self.beg,self.end)

  def __contains__(self,other):
    if (self.chr != other.chr ) : 
      return False
    if ( other.end < self.beg or self.end < other.beg ) :
      return False
    return True
  
  def __sub__(self,other) :
    if (self.chr != other.chr ) : 
      return [ self ]
    if ( other.beg > self.beg and other.beg <= self.end ) :
      if ( other.end < self.end ) : 
        return [ GInterval(self.chr,self.beg,other.beg-1),GInterval(self.chr,other.end+1,self.end) ] 
      if ( other.end >= self.end) :
        return [ GInterval(self.chr,self.beg, other.beg -1 ) ]
    if ( other.beg <= self.beg ) : 
      if ( other.end < self.end and other.end >= self.beg ) :
        return [ GInterval(self.chr,other.end + 1, self.end) ]
      if ( other.end >= self.end ) :
        return []
    return [self]

  def __add__(self,other) :
    if ( other.chr != self.chr ) :
      return [self, other]
    if ( other.end +1 == self.beg or self.end +1 == other.beg ) :
      return [ GInterval(self.chr, min([self.beg,other.beg]),max([self.end,other.end])) ] 
    if ( other.end < self.beg or self.end < other.beg ) :
      return [ self, other]
    if ( self.beg <= other.end or other.beg <= self.end) :
      return [ GInterval(self.chr, min([self.beg,other.beg]),max([self.end,other.end])) ]

  ## define to allow sorting of intervals

  #sort order for chromosomes
  #Compares two chromosome names with the order being -
  #     1. Simple names before long names
  #     2. Names with 'chr' then number (or just bare number) before Names with 'chr' then letter(s)
  #     3. 'chrX', 'X', 'chrY', 'Y', 'M', 'MT' after 2.
  #     4. Others
  def chr_isless( x, y) :
    if ( x == y ) : return False

    xc = re.sub('^chr','', x)
    yc = re.sub('^chr','', y)
    
    if ( xc.isdigit() ) :
      if ( yc.isdigit() ) :
        return int(xc) < int(yc)
      else :
        return True # rule2
    elif ( yc.isdigit() ) :
      return False # rule2
    elif ( ( xc == "M" or xc == "MT" ) and ( yc == "X" or yc == "Y") ):
      return False # rule 3
    elif ( ( xc == "X" or xc == "Y" ) and ( yc == "M" or yc == "MT") ):
      return True # rule 3
    elif ( ( xc == "M" or xc == "MT" ) ):
      return True # rule 3
    elif ( ( yc == "M" or yc == "MT") ) :
      return False # rule 3
    else :
      return xc < yc
      
   
  def __eq__(self,other) :
    if ( self.chr == other.chr ) :
      if ( self.beg == other.beg ) :
        return True
    return False

  def __ne__(self,other) :
    if ( self.chr == other.chr ) :
      if ( self.beg == other.beg ) :
        return False
    return True
  
  def __gt__(self,other): 
    if ( GInterval.chr_isless(other.chr,self.chr  ) ) :
      return True
    if ( GInterval.chr_isless(self.chr, other.chr) ) :
      return False
    if ( self.beg > other.beg ) :
      return True
    else :
      return False

  def __ge__(self,other): 
    if ( GInterval.chr_isless(other.chr, self.chr )) :
      return True
    if ( GInterval.chr_isless(self.chr, other.chr )) :
      return False
    if ( self.beg >= other.beg ) :
      return True
    else :
      return False

  def __lt__(self,other): 
    if ( GInterval.chr_isless(self.chr, other.chr )) :
      return True
    if ( GInterval.chr_isless(other.chr,self.chr )) :
      return False
    if ( self.beg < other.beg ) :
      return True
    else :
      return False

  def __le__(self,other): 
    if ( GInterval.chr_isless(self.chr, other.chr )) :
      return True
    if ( GInterval.chr_isless(other.chr, self.chr )) :
      return False
    if ( self.beg <= other.beg ) :
      return True
    else :
      return False

# compute coverage histogram
def runBamcov(bam,mapQ, window, ichr, out) :
  out = out+".cov.bed"
  p = subprocess.Popen(['bamcov','-q',mapQ,'-w',window,'-r',ichr,'-o',out, bam],stderr=subprocess.PIPE)
  p.wait()
  if (p.returncode != 0 ) : 
    error = p.stderr.read().decode()
    print("\n!Error running bamcov for chr %s (cov): %s\n"%(ichr,error))
    raise ValueError(error)
  p = subprocess.Popen(['gzip', '-1', '-f', out],stderr=subprocess.PIPE)
  p.wait()
  if (p.returncode != 0 ) : 
    error = p.stderr.read().decode()
    print("\n!Error while compressing %s with gzip (cov) : %s\n"%(out,error))
    raise ValueError(error)
  outdone = re.sub('cov.bed$','done',out)
  open(outdone,'w').close()
  return

#run an external command
def runCommand(command) :
  for ijob in command.rstrip(';').split(';') :
    print("\nExecuting: %s\n"%ijob)
    p = subprocess.Popen(ijob,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()
    if (p.returncode != 0 ) :
      error = p.stderr.read().decode()
      print("\n!Error processing:  %s\n"%ijob )
      raise ValueError(error)
  return

print("command arguments :")

print(args.__dict__)
print()
# create all directory tree for the tmp files
tmpDir = args.out + "/tmpNanoSeq"
if (args.subcommand == 'cov'):
  if ( args.index == None or args.index == 1 ) :
    if (not os.path.isdir(tmpDir+'/dsa') ) :
      os.makedirs(tmpDir+'/dsa')
    if (not os.path.isdir(tmpDir+'/var') ) :
      os.makedirs(tmpDir+'/var')
    if (not os.path.isdir(tmpDir+'/post') ) :
      os.makedirs(tmpDir+'/post')
    if (not os.path.isdir(tmpDir+'/cov') ) :
      os.makedirs(tmpDir+'/cov')
    if (not os.path.isdir(tmpDir+'/indel') ) :
      os.makedirs(tmpDir+'/indel')

with open("%s/%s/args.json"%(tmpDir,args.subcommand), "w") as jsonOut :
    json.dump(args.__dict__, jsonOut)

#coverage section
if (args.subcommand == 'cov'):
  # build the chromosome dictionary, list and intervals
  excludes = [ re.compile(istr + "$") for istr in args.exclude.replace("%",".+").split(',') ]
  chrList = []
  rnames = {}
  with open(args.ref + '.fai','r') as iofile :
    for iline in iofile :
      ichr = iline.split('\t')[0]
      if any( iregx.match(ichr) for iregx in excludes ): continue
      ilength = iline.split('\t')[1]
      chrList.append(ichr)
      rnames[ichr] = int(ilength)

  gintervals = []
  for ichr in chrList :
    gintervals.append( GInterval(ichr,1,rnames[ichr]) )
  gintervals.sort()
  reorderchr = []
  for iint in gintervals :
    reorderchr.append( iint.chr)
  chrList = reorderchr

  #remove regions to ignore from exclude BED
  if ( args.excludeBED != None ) :
    xIntervals = []
    with open(args.excludeBED,'r') as iofile :
      for iline in iofile :
        ichr = str(iline.split('\t')[0])
        ib = int(iline.split('\t')[1])
        ie = int(iline.split('\t')[2])
        xIntervals.append( GInterval(ichr,ib+1,ie)) #convert to ginterval 

    xIntervals.sort()
    with open("%s/cov/%s"%(tmpDir,'gIntervalsExclude.dat'), 'wb') as iofile :
      pickle.dump(xIntervals,iofile)
        
    iiresult= []
    for ii in gintervals :
      ifrag = ii
      for jj in xIntervals :
        if ( ifrag.chr != jj.chr ) :continue
        diff = ifrag - jj
        if len(diff) == 2 :
          iiresult.append( diff[0] )
          ifrag = diff[1]
        else :
          ifrag = diff[0]
      iiresult.append( ifrag )

    gintervals = iiresult

  print("Starting coverage calculation\n")
  with open("%s/cov/nfiles"%(tmpDir), "w") as iofile :
    iofile.write(str(len(chrList) ))
  inputs = []
  for ii,ichr in enumerate(chrList) :
    inputs.append( (args.tumour,str(args.Q),str(args.win),ichr,"%s/cov/%s"%(tmpDir,ii) ))
  if ( args.index == None ) :
    #multiple threads are available
    with open("%s/cov/%s"%(tmpDir,'gIntervals.dat'), 'wb') as iofile :
      pickle.dump(gintervals,iofile)
    with Pool( args.threads ) as p:
      p.starmap(runBamcov, inputs)
    print("Completed coverage calculation\n")

  else :
    #single thread execution
    if (args.index == 1) :
      with open("%s/cov/%s"%(tmpDir,'gIntervals.dat'), 'wb') as iofile :
        pickle.dump(gintervals,iofile)

    jobsPerCPU = math.ceil( len(inputs)/args.max_index )
    for ii in range(jobsPerCPU * (args.index-1), jobsPerCPU*args.index) :
      if ( ii >= len( inputs) ) : break
      runBamcov(inputs[ii][0], inputs[ii][1], inputs[ii][2], inputs[ii][3], inputs[ii][4])
    print("Completed coverage calculation for this job\n")

#dsa section
if (args.subcommand == 'dsa' ) :
  if (not os.path.isfile(tmpDir+'/cov/args.json') ):
    sys.exit("\nMust run cov submmand prior to dsa\n")
  else :
    with open(tmpDir+'/cov/args.json') as iofile :
      oargs = json.load( iofile )

  if (not os.path.isfile(tmpDir+'/cov/nfiles') ):
    sys.exit(tmpDir+'/cov/nfiles not found!\n')
  else :
    with open(tmpDir+'/cov/nfiles') as iofile :
      nfiles = int(iofile.readline())

  for i in range(nfiles ) :
    if ( len(glob.glob(tmpDir+"/cov/%s.done"%i)) != 1 ) :
      sys.exit("\ncov job %s did not complete correctly\n"%i)

  #try to stagger file access for array execution
  if ( args.index != None ) : time.sleep( 3 * args.index )

  with open(tmpDir+"/cov/gIntervals.dat", 'rb') as iofile :
    gIntervals = pickle.load(iofile)

  coverage = []
  cctotal = 0
  chrOffset = {}
  print("\nParsing coverage files\n")
  for i in range(nfiles) :
    with gzip.open( tmpDir+"/cov/%s.cov.bed.gz"%i,'rt') as iofile :
      for iline in iofile :
        ichr = str(iline.split('\t')[0])
        ib = int(iline.split('\t')[1])
        cc = int(iline.split('\t')[3])
        cctotal += cc
        if (ib == 0 ) : chrOffset[str(ichr)] = len( coverage )
        coverage.append( [ichr,ib,cc]) #convert to ginterval 
  print("\nCompleted parsing coverage files\n")

  #correct cctotal if there are exclusion regions
  if ( oargs['excludeBED'] != None ) :
    with open(tmpDir+"/cov/gIntervalsExclude.dat", 'rb') as iofile :
      xIntervals = pickle.load(iofile)
    xSumCov = 0
    for iinterval in xIntervals :
      ichar =iinterval.chr
      ibeg = iinterval.beg - 1
      iend = iinterval.end
      for i in range(math.floor(ibeg/oargs['win']), math.floor(iend/oargs['win'])+ 1 ):
        j = i + chrOffset[ ichar ]
        xSumCov += coverage[j][2]
    cctotal = cctotal - xSumCov

  #determine the genomic intervals to give to each job so that each one roughly
  #has the same amount of coverage
  basesPerCPU = 0
  if ( args.index == None ) :
    njobs = args.threads
  else :
    njobs = args.max_index
  basesPerCPU = cctotal/ njobs
  print( "\nPartitioning %s jobs with %s bases/task\n"%(njobs,basesPerCPU))
  with open("%s/dsa/nfiles"%(tmpDir), "w") as iofile :
    iofile.write(str(njobs ))

  sumCov = 0
  oIntervals = []
  intervalsPerCPU = []
  while ( len(gIntervals) > 0 ):
    iinterval = gIntervals.pop(0)
    ichar =iinterval.chr
    ibeg = iinterval.beg - 1
    iend = iinterval.end
    for i in range(math.floor(ibeg/oargs['win']), math.floor(iend/oargs['win'])+ 1 ):
      j = i + chrOffset[ ichar ]
      sumCov += coverage[j][2]
      if ( sumCov > basesPerCPU ) :
        oIntervals.append(GInterval(ichar, ibeg + 1, coverage[j][1] + oargs['win']))
        intervalsPerCPU.append( oIntervals)
        oIntervals = [] 
        sumCov = 0
        ibeg=  coverage[j][1] + oargs['win'] 
    oIntervals.append(GInterval(ichar, ibeg + 1, iend))
  if ( len(oIntervals) > 0 ) : intervalsPerCPU.append( oIntervals)
  
  commands = []
  for i in range(njobs) :
    #for restarts remove job files that didn't complete correctly
    if ( not os.path.isfile("%s/dsa/%s.done"%(tmpDir,i)) ) :
      if ( os.path.isfile("%s/dsa/%s.bed"%(tmpDir,i) ) ) :  os.remove( "%s/dsa/%s.bed"%(tmpDir,i))
      if ( os.path.isfile("%s/dsa/%s.bed.gz"%(tmpDir,i) ) ) :  os.remove( "%s/dsa/%s.bed.gz"%(tmpDir,i))
    else :
      continue

    #construct dsa commands
    cmd = ""
    for iinterval in intervalsPerCPU[i] :
      dsaInt = iinterval.convert2DSAInput()
      cmd += "dsa -A %s -B %s -C %s -D %s -R %s -d %s -Q %s -M %s -r %s -b %s -e %s >> %s/dsa/%s.dsa.bed;" \
              %(args.normal, args.tumour, args.snp, args.mask, args.ref, args.d, args.q, oargs['Q'],
                dsaInt.chr, dsaInt.beg, dsaInt.end,tmpDir,i)
    cmd += "gzip -1 %s/dsa/%s.dsa.bed; sleep 2; gzip -t %s/dsa/%s.dsa.bed.gz;"%(tmpDir,i,tmpDir,i)
    cmd += "touch %s/dsa/%s.done"%(tmpDir,i)
    commands.append( ( cmd , ))

  #execute dsa commans
  print("Starting dsa calculation\n")
  if ( args.index == None ) :
    #multithread
    with Pool( args.threads ) as p :
      p.starmap( runCommand, commands )
  else :
    #array execution
    runCommand( commands[args.index][0] )