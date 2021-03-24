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
import pickle
import os
import getopt
import math
import glob
import shutil
import time
import re

version='1.5.2'

usage = """
  runBotSeq.py  [general options] dsa-options variantcall-options

     general paramameters :

      -o                location of outputs and temporary files ( pwd )
      -j                one-based job index (array jobs), must set -t to max index value
      -t                number of threads (default 1)
      -s                only do the coverage calculation step
      -R                Reference sequence file (required )
      -w                Bin size of precalculation coverage ( 20000 )
      -X                BED file with regions to be excluded
      -H                List of contigs from the reference to exclude. 
                          Comma separated using '%' as wildcard.   ( MT,GL%,NC_%,hs37d5 )
      -T                Archive table and variants directories
   --no-post            Skip post-analysis steps
   --version            Version of the code

     dsa  required parameters:

      -B, --duplexBAM   duplex BAM file aka BotSeq
      -A, --bulkBAM     bulk BAM file aka Normal
      -C                SNP BED file name
      -D                Mask BED file name

    dsa optional parameters:

      -d                Minimum duplex depth ( 2 )
      -Q                Minimum base quality for bulk sequencing ( 30 )
      -M                Minimum MAPQ of the duplex reads ( 0 )

    variantcaller optional parameters:

      -a                minimum AS-XS (50)
      -b                minimum number of bulk reads per strand (5)
      -z                minimum number of bulk reads in total (12)
      -c                maximum number of clips (0)
      -f                minimum fraction of reads for consensus (0.9)
      -i                maximum fraction of reads with an indel (1)
      -n                maximun number of mismatches (3)
      -p                minimum fraction of reads that are proper-pairs (0)
      -q                minimum consensus base quality (60)
      -r                read length (after 5' trimming) (144)
      -m                minimum cycle number (8)
      -x                maximum cycle number (8)
      -v                maximum bulk VAF (0.01)

"""
argv = sys.argv[1:]
if ( len(argv) == 0 ):
  print(usage)
  sys.exit(0)

# count the number of reads in an interval with samtools
def countReads(bam,mapQ, intervals) :
  strIntervals = [ str(i) for i in intervals ]
  p = subprocess.Popen(['samtools','view','-q',mapQ,'-Mc',bam] + strIntervals,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  p.wait()
  if (p.returncode != 0 ) : 
    error = p.stderr.read().decode()
    print("\n!Error with samtools processing : %s\n"%error)
    raise ValueError(error)
  nReads = int(p.stdout.read().decode().rstrip('\n'))
  return nReads

#run an external command
def runCommand(command) :
  for ijob in command.rstrip('\n').split('\n') :
    p = subprocess.Popen(ijob,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()
    if (p.returncode != 0 ) :
      error = p.stderr.read().decode()
      print("\n!Error processing:  %s\n"% command)
      raise ValueError(error)
  return

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
    if ( self.chr > other.chr ) :
      return True
    if ( self.chr < other.chr ) :
      return False
    if ( self.beg > other.beg ) :
      return True
    else :
      return False

  def __ge__(self,other): 
    if ( self.chr > other.chr ) :
      return True
    if ( self.chr < other.chr ) :
      return False
    if ( self.beg >= other.beg ) :
      return True
    else :
      return False

  def __lt__(self,other): 
    if ( self.chr < other.chr ) :
      return True
    if ( self.chr > other.chr ) :
      return False
    if ( self.beg < other.beg ) :
      return True
    else :
      return False

  def __le__(self,other): 
    if ( self.chr < other.chr ) :
      return True
    if ( self.chr > other.chr ) :
      return False
    if ( self.beg <= other.beg ) :
      return True
    else :
      return False


opts, args = getopt.getopt(argv, 'hsj:t:A:B:C:D:R:d:Q:M:H:TX:a:c:n:m:o:x:b:f:i:p:q:r:v:w:z:', ['duplexBAM=','bulkBAM=','no-post','version'])
#default arguments for dsa
dsaOpts = { '-d': 2, '-Q': 30, '-M': 0 }
#variantcaller params from one of Fede's scripts
varOpts = { '-a': 50, '-c': 0, '-n': 3, '-m': 8, '-x': 8, '-b': 5, '-f': 0.9, '-i': 1, '-p': 0, '-q': 60, '-r': 144, '-v': 0.01,'-z':12}
workDir = os.getcwd()
jobIndex = None
onlySplitStep = False
duplexBAM = None
bulkBAM = None
refFile = None
minMapQ = 0
nCPU = 1
window =20000
excludeBED = None
archive=False
noPost=False
excludeChr = "MT,GL%,NC_%,hs37d5"

for iopt, iarg in opts:
  if iopt in ('-A','--bulkBAM'):
    bulkBAM = iarg
    dsaOpts['-A'] = iarg
  elif iopt in ('-B','--duplexBAM'):
    duplexBAM = iarg
    dsaOpts['-B'] = iarg
  elif iopt == '-j':
    jobIndex = int(iarg)
  elif iopt == '-s':
    onlySplitStep = True
  elif iopt == '-t':
    nCPU = int(iarg)
  elif iopt == '-o':
    workDir = os.path.abspath(iarg)
  elif iopt == '-C':
    dsaOpts['-C'] = iarg
  elif iopt == '-D':
    dsaOpts['-D'] = iarg
  elif iopt == '-R':
    refFile = iarg
  elif iopt == '-d':
    dsaOpts['-d'] = iarg
  elif iopt == '-Q':
    dsaOpts['-Q'] = iarg
  elif iopt == '-M':
    dsaOpts['-M'] = iarg
    minMapQ = int(iarg)
  elif iopt == '-b':
    varOpts['-b'] = iarg
  elif iopt == '-z':
    varOpts['-z'] = iarg
  elif iopt == '-c':
    varOpts['-c'] = iarg
  elif iopt == '-f':
    varOpts['-f'] = iarg
  elif iopt == '-i':
    varOpts['-i'] = iarg
  elif iopt == '-n':
    varOpts['-n'] = iarg
  elif iopt == '-p':
    varOpts['-p'] = iarg
  elif iopt == '-q':
    varOpts['-q'] = iarg
  elif iopt == '-r':
    varOpts['-r'] = iarg
  elif iopt == '-m':
    varOpts['-m'] = iarg
  elif iopt == '-x':
    varOpts['-x'] = iarg
  elif iopt == '-H':
    excludeChr = iarg
  elif iopt == '-X':
    excludeBED = iarg
  elif iopt == '-v':
    varOpts['-v'] = iarg
  elif iopt == '-w':
    window = int(iarg)
  elif iopt == '-T':
    archive = True
  elif iopt == '--no-post':
    noPost = True 
  elif iopt == '--version':
    print(version)
    sys.exit(0)
  elif iopt == '-h':
    print(usage)
    sys.exit(0)  

#check if samtools is in path
if ( shutil.which('samtools') == None ) : 
  raise ValueError('samtools not found found in path!')

#check if dsa is in path
if ( shutil.which('dsa') == None ) :
  raise ValueError('dsa not found found in path!')

#check if samtools is in path
if ( shutil.which('variantcaller') == None ) :
  raise ValueError('variantcaller not found found in path!')

#check if merge script is in path
if ( shutil.which('variantcaller_merge.py') == None ) :
  raise ValueError('variantcaller_merge.py not found found in path!')

#check if summary script is in path
if ( shutil.which('variantcaller.R') == None ) :
  raise ValueError('variantcaller.R not found found in path!')

#check validity of the job index
if ( jobIndex != None ) :
  if ( jobIndex < 1 or jobIndex > nCPU ) :
    raise ValueError( "Job index must be a value between 1 - nCPU ")

#These parameters appear to be the same. Need to check with Fede
varOpts['-d'] = dsaOpts['-d']
dsaOpts['-R'] = refFile

for iopt in ['-C','-D']:
  if ( not iopt in dsaOpts ) :
    raise ValueError("Require definition of dsa parameter %s"%(iopt))
  if ( not os.path.exists( dsaOpts[iopt] ) ) :
    raise ValueError("Option %s , file %s not found"%(iopt,dsaOpts[iopt]))

if( refFile == None) :
  raise ValueError("Need to a reference file fasta ( -R )")
if ( not os.path.exists( refFile ) ) :
  raise ValueError("Reference file %s was not found!"%(refFile))

if ( not os.path.exists( refFile + '.fai' ) ) :
  #try to build index file if not found
  runCommand("samtools faidx %s"% refFile)

if( duplexBAM == None) :
  raise ValueError("Need to specify a duplex BAM file (--duplexBAM )")
if ( not os.path.exists( duplexBAM ) ) :
  raise ValueError("Duplex BAM %s was not found!"%(duplexBAM))

if( bulkBAM == None) :
  raise ValueError("Need to specify a bulk BAM file (--bulkBAM )")
if ( not os.path.exists( bulkBAM ) ) :
  raise ValueError("Duplex BAM %s was not found!"%(bulkBAM))

if ( excludeBED != None ) :
  if ( not os.path.exists( excludeBED ) ) :
    raise ValueError("Exclude BED file %s not found!"%(exludeBED))

if( nCPU == 0 ) :
  raise ValueError("Need number of execution threads")


# build the chromosome dictionary, list and intervals
excludes = [ re.compile(istr + "$") for istr in excludeChr.replace("%",".+").split(',') ]
chrList = []
rnames = {}
with open(refFile + '.fai','r') as iofile :
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

#remove regions to ignore from exclude BED
if ( excludeBED != None ) :
  xIntervals = []
  with open(excludeBED,'r') as iofile :
    for iline in iofile :
      ichr = str(iline.split('\t')[0])
      ib = int(iline.split('\t')[1])
      ie = int(iline.split('\t')[2])
      xIntervals.append( GInterval(ichr,ib+1,ie)) #convert to ginterval
  
  xIntervals.sort()
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

#check that tables and variants directories are present
if ( jobIndex == None or jobIndex == 1 ) :
  if (not os.path.isdir(workDir+'/tables') ) :
    os.makedirs(workDir+'/tables/tmp')
  if (not os.path.isdir(workDir+'/variants') ) :
    os.makedirs(workDir+'/variants/tmp')
  if (not os.path.isdir(workDir+'/summary') ) :
    os.makedirs(workDir+'/summary')
else :
  time.sleep(6) #give timeout for directory setup for multijob exec

#generate intervals for coverage calculation of duplex BAM
newIntervals = []
last = 0
iacc = 0
ii  = 0
for iinterval in gintervals :
  ipos = iinterval.beg+  last
  for ipos in range(iinterval.beg + window + last, iinterval.end , window):
    iacc =0
    last = 0
    newIntervals.append(GInterval(iinterval.chr, ipos - window, ipos - 1) )
  if ipos <= iinterval.end:
    newIntervals.append(GInterval(iinterval.chr, ipos, iinterval.end  ) )
    iacc += iinterval.end - ipos + 1
    last = 0
  if  ii < len(gintervals ) - 1 :
    if  (iacc + gintervals[ii+1].l > window) :
      newIntervals.append(GInterval(gintervals[ii + 1].chr, gintervals[ii+1].beg, gintervals[ii+1].beg + window - iacc -1))
      last =  window  - iacc 
      iacc = 0
  ii += 1

#group intervals so each job covers a window
#samtools intervals are closed 1-based
ilist = []
samIntervals = []
isum  = 0
nloci = len(newIntervals)
for ii in range(nloci) :
  iinterval = newIntervals[ii]
  isum += iinterval.l
  if ( isum < window  ) :
    ilist.append( iinterval)
  else :
    isum = 0 
    ilist.append( iinterval)
    samIntervals.append(ilist)
    ilist = [] 
  if ( ii == nloci -1 ) :
    samIntervals.append(ilist)
#newIntervals =[] #dont need anymore

if ( not os.path.exists( "%s/tables/tmp/%s"%(workDir,'allCounts.done') ) ) :

  print("\nStarted coverage calculation\n")

  #prepare samtools inputs
  inputs = []
  for iinterval in samIntervals :
    inputs.append( (duplexBAM,str(minMapQ),iinterval ))

  if ( jobIndex == None ) :
    #multiple threads are available
    with Pool( nCPU ) as p:
      counts =  p.starmap(countReads, inputs)
    with open("%s/tables/tmp/%s"%(workDir,'allCounts.dat'), 'wb') as iofile :
      pickle.dump(counts,iofile)
    outIntervals = []
    for ii in samIntervals :
      ilist = [ jj.convert2DSAInput() for jj in ii ]
      outIntervals.append(ilist)
    with open("%s/tables/tmp/%s"%(workDir,'gIntervals.dat'), 'wb') as iofile :
      pickle.dump(outIntervals,iofile)
    command = "touch %s/tables/tmp/allCounts.done"%(workDir)
    runCommand( command )

    print("Completed coverage calculation\n")

  else :
    #single thread execution
    jobsPerCPU = math.ceil( len(inputs)/nCPU )
    iicounts = []
    for ii in range(jobsPerCPU * (jobIndex-1), jobsPerCPU*jobIndex) :
      if ( ii >= len( inputs) ) : break
      iicounts.append( countReads(inputs[ii][0],inputs[ii][1],inputs[ii][2]  ))
    with open("%s/tables/tmp/%s"%(workDir,"count.%s.dat"%jobIndex), 'wb') as iofile :
      pickle.dump(iicounts,iofile)
    command = "touch %s/tables/tmp/count.%s.done"%(workDir,jobIndex)
    runCommand( command )

    #wait until all the jobs complete
    while ( not( len(glob.glob("%s/tables/tmp/count.*.done"%(workDir)) ) == nCPU or 
             len(glob.glob("%s/tables/tmp/allCounts.done"%(workDir)) ) == 1 ) ) :
      time.sleep(5)

    #combine the results of all the jobs 
    counts =[]
    if ( jobIndex == 1 ) :
      for ii in range(nCPU) :
        ifname = "%s/tables/tmp/count.%s.dat"%(workDir,ii + 1)
        jfname = "%s/tables/tmp/count.%s.done"%(workDir,ii + 1)
        with open(ifname, 'rb') as iofile :
          counts.extend( pickle.load(iofile) )
        os.remove(ifname)
        os.remove(jfname)
      with open("%s/tables/tmp/%s"%(workDir,'allCounts.dat'), 'wb') as iofile :
        pickle.dump(counts,iofile)
      outIntervals = []
      for ii in samIntervals :
        ilist = [ jj.convert2DSAInput() for jj in ii ]
        outIntervals.append(ilist)
      with open("%s/tables/tmp/%s"%(workDir,'gIntervals.dat'), 'wb') as iofile :
        pickle.dump(outIntervals,iofile)
      command = "touch %s/tables/tmp/allCounts.done"%(workDir)
      runCommand( command )

      print("Completed coverage calculation\n")

    else :
      while ( len(glob.glob("%s/tables/tmp/allCounts.done"%(workDir)) ) != 1 )  :
        time.sleep(5)

if ( onlySplitStep ) : sys.exit(0)

#load the counts and intervals so that all jobs have them
with open("%s/tables/tmp/%s"%(workDir,'allCounts.dat'), 'rb') as iofile :
  counts = pickle.load(iofile)
with open("%s/tables/tmp/%s"%(workDir,'gIntervals.dat'), 'rb') as iofile :
  gIntervals = pickle.load(iofile)

#check widown value didn't change when rerunning
isum =0
for ii in gIntervals[0] :
  isum += ii.l
if ( isum != window ) :
  raise ValueError("Specified window %s doesn't match precomputed intervals")

#compute the intervals used for dsa and variantcaller
totalReads = sum(counts)
readsPerCPU = totalReads/nCPU
dsaIntervals = [ [] for i in range(nCPU) ]
ii = 0
accumulatedReads = 0
icpu = 0
ranges = []
nreads= [0] * nCPU 
for ii in range(len(counts)) :
  accumulatedReads += counts[ii]
  if (icpu < nCPU ) : nreads[icpu] += counts[ii]
  ranges.extend( [ jj for jj in gIntervals[ii] ])
  if ( accumulatedReads  >= (icpu+1)*readsPerCPU ) :
    merged = [ ranges[0] ]
    for iranges in ranges :
      merged.extend(merged.pop() + iranges )
    dsaIntervals[icpu].extend( merged )
    ranges = []
    icpu += 1
if ( len(ranges) > 0) :
  merged = [ dsaIntervals[-1].pop() ]
  for iranges in ranges :
    merged.extend(merged.pop() + iranges )
  dsaIntervals[icpu-1].extend( merged )

#construct the command strings to be run in each thread
jobs = []
nMarkerFiles = 0
nTask2Run = 0
for ii  in dsaIntervals :
  commands = "set -e;"
  for idsa in ii :
    #skip commands that have already completed succesfully (markerFile is present)
    markerFile1 = "%s/tables/tmp/%s.dsa.done"%(workDir,nMarkerFiles)
    markerFile2 = "%s/variants/tmp/%s.var.done"%(workDir,nMarkerFiles)
    nMarkerFiles += 1
    if ( not os.path.exists( markerFile1 )) :
      nTask2Run += 1
      #build table step (most time consuming)
      commands += "dsa -A %s -B %s -C %s -D %s -Q %s -M %s -R %s -d %s -r %s -b %s -e %s -O %s/tables/%s.%s.%s.tbl;"%(dsaOpts['-A'], dsaOpts['-B'], dsaOpts['-C'], dsaOpts['-D'], dsaOpts['-Q'],dsaOpts['-M'], dsaOpts['-R'], dsaOpts['-d'], idsa.chr,idsa.beg,idsa.end,workDir,idsa.chr,idsa.beg,idsa.end)
      commands += "touch %s;"%(markerFile1)

    if ( not os.path.exists( markerFile2 )) :
      nTask2Run += 1
      #variantcaller step
      commands +="variantcaller -B %s/tables/%s.%s.%s.tbl.gz -O %s/variants/%s.%s.%s.var -a %s -b %s -c %s -d %s -f %s -i %s -m %s -n %s -p %s -q %s -r %s -v %s -x %s -z %s;"%(workDir,idsa.chr,idsa.beg,idsa.end,workDir,idsa.chr,idsa.beg,idsa.end,varOpts['-a'],varOpts['-b'],varOpts['-c'],varOpts['-d'],varOpts['-f'],varOpts['-i'],varOpts['-m'],varOpts['-n'],varOpts['-p'],varOpts['-q'],varOpts['-r'],varOpts['-v'],varOpts['-x'],varOpts['-z'])
      commands += "touch %s;\n"%(markerFile2)
  jobs.append(commands)


if ( jobIndex == None ) :
  #multiple threads are available
  with Pool( nCPU ) as p:
    p.map(runCommand, jobs)
else :
  print( "Running command:" )
  print( jobs[ jobIndex -1] )
  runCommand( jobs[ jobIndex - 1 ] )

if ( noPost ) :
  print( "\nCompleted analysis with no-postprocessing option set\n")
  sys.exit(0)

#when using job indexes only last job must keep running after this
if ( jobIndex != None ) :
  if ( nTask2Run == 0 and jobIndex != 1 ) :
    sys.exit(0)
  if ( len(glob.glob("%s/tables/tmp/*.dsa.done"%(workDir)) ) + len(glob.glob("%s/variants/tmp/*.var.done"%(workDir)) )  != 2* nMarkerFiles ) :
    sys.exit(0)


#merge the variantcall results, generate summary

print( "\nGenerating summary\n" )
command = "set -e; variantcaller_merge.py -d %s/variants"%workDir
runCommand( command )

command = "set -e; cd %s;variantcaller.R ./variants > %s/variants/summary.txt"%(workDir,workDir)
try :
  runCommand( command )
except :
  pass

for ifile in glob.glob("%s/variants/*.csv"%workDir) : shutil.move(ifile,workDir+"/summary/"+os.path.basename(ifile))
for ifile in glob.glob("%s/variants/summary.txt"%workDir) : shutil.move(ifile,workDir+"/summary/"+os.path.basename(ifile))
for ifile in glob.glob("%s/variants/*.pdf"%workDir) : shutil.move(ifile,workDir+"/summary/"+os.path.basename(ifile))

command = "set -e; cd %s;botseq_results_plotter.R . results"%(workDir+'/summary')
try :
  runCommand( command )
except :
  pass

if ( os.path.exists( "%s/summary/variants.csv"%(workDir) ) ) : shutil.copy( "%s/summary/variants.csv"%(workDir), workDir )
if ( os.path.exists( "%s/summary/mismatches.csv"%(workDir) ) ) : shutil.copy( "%s/summary/mismatches.csv"%(workDir), workDir )
if ( os.path.exists( "%s/summary/results.muts.vcf"%(workDir) ) ) : shutil.copy( "%s/summary/results.muts.vcf"%(workDir) , workDir )

print( "\nArchiving results\n" )
#archive summary directories
command = "set -e;tar -cz -C %s -f %s/summary.tar.gz summary"%(workDir, workDir)
runCommand( command )
shutil.rmtree("%s/summary"%(workDir))

if ( archive ) :
  command = "set -e;tar -cz -C %s -f %s/variants.tar.gz variants"%(workDir, workDir)
  runCommand( command )
  shutil.rmtree("%s/variants"%(workDir))

  command = "set -e;tar -c  -C %s -f %s/tables.tar tables"%(workDir, workDir)
  runCommand( command )
  shutil.rmtree("%s/tables"%(workDir))

