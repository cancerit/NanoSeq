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
import copy

version = '3.5.7'

parser = argparse.ArgumentParser()
# arguments for all subcommands
parserO = parser._action_groups.pop()  # need this to get opt/req args correctly
# need this to get opt/req args correctly
parserR = parser.add_argument_group('required arguments')
parserO.add_argument('--out',   action='store', default='.',
                     help='path of the output files and scratch directory (.)')
parserO.add_argument('-j', '--index', type=int, action='store',
                     help='index of the LSF job array. One based')
parserO.add_argument('-k', '--max_index', type=int,
                     action='store', help='maximum index of the LSF job array')
parserO.add_argument('-t', '--threads', type=int,
                     action='store', default=1, help='number of threads (1)')
parserR.add_argument('-R', '--ref', action='store',
                     required=True, help="referene sequence")
parserR.add_argument('-A', '--normal', action='store',
                     required=True, help="normal BAM / CRAM")
parserR.add_argument('-B', '--duplex', action='store',
                     required=True, help="duplex (tumour) BAM / CRAM")
parserO.add_argument('-v', '--version', action='version', version=version)
# need this to get opt/req args correctly
parser._action_groups.append(parserO)

# subcommands and their arguments
# coverage
subparsers = parser.add_subparsers(dest='subcommand', help='subcommands')
subparsers.required = True  # work around for older python versions
parser_cov = subparsers.add_parser('cov', help='coverage calculation')
parser_covO = parser_cov._action_groups.pop()
parser_covR = parser_cov.add_argument_group('required arguments')
parser_covO.add_argument('--exclude', action='store', default='MT,GL%%,NC_%,hs37d5',
                         help='List of contigs to exclude. Comma separated, %% acts as a wild card. (MT,GL%%,NC_%%,hs37d5)')
parser_covO.add_argument('--include', action='store',
                         help='Only include these contigs. Comma separated, %% acts as a wild card.')
parser_covO.add_argument('--larger', type=int, action='store',
                         default=0, help='Only include contigs larger than this size. (0)')
parser_covO.add_argument('-w', '--win', type=int, action='store',
                         default=100, help='bin size for coverage distribution (100)')
parser_covO.add_argument('-Q', type=int, action='store',
                         default=0, help="minimum mapQ to include a duplex read (0)")
parser_cov._action_groups.append(parser_covO)

parser_part = subparsers.add_parser(
    'part', help='partition intervals into n jobs using coverage information')
parser_partO = parser_part._action_groups.pop()
parser_partR = parser_part.add_argument_group('required arguments')
parser_partR.add_argument('-n', '--jobs', type=int, action='store',
                          required=True, help='partition dsa,var,indel to this many tasks')
parser_partO.add_argument('--excludeBED', action='store',
                          help='BED (gz) file with regions to exclude from analysis.')
parser_partO.add_argument('--excludeCov', type=int, action='store',
                          help='Exclude regions with coverage values higher than this')
parser_part._action_groups.append(parser_partO)

# dsa
parser_dsa = subparsers.add_parser('dsa', help='compute tables')
parser_dsaO = parser_dsa._action_groups.pop()
parser_dsaR = parser_dsa.add_argument_group('required arguments')
parser_dsaR.add_argument('-C', '--snp', action='store',
                         help="SNP BED (gz) file")
parser_dsaR.add_argument('-D', '--mask', action='store',
                         help="mask BED (gz) file")
parser_dsaO.add_argument('-d', type=int, action='store',
                         default=2, help="minimum duplex depth (2)")
parser_dsaO.add_argument('-q', type=int, action='store',
                         default=30, help="minimum base quality for normal (30)")
parser_dsaO.add_argument('--no_test', action='store_true',
                         help="skip BAM format tests, use with caution")
parser_dsa._action_groups.append(parser_dsaO)

# variant caller
parser_var = subparsers.add_parser('var', help='variant caller')
parser_varO = parser_var._action_groups.pop()
parser_varR = parser_var.add_argument_group('required arguments')
parser_varO.add_argument('-a', type=int, action='store',
                         default=50, help="minimum AS-XS (50)")
parser_varO.add_argument('-b', type=int, action='store', default=0,
                         help="minimum matched normal reads per strand (0)")
parser_varO.add_argument('-c', type=float, action='store',
                         default=0.02, help="fraction of clips (0.02)")
parser_varO.add_argument('-d', type=int, action='store',
                         default=2, help="minimum duplex depth (2)")
parser_varO.add_argument('-f', type=float, action='store', default=0.9,
                         help="minimum fraction of reads for consensus (0.9)")
parser_varO.add_argument('-i', type=float, action='store', default=1.0,
                         help="maximum fraction of reads with an indel (1.0)")
parser_varO.add_argument('-m', type=int, action='store',
                         default=8, help="minimum cycle number (8)")
parser_varO.add_argument('-n', type=int, action='store',
                         default=3, help="maximun number of mismatches (3)")
parser_varO.add_argument('-p', type=int, action='store', default=0,
                         help="minimum fraction of reads that are proper-pairs (0)")
parser_varO.add_argument('-q', type=int, action='store',
                         default=45, help="minimum consensus base quality (45)")
parser_varO.add_argument('-r', type=int, action='store',
                         default=144, help="read length (after 5' trimming) (144)")
parser_varO.add_argument('-v', type=float, action='store',
                         default=0.01, help="maximum normal VAF (0.01)")
parser_varO.add_argument('-x', type=int, action='store',
                         default=8, help="maximum cycle number (8)")
parser_varO.add_argument('-z', type=int, action='store',
                         default=15, help="minimum normal coverage (15)")
parser_var._action_groups.append(parser_varO)

# compute indels
parser_indel = subparsers.add_parser('indel', help='indel caller')
parser_indelO = parser_indel._action_groups.pop()
parser_indelR = parser_indel.add_argument_group('required arguments')
parser_indelO.add_argument('-s', '--sample', action='store',
                           default='sample_1', help="sample name in output vcf (sample_1)")
parser_indelO.add_argument('--rb', type=int, action='store',
                           default=2, help="minimum reads in a bundle. (2)")
parser_indelO.add_argument('--t3', type=int, action='store', default=136,
                           help="excess bases above this value are trimmed from 3' (136)")
parser_indelO.add_argument('--t5', type=int, action='store',
                           default=8, help="bases to trim from 5' reads (8)")
parser_indelO.add_argument('-a', type=int, action='store',
                           default=50, help="minimum AS-XS (50)")
parser_indelO.add_argument('-c', type=float, action='store',
                           default=0.02, help="fraction of clips (0.02)")
parser_indelO.add_argument(
    '-z', type=int, action='store', default=15, help="minimum normal coverage (mc) (15)")
parser_indelO.add_argument(
    '-v', type=float, action='store', default=0.01, help="maximum normal VAF (0.01)")
parser_indel._action_groups.append(parser_indelO)

# carry out gather operations and compute summaries
parser_post = subparsers.add_parser(
    'post', help='gather final files, compute summaries')
parser_postO = parser_post._action_groups.pop()
parser_postR = parser_post.add_argument_group('required arguments')
parser_postO.add_argument('--name', action='store',
                          default='results', help="name for output files (results)")
parser_postO.add_argument('--triNuc', action='store',
                          help="tri-nucleotide correction file")
parser_post._action_groups.append(parser_postO)

args = parser.parse_args()

# check job partition arguments
if ((args.index is None and args.max_index is not None) or
        (args.index is not None and args.max_index is None)):
    parser.error("Must specify index and max_index for array execution!")

if (args.index is None and args.max_index is None):
    jobArray = False
else:
    jobArray = True

if (jobArray and args.threads > 1):
    print("Warning: execution with a job array is single threaded so thread argument is ignored")

# check the location of the tmpNanoSeq dir
if not os.path.isdir(args.out):
    parser.error("Specified out directory %s is not accessible!" % args.out)

try:
    testfile = tempfile.TemporaryFile(dir=args.out)
    testfile.close()
except OSError:
    sys.exit("\nCan't write to out directory %s\n" % args.out)

# check required files


def file_chk(fn, idx_ext, msg_prefix):
    if not os.path.isfile(fn):
        parser.error(f"{msg_prefix} file {fn} was not found!")
    if not os.path.isfile(fn + idx_ext):
        parser.error(f"{msg_prefix} index file {fn}{idx_ext} was not found!")


if (hasattr(args, 'duplex')):
    ext = os.path.splitext(args.duplex)[1][0:-1] + "i"
    file_chk(args.duplex, ext, "BAM/CRAM")

if (hasattr(args, 'normal')):
    ext = os.path.splitext(args.normal)[1][0:-1] + "i"
    file_chk(args.normal, ext, "BAM/CRAM")

if (hasattr(args, 'ref')):
    file_chk(args.ref, ".fai", "Reference")

if (hasattr(args, 'snp')):
    if ( args.snp is not None ) :
      file_chk(args.snp, ".tbi", "SNP")

if (hasattr(args, 'mask')):
    if ( args.mask is not None ) :
      file_chk(args.mask, ".tbi", "Mask")

if (hasattr(args, 'triNuc') and args.triNuc is not None):
    if not os.path.isfile(args.triNuc):
        parser.error(
            "Trinucleotide correction file %s was not found!" % (args.triNuc))

# check that all the dependancies for the analysis are in PATH
scripts = ["Rscript",  # indel, post
           "bcftools",  # indel
           "samtools",  # indel
           "bamcov",  # cov
           "dsa",  # dsa
           "variantcaller",  # var
           "tabix",
           "bgzip",
           "variantcaller.R",  # post
           "indelCaller_step1.pl",  # indel
           "indelCaller_step2.pl",  # indel
           "indelCaller_step3.R",  # indel
           "efficiency_nanoseq.pl",  # post
           "efficiency_nanoseq.R"  # post
           ]

for icode in scripts:
    if (shutil.which(icode) is None):
        raise ValueError("%s was not found in path!" % icode)

# handle genomic intervals with a class


class GInterval:
    def __init__(self, chrr, beg, end):
        self.chr = chrr
        if (end < beg):
            raise ValueError("Interval %s: %s - %s is invalid!" %
                             (chrr, beg, end))
        self.beg = beg
        self.end = end
        self.l = end - beg + 1

    def convert2DSAInput(self):
        # zero based and inclusive of end
        return(GInterval(self.chr, self.beg-1, self.end-1))

    def __repr__(self):
        return repr((self.chr, self.beg, self.end))

    def __str__(self):
        return "%s:%s-%s" % (self.chr, self.beg, self.end)

    def __contains__(self, other):
        if (self.chr != other.chr):
            return False
        if (other.end < self.beg or self.end < other.beg):
            return False
        return True

    def __sub__(self, other):
        if (self.chr != other.chr):
            return [self]
        if (other.beg > self.beg and other.beg <= self.end):
            if (other.end < self.end):
                return [GInterval(self.chr, self.beg, other.beg-1), GInterval(self.chr, other.end+1, self.end)]
            if (other.end >= self.end):
                return [GInterval(self.chr, self.beg, other.beg - 1)]
        if (other.beg <= self.beg):
            if (other.end < self.end and other.end >= self.beg):
                return [GInterval(self.chr, other.end + 1, self.end)]
            if (other.end >= self.end):
                return []
        return [self]

    def __add__(self, other):
        if (other.chr != self.chr):
            return [self, other]
        if (other.end + 1 == self.beg or self.end + 1 == other.beg):
            return [GInterval(self.chr, min([self.beg, other.beg]), max([self.end, other.end]))]
        if (other.end < self.beg or self.end < other.beg):
            return [self, other]
        if (self.beg <= other.end or other.beg <= self.end):
            return [GInterval(self.chr, min([self.beg, other.beg]), max([self.end, other.end]))]

    # define to allow sorting of intervals

    # sort order for chromosomes
    # Compares two chromosome names with the order being -
    #     1. Simple names before long names
    #     2. Names with 'chr' then number (or just bare number) before Names with 'chr' then letter(s)
    #     3. 'chrX', 'X', 'chrY', 'Y', 'M', 'MT' after 2.
    #     4. Others
    def chr_isless(x, y):
        if (x == y):
            return False

        xc = re.sub('^chr', '', x)
        yc = re.sub('^chr', '', y)

        if (xc.isdigit()):
            if (yc.isdigit()):
                return int(xc) < int(yc)
            else:
                return True  # rule2
        elif (yc.isdigit()):
            return False  # rule2
        elif ((xc == "M" or xc == "MT") and (yc == "X" or yc == "Y")):
            return False  # rule 3
        elif ((xc == "X" or xc == "Y") and (yc == "M" or yc == "MT")):
            return True  # rule 3
        elif (xc == "X"):
            return True  # rule 3
        elif (yc == "X"):
            return False  # rule 3
        elif (xc == "Y"):
            return True  # rule 3
        elif (yc == "Y"):
            return False  # rule 3
        elif ((xc == "M" or xc == "MT")):
            return True  # rule 3
        elif ((yc == "M" or yc == "MT")):
            return False  # rule 3
        else:
            return xc < yc

    def __eq__(self, other):
        if (self.chr == other.chr):
            if (self.beg == other.beg):
                if (self.end == other.end):
                    return True
        return False

    def __ne__(self, other):
        if (self.chr == other.chr):
            if (self.beg == other.beg):
                return False
        return True

    def __gt__(self, other):
        if (GInterval.chr_isless(other.chr, self.chr)):
            return True
        if (GInterval.chr_isless(self.chr, other.chr)):
            return False
        if (self.beg > other.beg):
            return True
        else:
            return False

    def __ge__(self, other):
        if (GInterval.chr_isless(other.chr, self.chr)):
            return True
        if (GInterval.chr_isless(self.chr, other.chr)):
            return False
        if (self.beg >= other.beg):
            return True
        else:
            return False

    def __lt__(self, other):
        if (GInterval.chr_isless(self.chr, other.chr)):
            return True
        if (GInterval.chr_isless(other.chr, self.chr)):
            return False
        if (self.beg < other.beg):
            return True
        else:
            return False

    def __le__(self, other):
        if (GInterval.chr_isless(self.chr, other.chr)):
            return True
        if (GInterval.chr_isless(other.chr, self.chr)):
            return False
        if (self.beg <= other.beg):
            return True
        else:
            return False


def runCommand(command):
    if (command is None):
        return
    for ijob in command.rstrip(';').split(';'):
        print("\nExecuting: %s\n" % ijob)
        p = subprocess.Popen(
            ijob, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        std: tuple = p.communicate()
        if (p.returncode != 0):
            error = std[1].decode()
            sys.stderr.write("\n!Error processing:  %s\n" % ijob)
            raise ValueError(error)
    return

# compute coverage histogram


def runBamcov(bam, mapQ, window, ichr, out):
    if (bam is None):
        return
    out = out+".cov.bed"
    runCommand("bamcov -q %s -w %s -r \"%s\" -o %s %s" %
               (mapQ, window, ichr, out, bam))
    runCommand("bgzip -l 2 -f %s" % (out))
    outdone = re.sub('cov.bed$', 'done', out)
    open(outdone, 'w').close()
    return


def getBAMcontigs(bam):
    contigs = []
    if (bam is None):
        return contigs
    job = "samtools view -H %s" % bam
    p = subprocess.Popen(
        job, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        output = p.stdout.readline().decode()
        if output == '' and p.poll() is not None:
            break
        if (re.match("@SQ", output)):
            contigs.append(output.strip())
    rc = p.poll()
    if (rc != 0):
        error = p.stderr.read().decode()
        sys.stderr.write("\n!Error processing:  %s\n" % job)
        raise ValueError(error)
    return contigs


def vcfHeader(args):
    header = '##fileformat=VCFv4.2\n'
    header += '##source=NanoSeq pipeline\n'
    header += '##FILTER=<ID=PASS,Description="All filters passed">\n'
    header += "##reference=file://%s\n" % (args.ref)

    contigs = []
    with open(args.ref + ".fai", 'r') as iofile:
        for iline in iofile:
            ichr = iline.split('\t')[0]
            ilength = int(iline.split('\t')[1])
            contigs.append(GInterval(ichr, 1, ilength))
    contigs.sort()
    for ii in (contigs):
        ichr = ii.chr
        ilength = ii.end
        header += "##contig=<ID=%s,length=%s>\n" % (ichr, ilength)
    header += '##INFO=<ID=BTAG,Number=1,Type=String,Description="Read bundle tag (duplex barcode)">\n'
    header += '##INFO=<ID=BBEG,Number=1,Type=String,Description="Read bundle left breakpoint">\n'
    header += '##INFO=<ID=BEND,Number=1,Type=String,Description="Read bundle right breakpoint">\n'
    header += '##INFO=<ID=TRI,Number=1,Type=String,Description="Pyrimidine context, trinucleotide substitution">\n'
    header += '##INFO=<ID=QPOS,Number=1,Type=Integer,Description="Read position closest to 5-prime end">\n'
    header += '##INFO=<ID=DEPTH_FWD,Number=1,Type=Integer,Description="Read bundle forward reads depth">\n'
    header += '##INFO=<ID=DEPTH_REV,Number=1,Type=Integer,Description="Read bundle reverse reads depth">\n'
    header += '##INFO=<ID=DEPTH_NORM_FWD,Number=1,Type=Integer,Description="Matched normal forward reads depth">\n'
    header += '##INFO=<ID=DEPTH_NORM_REV,Number=1,Type=Integer,Description="Matched normal reverse reads depth">\n'
    header += '##INFO=<ID=DPLX_ASXS,Number=1,Type=Integer,Description="AS-XS for duplex">\n'
    header += '##INFO=<ID=DPLX_CLIP,Number=1,Type=Integer,Description="Clipping for duplex">\n'
    header += '##INFO=<ID=DPLX_NM,Number=1,Type=Integer,Description="Mismatches in duplex">\n'
    header += '##INFO=<ID=BULK_ASXS,Number=1,Type=Integer,Description="AS-XS for bulk">\n'
    header += '##INFO=<ID=BULK_NM,Number=1,Type=Integer,Description="Mismatches in bulk">\n'
    header += '##FILTER=<ID=dbsnp,Description="Common SNP site">\n'
    header += '##FILTER=<ID=shearwater,Description="Noisy site">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    return header


print("command arguments :")

print(args.__dict__)
print()

# for array execution try to stagger access for files
if (args.index is not None):
    time.sleep( args.index)

# create all directory tree for the tmp files
tmpDir = args.out + "/tmpNanoSeq"
if (args.subcommand == 'cov'):
    for idir in ('cov', 'part', 'dsa', 'var', 'indel', 'post'):
        if (not os.path.isdir(tmpDir+'/'+idir)):
            os.makedirs(tmpDir+'/'+idir, exist_ok=True)

if (args.index is None or args.index == 1):
    with open("%s/%s/args.json" % (tmpDir, args.subcommand), "w") as jsonOut:
        json.dump(args.__dict__, jsonOut)

# coverage section
if (args.subcommand == 'cov'):
    # build the chromosome dictionary, list and intervals
    if (args.exclude is None or args.exclude == ""):
        excludes = []  # exlcude None
    else:
        excludes = [re.compile(istr + "$")
                    for istr in args.exclude.replace("%", ".+").split(',')]
    if (args.include is None or args.include == ""):
        includes = [re.compile(".+")]  # include all
    else:
        includes = [re.compile(istr + "$")
                    for istr in args.include.replace("%", ".+").split(',')]
    chrList = []
    rnames = {}
    with open(args.ref + '.fai', 'r') as iofile:
        for iline in iofile:
            ichr = iline.split('\t')[0]
            if any(iregx.match(ichr) for iregx in includes):
                if any(iregx.match(ichr) for iregx in excludes):
                    continue
                ilength = int(iline.split('\t')[1])
                if (ilength <= args.larger):
                    continue
                chrList.append(ichr)
                rnames[ichr] = ilength

    # BAMs vs reference sanity checks for contigs
    bamNContigs = {}
    bamNOrder = []
    for iline in getBAMcontigs(args.normal):
        ichr = iline.split('\t')[1].replace('SN:', '')
        ilength = int(iline.split('\t')[2].replace('LN:', ''))
        bamNContigs[ichr] = ilength
        bamNOrder.append(ichr)

    bamTContigs = {}
    bamTOrder = []
    for iline in getBAMcontigs(args.duplex):
        ichr = iline.split('\t')[1].replace('SN:', '')
        ilength = int(iline.split('\t')[2].replace('LN:', ''))
        bamTContigs[ichr] = ilength
        bamTOrder.append(ichr)

    for icontig in rnames:
        if (not (icontig in bamNContigs)):
            sys.exit("Reference contig %s was not found in normal BAM" % icontig)
        if (bamNContigs[icontig] != rnames[icontig]):
            sys.exit("Length of contig %s in normal BAM doesn't match reference (%s vs %s)" % (
                icontig, bamNContigs[icontig], rnames[icontig]))

    for icontig in rnames:
        if (not (icontig in bamTContigs)):
            sys.exit("Reference contig %s was not found in duplex BAM" % icontig)
        if (bamTContigs[icontig] != rnames[icontig]):
            sys.exit("Length of contig %s in duplex BAM doesn't match reference (%s vs %s)" % (
                icontig, bamTContigs[icontig], rnames[icontig]))

    # if order of BAMs is not the same mpilupe bogus results
    for icontig in rnames:
        if (bamNOrder.index(icontig) != bamTOrder.index(icontig)):
            sys.exit(
                "Contigs in BAM files must have the same order ( check order in headers )")

    gintervals = []
    for ichr in chrList:
        gintervals.append(GInterval(ichr, 2, rnames[ichr]-1))
    gintervals.sort()

    reorderchr = []
    for iint in gintervals:
        reorderchr.append(iint.chr)
    chrList = reorderchr
    print("\nAnalysing contigs: %s\n" % chrList)
    print("Starting cov calculation\n")
    if (args.index is None or args.index == 1):
        with open("%s/cov/nfiles" % (tmpDir), "w") as iofile:
            iofile.write(str(len(chrList)))
    inputs = []
    for ii, ichr in enumerate(chrList):
        if (os.path.isfile("%s/cov/%s.done" % (tmpDir, ii + 1))):
            # restart, don't do anything
            inputs.append((None, None, None, None, None))
        else:
            inputs.append((args.duplex, str(args.Q), str(
                args.win), ichr, "%s/cov/%s" % (tmpDir, ii + 1)))
    if (args.index is None):
        # multiple threads are available
        with open("%s/cov/%s" % (tmpDir, 'gIntervals.dat'), 'wb') as iofile:
            pickle.dump(gintervals, iofile)
        with Pool(args.threads) as p:
            p.starmap(runBamcov, inputs)
    else:
        # single thread execution
        if (args.index == 1):
            with open("%s/cov/%s" % (tmpDir, 'gIntervals.dat'), 'wb') as iofile:
                pickle.dump(gintervals, iofile)

        commands = [[] for i in range(args.max_index)]
        jj = 0
        for ii in range(len(inputs)):
            if (ii % args.max_index == 0):
                jj = 0
            commands[jj].append(
                "runBamcov(inputs[%s][0], inputs[%s][1], inputs[%s][2], inputs[%s][3], inputs[%s][4])" % (ii, ii, ii, ii, ii))
            jj += 1
        for icmd in commands[args.index - 1]:
            exec(icmd)

    print("Completed cov calculation\n")

# merge coverage files, partition into desired number of jobs
if (args.subcommand == 'part'):
    if (os.path.isfile("%s/part/1.done" % (tmpDir))):
        exit(0)  # restart
    if (not os.path.isfile(tmpDir+'/cov/args.json')):
        sys.exit("\nMust run cov submmand prior to part\n")
    else:
        with open(tmpDir+'/cov/args.json') as iofile:
            oargs = json.load(iofile)
    if (not os.path.isfile(tmpDir+'/cov/nfiles')):
        sys.exit(tmpDir+'/cov/nfiles not found!\n')
    else:
        with open(tmpDir+'/cov/nfiles') as iofile:
            nfiles = int(iofile.readline())
    for i in range(nfiles):
        if (len(glob.glob(tmpDir+"/cov/%s.done" % (i+1))) != 1):
            sys.exit("\ncov job %s did not complete correctly\n" % (i+1))
        if (len(glob.glob(tmpDir+"/cov/%s.cov.bed.gz" % (i+1))) != 1):
            sys.exit("\ncov job %s did not complete correctly\n" % (i+1))

    if (args.index is not None and args.index > 1):
        print("\nWarning can only use 1 job of array\n")
        sys.exit(0)

    with open(tmpDir+"/cov/gIntervals.dat", 'rb') as iofile:
        gIntervals = pickle.load(iofile)

    coverage = []
    cctotal = 0
    chrOffset = {}
    tmpIntervals = []
    print("\nParsing coverage files\n")
    for i in range(nfiles):
        with gzip.open(tmpDir+"/cov/%s.cov.bed.gz" % (i+1), 'rt') as iofile:
            for iline in iofile:
                ichr = str(iline.split('\t')[0])
                ib = int(iline.split('\t')[1])
                ie = int(iline.split('\t')[2])
                cc = int(iline.split('\t')[3])
                cctotal += cc
                if (args.excludeCov is not None):
                    if (cc >= args.excludeCov):
                        tmpIntervals.append(GInterval(ichr, ib+1, ie))
                if (ib == 0):
                    chrOffset[str(ichr)] = len(coverage)
                coverage.append([ib, cc])
    print("\nCompleted parsing coverage files\n")

    # remove regions to ignore from exclude BED
    if (args.excludeBED is not None):
        with gzip.open(args.excludeBED, 'rt') as iofile:
            for iline in iofile:
                ichr = str(iline.split('\t')[0])
                ib = int(iline.split('\t')[1])
                ie = int(iline.split('\t')[2])
                # convert to ginterval
                tmpIntervals.append(GInterval(ichr, ib+1, ie))
        tmpIntervals.sort()

    if (len(tmpIntervals) > 0):
        # merge overlapping intervals
        xIntervals = [tmpIntervals.pop(0)]
        while (len(tmpIntervals) > 0):
            xIntervals.extend(xIntervals.pop() + tmpIntervals.pop(0))

        print("\nExcluding %s intervals, dumping to BED file : %s \n" %
              (len(xIntervals), "%s/part/exclude.bed" % tmpDir))
        if (args.index is None or args.index == 1):
            with open("%s/part/%s" % (tmpDir, 'exclude.bed'), 'w') as iofile:
                for ii in xIntervals:
                    iofile.write("%s\t%s\t%s\n" % (ii.chr, ii.beg-1, ii.end))
            cmd = "bgzip -f %s/part/%s; sleep 3; bgzip -t %s/part/%s; tabix %s/part/%s" % (
                tmpDir, 'exclude.bed', tmpDir, 'exclude.bed.gz', tmpDir, 'exclude.bed.gz')
            runCommand(cmd)

        # remove the excluded intervals
        iiresult = []
        for ii in gIntervals:
            ifrag = ii
            for jj in xIntervals:
                if (ifrag.chr != jj.chr):
                    continue
                diff = ifrag - jj
                if len(diff) == 2:
                    iiresult.append(diff[0])
                    ifrag = diff[1]
                elif len(diff) == 1:
                    ifrag = diff[0]
                else:
                    break
            else:
                iiresult.append(ifrag)
        gIntervals = iiresult

        # correct total coverage (cctotal)
        xSumCov = 0
        for iinterval in gIntervals:
            ichar = iinterval.chr
            ibeg = iinterval.beg - 1
            iend = iinterval.end - 1
            for i in range(math.floor(ibeg/oargs['win']), math.floor(iend/oargs['win']) + 1):
                if ( not ichar in chrOffset ) : break
                j = i + chrOffset[ichar]
                xSumCov += coverage[j][1]
        cctotal = xSumCov

    # determine the genomic intervals to give to each job so that each one roughly
    # has the same amount of coverage
    basesPerCPU = 0
    njobs = args.jobs
    basesPerCPU = cctotal / njobs
    print("\nPartitioning %s jobs with %s bases/task\n" % (njobs, basesPerCPU))

    sumCov = 0
    oIntervals = []
    intervalsPerCPU = []
    gIntervalsCopy = copy.deepcopy(gIntervals)
    while (len(gIntervals) > 0):
        iinterval = gIntervals.pop(0)
        ichar = iinterval.chr
        ibeg = iinterval.beg - 1
        iend = iinterval.end - 1
        for i in range(math.floor(ibeg/oargs['win']), math.floor(iend/oargs['win']) + 1):
            j = i + chrOffset[ichar]
            sumCov += coverage[j][1]
            if (sumCov > basesPerCPU):
                jend = min([coverage[j][0] + oargs['win'], iend])
                oIntervals.append(GInterval(ichar, ibeg + 1, jend + 1))
                intervalsPerCPU.append(oIntervals)
                oIntervals = []
                sumCov = 0
                ibeg = jend + 1
        if (iend >= ibeg):
            oIntervals.append(GInterval(ichar, ibeg + 1, iend + 1))
    if (len(oIntervals) > 0):
        intervalsPerCPU.append(oIntervals)
    while ( len(intervalsPerCPU) < njobs ) :
        intervalsPerCPU.append([])
    # check partitioning code is working as expected
    # compare merged partitioned intervals to original gIntervals (copy)
    print("\nChecking partition of intervals..", end='')
    flatInt = [item for sublist in intervalsPerCPU[0:njobs] for item in sublist]
    nbases1 = 0
    for i in flatInt:
        nbases1 += i.l
    nbases2 = 0
    for i in gIntervalsCopy:
        nbases2 += i.l
    if (nbases1 != nbases2):
        sys.exit("Internal check failed: interval lenght after partition doesent match origianl (%s, %s)\n"%(nbases2, nbases1))
    print("..", end='')
    mIntervals = [flatInt.pop(0)]
    while (len(flatInt) > 0):
        mIntervals.extend(mIntervals.pop() + flatInt.pop(0))
    for (i, ival) in enumerate(gIntervalsCopy):
        if (not (ival == mIntervals[i])):
            sys.exit("Internal check failed: mismatch after part (interval %s should be %s)\n" % (
                mIntervals[i], ival))
    print(" OK\n")
    with open("%s/part/%s" % (tmpDir, 'intervalsPerCPU.dat'), 'wb') as iofile:
        pickle.dump(intervalsPerCPU, iofile)
    cmd = "touch %s/part/1.done" % (tmpDir)
    runCommand(cmd)
    print("\nCompleted part job\n")

# dsa section
if (args.subcommand == 'dsa'):
    if (not os.path.isfile(tmpDir+'/part/args.json')):
        sys.exit("\nMust run cov and part submmands prior to dsa\n")
    else:
        with open(tmpDir+'/part/args.json') as iofile:
            oargs = json.load(iofile)
    njobs = oargs['jobs']

    if (len(glob.glob(tmpDir+"/part/1.done")) != 1):
        sys.exit("\npart job did not complete correctly\n")
    if (len(glob.glob(tmpDir+"/part/intervalsPerCPU.dat")) != 1):
        sys.exit("\npart job did not complete correctly\n")

    # make sure that number of jobs matches what was specified in part
    if (args.max_index is not None):
        # array execution
        if (args.max_index < njobs):
            sys.exit(
                "\nLSF array size must match number of jobs specified for part (%s)\n" % njobs)
        if (args.index > njobs):
            print(
                "\nWarning specified LSF array size is larger than jobs specified for part (%s)\n" % njobs)
            sys.exit(0)
    else:
        # multithread
        if (args.threads < njobs):
            sys.exit(
                "\nNumber of threads must match number of jobs specified for part (%s)\n" % njobs)
        if (args.threads > njobs):
            print(
                "\nWarning number of threads is larger than jobs specified for part (%s)\n" % njobs)

    with open(tmpDir+"/part/intervalsPerCPU.dat", 'rb') as iofile:
        intervalsPerCPU = pickle.load(iofile)

    mapQ = None
    with open(tmpDir+'/cov/args.json') as iofile:
        mapQ = json.load(iofile)['Q']

    commands = [(None, )] * njobs
    for i in range(njobs):
        # check for restarts
        if (os.path.isfile("%s/dsa/%s.done" % (tmpDir, i+1)) and
                os.path.isfile("%s/dsa/%s.dsa.bed.gz" % (tmpDir, i+1))):
            continue

        # construct dsa commands
        cmd = ""
        testOpt = ""
        snpOpt = ""
        maskOpt = ""
        if (args.no_test):
            testOpt = "-t"
        if ( args.snp is not None ) :
            snpOpt = "-C %s"%args.snp
        if ( args.mask is not None ) :
            maskOpt = "-D %s"%args.mask

        for (ii, iinterval) in enumerate(intervalsPerCPU[i]):
            dsaInt = iinterval.convert2DSAInput()
            pipe = ">" if ii == 0 else ">>"  # ensure first command overwrittes
            cmd += "dsa -A %s -B %s  %s %s -R %s -d %s -Q %s -M %s %s -r \"%s\" -b %s -e %s %s %s ;" \
                % (args.normal, args.duplex, snpOpt, maskOpt, args.ref, args.d, args.q, mapQ, testOpt,
                   dsaInt.chr, dsaInt.beg, dsaInt.end, pipe, "%s/dsa/%s.dsa.bed" % (tmpDir, i + 1))
        # check number of fields in the last line it has to have 45 fields
        cmd += "awk  \'END{  if (NF != 45)  print \"Truncated dsa output file for job %s !\" > \"/dev/stderr\"}{ if (NF != 45) exit 1 }\' %s/dsa/%s.dsa.bed;" % (i+1, tmpDir, i+1)
        cmd += "bgzip -f -l 2 %s/dsa/%s.dsa.bed; sleep 2; bgzip -t %s/dsa/%s.dsa.bed.gz;" % (
            tmpDir, i+1, tmpDir, i+1)
        cmd += "touch %s/dsa/%s.done" % (tmpDir, i+1)
        if ( len(intervalsPerCPU[i]) == 0 ) : cmd = "touch %s/dsa/%s.dsa.bed.gz;touch %s/dsa/%s.done" % (tmpDir, i+1,tmpDir, i+1)
        commands[i] = (cmd, )

    if (args.index is None or args.index == 1):
        with open("%s/dsa/nfiles" % (tmpDir), "w") as iofile:
            iofile.write(str(njobs))

    # execute dsa commans
    print("Starting dsa calculation\n")
    if (args.index is None):
        # multithread
        with Pool(args.threads) as p:
            p.starmap(runCommand, commands)
    else:
        # array execution
        runCommand(commands[args.index - 1][0])
    print("Completed dsa calculation\n")

# var section
if (args.subcommand == 'var'):
    if (not os.path.isfile(tmpDir+'/dsa/args.json')):
        sys.exit("\nMust run dsa submmand prior to var\n")
    else:
        with open(tmpDir+'/dsa/args.json') as iofile:
            oargs = json.load(iofile)

    if (not os.path.isfile(tmpDir+'/dsa/nfiles')):
        sys.exit(tmpDir+'/dsa/nfiles not found!\n')
    else:
        with open(tmpDir+'/dsa/nfiles') as iofile:
            nfiles = int(iofile.readline())
    njobs = nfiles

    for i in range(nfiles):
        if (len(glob.glob(tmpDir+"/dsa/%s.done" % (i+1))) != 1):
            sys.exit("\ndsa job %s did not complete correctly\n" % (i+1))
        if (len(glob.glob(tmpDir+"/dsa/%s.dsa.bed.gz" % (i+1))) != 1):
            sys.exit("\ndsa job %s did not complete correctly\n" % (i+1))

    # make sure that number of jobs matches what was specified in part
    if (args.max_index is not None):
        # array execution
        if (args.max_index < njobs):
            sys.exit(
                "\nLSF array size must match number of jobs specified for part (%s)\n" % njobs)
        if (args.index > njobs):
            print(
                "\nWarning specified LSF array size is larger than jobs specified for part (%s)\n" % njobs)
            sys.exit(0)
    else:
        # multithread
        if (args.threads < njobs):
            sys.exit(
                "\nNumber of threads must match number of jobs specified for part (%s)\n" % njobs)
        if (args.threads > njobs):
            print(
                "\nWarning number of threads is larger than jobs specified for part (%s)\n" % njobs)

    commands = [(None,)] * njobs
    for i in range(njobs):
        # check for restarts
        if (os.path.isfile("%s/var/%s.done" % (tmpDir, i+1)) and
            os.path.isfile("%s/var/%s.var" % (tmpDir, i+1)) and
                os.path.isfile("%s/var/%s.cov.bed.gz" % (tmpDir, i+1))):
            continue

        # construct variantcaller commands
        cmd = "variantcaller -B %s -U %s -O %s -D %s -a %s -b %s -c %s -d %s -f %s -i %s -m %s -n %s -p %s -q %s -r %s -v %s -x %s -z %s ;" \
            % ("%s/dsa/%s.dsa.bed.gz" % (tmpDir, i+1), "%s/var/%s.cov.bed" % (tmpDir, i+1), "%s/var/%s.var" % (tmpDir, i+1), "%s/var/%s.discarded_var" % (tmpDir, i+1),
               args.a, args.b, args.c, args.d, args.f, args.i, args.m, args.n, args.p, args.q, args.r, args.v, args.x, args.z)
        cmd += "touch %s/var/%s.done" % (tmpDir, i+1)
        if ( os.stat("%s/dsa/%s.dsa.bed.gz"%(tmpDir, i+1)).st_size == 0) : cmd = "touch %s/var/%s.var; touch %s/var/%s.discarded_var; touch %s/var/%s.cov.bed.gz; touch %s/var/%s.done" % (tmpDir, i+1,tmpDir, i+1,tmpDir, i+1,tmpDir, i+1)
        commands[i] = (cmd, )

    if (args.index is None or args.index == 1):
        with open("%s/var/nfiles" % (tmpDir), "w") as iofile:
            iofile.write(str(njobs))

    # execute variantcaller commans
    print("Starting var calculation\n")
    if (args.index is None):
        # multithread
        with Pool(args.threads) as p:
            p.starmap(runCommand, commands)
    else:
        # array execution
        runCommand(commands[args.index - 1][0])
    print("Completed var calculation\n")

# indel section
if (args.subcommand == 'indel'):
    if (not os.path.isfile(tmpDir+'/dsa/args.json')):
        sys.exit("\nMust run dsa submmand prior to indel\n")
    else:
        with open(tmpDir+'/dsa/args.json') as iofile:
            oargs = json.load(iofile)

    if (not os.path.isfile(tmpDir+'/dsa/nfiles')):
        sys.exit(tmpDir+'/dsa/nfiles not found!\n')
    else:
        with open(tmpDir+'/dsa/nfiles') as iofile:
            nfiles = int(iofile.readline())
    njobs = nfiles

    for i in range(nfiles):
        if (len(glob.glob(tmpDir+"/dsa/%s.done" % (i+1))) != 1):
            sys.exit("\ndsa job %s did not complete correctly\n" % (i+1))
        if (len(glob.glob(tmpDir+"/dsa/%s.dsa.bed.gz" % (i+1))) != 1):
            sys.exit("\ndsa job %s did not complete correctly\n" % (i+1))

    # make sure that number of jobs matches what was specified in part
    if (args.max_index is not None):
        # array execution
        if (args.max_index < njobs):
            sys.exit(
                "\nLSF array size must match number of jobs specified for part (%s)\n" % njobs)
        if (args.index > njobs):
            print(
                "\nWarning specified LSF array size is larger than jobs specified for part (%s)\n" % njobs)
            sys.exit(0)
    else:
        # multithread
        if (args.threads < njobs):
            sys.exit(
                "\nNumber of threads must match number of jobs specified for part (%s)\n" % njobs)
        if (args.threads > njobs):
            print(
                "\nWarning number of threads is larger than jobs specified for part (%s)\n" % njobs)

    commands = [(None, )] * njobs
    for i in range(njobs):
        # check for restarts
        if (os.path.isfile("%s/indel/%s.done" % (tmpDir, i+1)) and
                os.path.isfile("%s/indel/%s.indel.filtered.vcf.gz" % (tmpDir, i+1))):
            continue

        # construct the indel commands ( 3 steps)
        cmd = "indelCaller_step1.pl -o %s -rb %s -t3 %s -t5 %s -mc %s -vaf %s -a %s -c %s %s ;"\
            % ("%s/indel/%s.indel.bed.gz" % (tmpDir, i+1), args.rb, args.t3, args.t5, args.z, args.v, args.a, args.c,
               "%s/dsa/%s.dsa.bed.gz" % (tmpDir, i+1))
        cmd += "indelCaller_step2.pl -t -o %s -r %s -b %s %s ;"\
            % ("%s/indel/%s.indel" % (tmpDir, i+1), args.ref, args.duplex,
               "%s/indel/%s.indel.bed.gz" % (tmpDir, i+1))
        cmd += "indelCaller_step3.R %s %s %s %s;"\
            % (args.ref, "%s/indel/%s.indel.vcf.gz" % (tmpDir, i+1), args.normal, args.v)
        cmd += "touch %s/indel/%s.done" % (tmpDir, i+1)
        if ( os.stat("%s/dsa/%s.dsa.bed.gz"%(tmpDir, i+1)).st_size == 0) :
            cmd = "touch %s/indel/%s.indel.bed.gz; touch %s/indel/%s.indel.vcf.gz;" \
                "touch %s/indel/%s.indel.filtered.vcf.gz; touch %s/indel/%s.indel.filtered.vcf.gz.tbi;" \
                "touch %s/indel/%s.done" % (tmpDir, i+1,tmpDir, i+1,tmpDir, i+1,tmpDir, i+1,tmpDir, i+1)
        commands[i] = (cmd, )

    if (args.index is None or args.index == 1):
        with open("%s/indel/nfiles" % (tmpDir), "w") as iofile:
            iofile.write(str(njobs))

    # execute indeliantcaller commans
    print("Starting indel calculation\n")
    if (args.index is None):
        # multithread
        with Pool(args.threads) as p:
            p.starmap(runCommand, commands)
    else:
        # array execution
        runCommand(commands[args.index - 1][0])
    print("Completed indel calculation\n")

# post section
if (args.subcommand == 'post'):
    if (os.path.isfile(tmpDir+'/post/1.done')):
        sys.exit(0)
    if (args.index is not None and args.index > 1):
        print("\nWarning can only use 1 job of array\n")
        sys.exit(0)

    # check dsa
    if (not os.path.isfile(tmpDir+'/dsa/nfiles')):
        sys.exit(tmpDir+'/dsa/nfiles not found!\n')
    else:
        with open(tmpDir+'/dsa/nfiles') as iofile:
            nfiles = int(iofile.readline())

    for i in range(nfiles):
        if (len(glob.glob(tmpDir+"/dsa/%s.done" % (i+1))) != 1):
            sys.exit("\ndsa job %s did not complete correctly\n" % (i+1))
        if (len(glob.glob(tmpDir+"/dsa/%s.dsa.bed.gz" % (i+1))) != 1):
            sys.exit("\ndsa job %s did not complete correctly\n" % (i+1))

    # check var
    did_var = False
    if (os.path.isfile(tmpDir+'/var/nfiles')):
        did_var = True
        with open(tmpDir+'/dsa/nfiles') as iofile:
            nfiles = int(iofile.readline())
        for i in range(nfiles):
            if (len(glob.glob(tmpDir+"/var/%s.done" % (i+1))) != 1):
                sys.exit("\nvar job %s did not complete correctly\n" % (i+1))
            if (len(glob.glob(tmpDir+"/var/%s.var" % (i+1))) != 1):
                sys.exit("\nvar job %s did not complete correctly\n" % (i+1))
            if (len(glob.glob(tmpDir+"/var/%s.discarded_var" % (i+1))) != 1):
                sys.exit("\nvar job %s did not complete correctly\n" % (i+1))
            if (len(glob.glob(tmpDir+"/var/%s.cov.bed.gz" % (i+1))) != 1):
                sys.exit("\nvar job %s did not complete correctly\n" % (i+1))

    # check indel
    did_indel = False
    if (os.path.isfile(tmpDir+'/indel/nfiles')):
        did_indel = True
        with open(tmpDir+'/indel/nfiles') as iofile:
            nfiles = int(iofile.readline())
        for i in range(nfiles):
            if (len(glob.glob(tmpDir+"/indel/%s.done" % (i+1))) != 1):
                sys.exit("\nindel job %s did not complete correctly\n" % (i+1))
            if (len(glob.glob(tmpDir+"/indel/%s.indel.filtered.vcf.gz" % (i+1))) != 1):
                sys.exit("\nindel job %s did not complete correctly\n" % (i+1))

    if (did_var):
        # generate csv files
        print("\nGenerating CSV files for var\n")
        csvFiles = ['Coverage', 'CallVsQpos', 'PyrVsMask',
                    'ReadBundles', 'Burdens', 'Variants', 'DiscardedVariants', 'Mismatches']
        csvIO = {}
        for ifile in csvFiles:
            csvIO[ifile] = open('%s/post/%s.csv' %
                                (tmpDir, ifile.lower()), 'w')

        # write headers
        csvIO['Coverage'].write('count\n')
        csvIO['CallVsQpos'].write('base,qpos,ismasked,count\n')
        csvIO['PyrVsMask'].write('pyrcontext,ismasked,count\n')
        csvIO['ReadBundles'].write('fwd,rev,ismasked,isvariant,count\n')
        csvIO['Burdens'].write('ismasked,isvariant,count\n')
        csvIO['Variants'].write('chrom,chromStart,context,commonSNP,'
                                'shearwater,bulkASXS,bulkNM,bulkForwardA,bulkForwardC,bulkForwardG,'
                                'bulkForwardT,bulkForwardIndel,bulkReverseA,bulkReverseC,bulkReverseG,'
                                'bulkReverseT,bulkReverseIndel,dplxBreakpointBeg,dplxBreakpointEnd,'
                                'bundleType,dplxASXS,dplxCLIP,dplxNM,dplxfwdA,dplxfwdC,dplxfwdG,dplxfwdT,'
                                'dplxfwdIndel,dplxrevA,dplxrevC,dplxrevG,dplxrevT,dplxrevIndel,'
                                'dplxCQfwdA,dplxCQfwdC,dplxCQfwdG,dplxCQfwdT,dplxCQrevA,'
                                'dplxCQrevC,dplxCQrevG,dplxCQrevT,bulkForwardTotal,bulkReverseTotal,'
                                'dplxfwdTotal,dplxrevTotal,left,right,qpos,call,isvariant,pyrcontext,'
                                'stdcontext,pyrsub,stdsub,ismasked,dplxBarcode\n')
        csvIO['DiscardedVariants'].write('chrom,chromStart,context,commonSNP,'
                                'shearwater,bulkASXS,bulkNM,bulkForwardA,bulkForwardC,bulkForwardG,'
                                'bulkForwardT,bulkForwardIndel,bulkReverseA,bulkReverseC,bulkReverseG,'
                                'bulkReverseT,bulkReverseIndel,dplxBreakpointBeg,dplxBreakpointEnd,'
                                'bundleType,dplxASXS,dplxCLIP,dplxNM,dplxfwdA,dplxfwdC,dplxfwdG,dplxfwdT,'
                                'dplxfwdIndel,dplxrevA,dplxrevC,dplxrevG,dplxrevT,dplxrevIndel,'
                                'dplxCQfwdA,dplxCQfwdC,dplxCQfwdG,dplxCQfwdT,dplxCQrevA,'
                                'dplxCQrevC,dplxCQrevG,dplxCQrevT,bulkForwardTotal,bulkReverseTotal,'
                                'dplxfwdTotal,dplxrevTotal,left,right,qpos,call,isvariant,pyrcontext,'
                                'stdcontext,pyrsub,stdsub,ismasked,dplxBarcode,'
                                'dplx_clip_filter,alignment_score_filter,mismatch_filter,matched_normal_filter,'
                                'duplex_filter,consensus_base_quality_filter,indel_filter,five_prime_trim_filter,'
                                'three_prime_trim_filter,proper_pair_filter,vaf_filter\n')
        csvIO['Mismatches'].write('chrom,chromStart,context,commonSNP,'
                                  'shearwater,bulkASXS,bulkNM,bulkForwardA,bulkForwardC,bulkForwardG,'
                                  'bulkForwardT,bulkForwardIndel,bulkReverseA,bulkReverseC,bulkReverseG,'
                                  'bulkReverseT,bulkReverseIndel,dplxBreakpointBeg,dplxBreakpointEnd,'
                                  'dplxASXS,dplxCLIP,dplxNM,dplxfwdA,dplxfwdC,dplxfwdG,dplxfwdT,'
                                  'dplxfwdIndel,dplxrevA,dplxrevC,dplxrevG,dplxrevT,dplxrevIndel,'
                                  'dplxCQfwdA,dplxCQfwdC,dplxCQfwdG,dplxCQfwdT,dplxCQrevA,'
                                  'dplxCQrevC,dplxCQrevG,dplxCQrevT,bulkForwardTotal,bulkReverseTotal,'
                                  'dplxfwdTotal,dplxrevTotal,left,right,qpos,mismatch,ismasked,dplxBarcode\n')

        # wirte body
        for i in range(nfiles):
            # Write variants
            ifile = "%s/var/%s.var" % (tmpDir, i+1)
            for row in open(ifile, 'r'):
                if (row[0] == '#'):
                    continue
                arow = row.strip().split('\t')
                if csvIO.get(arow[0], None):
                    csvIO[arow[0]].write('%s\n' % ','.join(arow[1:]))
            
            # Write discarded variants
            dfile = "%s/var/%s.discarded_var" % (tmpDir, i+1)
            for row in open(dfile, 'r'):
                if (row[0] == '#'):
                    continue
                arow = row.strip().split('\t')
                if csvIO.get(arow[0], None):
                    csvIO[arow[0]].write('%s\n' % ','.join(arow[1:]))


        for ifile in csvIO.values():
            ifile.close()

        print("\nMerge coverage files for var\n")
        # merge coverage files
        outFile = "%s/post/%s.cov.bed" % (tmpDir, args.name)
        cmd = "rm -f %s;" % (outFile)
        for i in range(nfiles):
            ifile = "%s/var/%s.cov.bed.gz" % (tmpDir, i+1)
            if ( os.stat(ifile).st_size == 0 ) : continue
            cmd += "bgzip -dc %s >> %s ;" % (ifile, outFile)
        cmd += "bgzip -@ %s -f %s; sleep 3; bgzip -@ %s -t %s.gz ;" % (
            args.threads, outFile, args.threads, outFile)
        cmd += "tabix -f %s.gz" % (outFile)
        runCommand(cmd)

        print("\nCompute summaries for var\n")
        # do the summary
        cmd = "variantcaller.R %s/post/ > %s/post/summary.txt" % (
            tmpDir, tmpDir)
        runCommand(cmd)

        cmd = "nanoseq_results_plotter.R %s/post %s/post/%s %s" % (
            tmpDir, tmpDir, args.name, args.triNuc or "")
        runCommand(cmd)

        print("\nGenerate vcf file from variants.csv\n")

        var = {}
        nVariants = 0
        with open("%s/post/variants.csv" % tmpDir, 'r') as iofile:
            iline = iofile.readline().rstrip('\n')
            fields = iline.split(',')
            for ifield in fields:
                var[ifield] = []
            for iline in iofile:
                nVariants += 1
                for (i, ival) in enumerate(iline.rstrip('\n').split(',')):
                    var[fields[i]].append(ival)

        #added so code will not break with older files
        if ( len(var['dplxBarcode']) == 0 ) : var['dplxBarcode'] = nVariants * [ '' ]

        header = vcfHeader(args)

        with open("%s/post/%s.muts.vcf" % (tmpDir, args.name), "w") as iofile:
            iofile.write(header)
            for i in range(nVariants):
                iline = "%s\t%s\t%s\t%s\t%s\t.\t" % \
                    (var['chrom'][i], int(var['chromStart'][i]) +
                     1, '.', var['context'][i][1], var['call'][i])
                ifilter = "PASS"
                if (var['shearwater'][i] == '1'):
                    ifilter = "shearwater"
                if (var['commonSNP'][i] == '1'):
                    ifilter = "dbsnp"
                iline += "%s\t" % ifilter
                iline += "BTAG=%s;BBEG=%s;BEND=%s;TRI=%s;QPOS=%s;DEPTH_FWD=%s;DEPTH_REV=%s;DEPTH_NORM_FWD=%s;DEPTH_NORM_REV=%s;DPLX_ASXS=%s;DPLX_CLIP=%s;DPLX_NM=%s;BULK_ASXS=%s;BULK_NM=%s\n" % \
                    (var['dplxBarcode'][i], var['dplxBreakpointBeg'][i], var['dplxBreakpointEnd'][i], var['pyrsub'][i],
                     var['qpos'][i], var['dplxfwdTotal'][i], var['dplxrevTotal'][i], var['bulkForwardTotal'][i],
                     var['bulkReverseTotal'][i], var['dplxASXS'][i], var['dplxCLIP'][i], var['dplxNM'][i], var['bulkASXS'][i],
                     var['bulkNM'][i])
                iofile.write(iline)
        cmd = "bgzip -@ %s -f %s/post/%s.muts.vcf; sleep 3; bgzip -@ %s -t %s/post/%s.muts.vcf.gz;" % (
            args.threads, tmpDir, args.name, args.threads, tmpDir, args.name)
        cmd += "bcftools index -t -f %s/post/%s.muts.vcf.gz " % (
            tmpDir, args.name)
        runCommand(cmd)

    if (did_indel):
        print("\nMerging vcf files for indel\n")
        vcf2Merge = []
        for i in range(nfiles):
            ifile = tmpDir+"/indel/%s.indel.filtered.vcf.gz" % (i+1)
            if ( os.stat(ifile).st_size == 0 ) : continue
            vcf2Merge.append(ifile)
        cmd = "cp %s %s/indel/merged.vcf.gz ;" % (vcf2Merge.pop(0), tmpDir)
        for ifile in vcf2Merge[0:]:
            cmd += "bcftools concat --no-version -Oz -o %s/indel/tmp.xxx.vcf.gz %s/indel/merged.vcf.gz %s; mv %s/indel/tmp.xxx.vcf.gz %s/indel/merged.vcf.gz; " % (tmpDir, tmpDir, ifile, tmpDir, tmpDir)
        cmd += "bcftools sort -Oz -o %s/post/%s.indel.vcf.gz %s/indel/merged.vcf.gz;" % (tmpDir, args.name, tmpDir)
        cmd += "bcftools index -t -f %s/post/%s.indel.vcf.gz " % (
            tmpDir, args.name)
        runCommand(cmd)

    cmd = "touch %s/post/1.done" % (tmpDir)
    runCommand(cmd)
