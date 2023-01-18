#!/usr/bin/env Rscript

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

#    Rscript to check whether each of the potential indel calls
#    happen at loci rich in indels in the matched normal. Those would be flagged as 'NEI_IND'.
#    Reliable indels are flagged with 'PASS'.
#    Three parameters are expected:
#      1. Reference genome
#      2. The name of the output file from step 2 of indelCaller
#      3. The BAM file containing the matched normal

suppressPackageStartupMessages({
  library(deepSNV)
  library(vcfR)
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
})
options(warn=2) #turn warnings into errors

args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0 ) {
  cat("indelCaller_step3.R  reference  vcf  bam  max_vaf\n\n Check identified indels against matched normal.\n\n Must specify: reference, VCF from indelCaller_step2, BAM/CRAM file for the matched normal\n\n")
  quit(save = "no", status = 0)
}

if (length(args) != 4) {
  cat("indelCaller_step3.R  reference  vcf  bam  max_vaf\n\n Check identified indels against matched normal.\n\n Must specify: reference, VCF from indelCaller_step2, BAM/CRAM file for the matched normal\n\n")
  quit(save = "no", status = 1)
}

genomeFile = args[1]
vcf_file = args[2]
bam_file = args[3]
max_vaf  = args[4]

if (!file.exists(genomeFile)) {
  stop("Reference file not found : ", genomeFile, call. = FALSE)
}
if (!file.exists(vcf_file)) {
  stop("VCF file not found : ", vcf_file, call. = FALSE)
}
if (!file.exists(bam_file)) {
  stop("Matched BAM (CRAM) file not found : ", bam_file, call. = FALSE)
}

if (length(grep("\\.gz", vcf_file)) > 0) {
  system(paste("gzip -t ", vcf_file), intern= TRUE)
}

FLANK = 5

vcf <- read.vcfR(vcf_file, verbose = FALSE)
out_vcf_file_tmp = gsub(".vcf", ".filtered.tmp.vcf", vcf_file)
out_vcf_file = gsub(".vcf", ".filtered.vcf", vcf_file)

# Create regions:
#  * For deletions: get pos + length(deletion)   // length(deletion) = length(ref)-1
#  * For insertions: get pos
#  * Then, sum +/-5 to each side to count indels in the vicinity
cat("processing: ", vcf_file, "\n")
if (nrow(vcf@fix) != 0) {
  for (i in c(1:nrow(vcf@fix))) {
    pos = strtoi(vcf@fix[i, "POS"])
    chr = vcf@fix[i, "CHROM"]
    len = max(1, length(vcf@fix[i, "REF"]) - 1)
    start = pos - FLANK
    end = pos + len + FLANK
    kk = bam2R(bam_file, chr, start, end, q = -100, mask = 3844, mq = 10)
    # dont want a filter on BQ because in some bams BQ of indels have -1
    n_bases = sum(kk[, c("A", "C", "G", "T", "a", "c", "g", "t")])
    n_indels = sum(kk[, c("-", "INS", "_", "ins")]) # Number of reads with an indel around the mutation
    cat(chr, pos, len, n_bases, n_indels, "\n");
    cat(n_bases, "/", n_indels, "\n")
    cat(i, "\n")
    max_per_site = max(apply(kk[, c("-", "INS", "_", "ins")],1,sum) / (apply(kk[, c("-", "INS", "_", "ins")],1,sum) + apply(kk[, c("A", "C", "G", "T", "a", "c", "g", "t")],1,sum)))
    if (n_bases == 0) {
      vcf@fix[i, "FILTER"] = "MISSINGBULK"
      vcf@fix[i, "INFO"] = paste(vcf@fix[i, "INFO"], ";NN=[", n_indels, ":", n_bases, ":", max_per_site,"]", sep = "")

    #} else if (n_indels/(n_bases + n_indels) > max_vaf) {
    #  vcf@fix[i, "FILTER"] = "NEI_IND"
    #  vcf@fix[i, "INFO"] = paste(vcf@fix[i, "INFO"], ";NN=[", n_indels, "/", n_bases, "]", sep = "")
    } else if (max_per_site > max_vaf ) { 
      vcf@fix[i, "FILTER"] = "NEI_IND"
      vcf@fix[i, "INFO"] = paste(vcf@fix[i, "INFO"], ";NN=[", n_indels, ":", n_bases, ":", max_per_site,"]", sep = "")
    } else {
      vcf@fix[i, "INFO"] = paste(vcf@fix[i, "INFO"], ";NN=[", n_indels, ":", n_bases, ":", max_per_site,"]", sep = "")
    }
    sequence = as.vector(scanFa(genomeFile, GRanges(chr, IRanges(start - 3, end + 3))))
    vcf@fix[i, "INFO"] = paste(vcf@fix[i, "INFO"], ";SEQ=", sequence, sep = "")
  }

  #Add new fields to the header
  newhead = c()
  once = TRUE
  for (i in vcf@meta) {
    if (once & grepl("##FILTER=<ID=MASKED", i, fixed = TRUE)) {
      i = paste("##FILTER=<ID=MISSINGBULK,Description=\"Site was not found in the matched normal\">\n", i, sep = "")
      i = paste("##FILTER=<ID=NEI_IND,Description=\"Site was found in an indel rich region of the matched normal\">\n", i, sep = "")
      i = paste("##INFO=<ID=NN,Number=1,Type=String,Description=\"n indels / n bases\">\n", i, sep = "")
      i = paste("##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Sequence of indel plus flanking sequences\">\n", i, sep = "")
      once = FALSE
    }
    newhead = c(newhead, i)
  }
  vcf@meta = newhead
}
write.vcf(vcf, file = out_vcf_file_tmp)
bgzip(out_vcf_file_tmp, dest = out_vcf_file, overwrite = TRUE)
unlink(out_vcf_file_tmp)
indexTabix(out_vcf_file, format = "vcf")


