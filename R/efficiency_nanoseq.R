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

suppressPackageStartupMessages({
  library(VGAM)
  library(Biostrings)
  library(Rsamtools)
  library(GenomicRanges)
  library(stringr)
  library(seqinr)
})

args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("efficiency_nanoseq.R rb_file genome\n\nMust specify a read-bundle file and a genome.\n\n")
  quit(save = "no", status = 0)
}

if (length(args) != 2) {
  cat("efficiency_nanoseq.R rb_file genome\n\nMust specify a read-bundle file and a genome.\n\n")
  quit(save = "no", status = 1)
}

rb_file = args[1]
genomeFile = args[2]

if (!file.exists(genomeFile)) {
  stop("Reference file not found : ", genomeFile, call. = FALSE)
}

if (!file.exists(rb_file)) {
  stop("readbundle file not found : ", rb_file, call. = FALSE)
}

pdf(paste(rb_file, ".pdf", sep = ""), width = 6, height = 6)
par(mar = c(2, 2, 2, 2))
cat("Reading file:", rb_file, "...\n")
rbs = read.table(rb_file, sep = "\t", header = F, row.names = 1)

rbs$x_tmp = pmin(rbs[, 1], 10)
rbs$y_tmp = pmin(rbs[, 2], 10)
rbs$x = rbs[,1]
rbs$y = rbs[,2]
rbs_bck = rbs;
kk = as.data.frame(table(rbs[, c("x", "y")]))
kk$point_size = kk$Freq / max(kk$Freq)
reads_per_RB = sum(c(rbs$x_tmp, rbs$y_tmp)) / nrow(rbs)
plot(as.numeric(kk$x),
       as.numeric(kk$y),
       pch = 20,
#cex=kk$point_size*10+.5,
       cex = kk$point_size * 10,
       xlab = "",
       ylab = "",
       main = paste(rb_file, ":", round(reads_per_RB, 3), sep = ""))
# Count how many do we have with size >= 4
# Define sizes 4, 5..., >=10
# For each size estimate how many we expect to have one strad ( 12.5% for 4, 6.25 for 5...)
rbs$size = rbs$x + rbs$y
rbs$size = pmin(rbs$size, 10)
total_missed = 0
for (size in c(4:10)) {
  exp_orphan = (0.5 ** size) * 2
  total_this_size = nrow(rbs[which(rbs$size == size),])
  if (total_this_size > 0) {
    with_both_strands = nrow(rbs[which(rbs$size == size & rbs$x > 0 & rbs$y > 0),])
    obs_orphan = 1 - with_both_strands / total_this_size
    missed = (obs_orphan - exp_orphan) * total_this_size
    total_missed = total_missed + missed
  }
}
total_missed_fraction = total_missed / nrow(rbs[which(rbs$size >= 4),])
rbs = rbs[which(rbs$size >= 2),]
rbs = rbs[sample(1:nrow(rbs), min(100000, nrow(rbs))),]
rbs$min = apply(rbs[, c("x", "y")], 1, min)
rbs$max = apply(rbs[, c("x", "y")], 1, max)
rbs$y = rbs$min
rbs$x = rbs$max
dev.off()

cat("READS_PER_RB\t", reads_per_RB, "\n", sep = "")
cat("F-EFF\t", total_missed_fraction, "\n", sep = "")
cat("OK_RBS\t", nrow(rbs_bck[which(rbs_bck$x >= 2 & rbs_bck$y >= 2),]), "\n", sep = "")
cat("TOTAL_RBS\t", nrow(rbs_bck), "\n", sep = "")
cat("TOTAL_READS\t", sum(c(rbs_bck$x, rbs_bck$y)) * 2, "\n", sep = "") #by 2 because we look only at one of the mates

##########################################################################################
# Now GC content:
#rbs = rbs_bck

colnames(rbs)[1:2] = c("plus", "minus")
rbs$id = rownames(rbs)
kk = str_split_fixed(rbs$id, "[:,]", 7)
rbs$chr = kk[, 3]
rbs$start = as.numeric(kk[, 4])
rbs$end = as.numeric(kk[, 5])

# Remove RBs with end breakpoints beyond chr/contig size
chr_coords = as.data.frame(scanFaIndex(genomeFile))
rownames(chr_coords) = chr_coords$seqnames
rbs$chr_end = chr_coords[rbs$chr, "end"]
rbs = rbs[which(rbs$end < rbs$chr_end),]

rbs_both = rbs[which(rbs$minus + rbs$plus >= 4 & rbs$minus >= 2 & rbs$plus >= 2),]
rbs_both = rbs_both[sample(1:nrow(rbs_both), min(10000, nrow(rbs_both))),] # 1000 random

rbs_single = rbs[which(rbs$minus + rbs$plus > 4 & (rbs$minus == 0 | rbs$plus == 0)),]
rbs_single = rbs_single[sample(1:nrow(rbs_single), min(10000, nrow(rbs_single))),] # 1000 random

seqs_both = as.vector(scanFa(genomeFile, GRanges(rbs_both$chr, IRanges(start = rbs_both$start, end = rbs_both$end))))
seqs_single = as.vector(scanFa(genomeFile, GRanges(rbs_single$chr, IRanges(start = rbs_single$start, end = rbs_single$end))))

seqs_both_collapsed = paste(seqs_both, collapse = "")
seqs_single_collapsed = paste(seqs_single, collapse = "")

tri_both = trinucleotideFrequency(DNAString(seqs_both_collapsed))
tri_single = trinucleotideFrequency(DNAString(seqs_single_collapsed))

tri_both_freqs = tri_both / sum(tri_both)
tri_single_freqs = tri_single / sum(tri_single)
if ( seqs_both_collapsed == "" ) {
  gc_both = 0
} else {
  gc_both = GC(s2c(seqs_both_collapsed))
}
if ( seqs_single_collapsed == "" ) {
  gc_single = 0
} else {
  gc_single = GC(s2c(seqs_single_collapsed))
}

cat("GC_BOTH\t", gc_both, "\n", sep = "")
cat("GC_SINGLE\t", gc_single, "\n", sep = "")

##########################################################################################
# RB sizes vs GC content & insert sizes:
rbs = rbs_bck


colnames(rbs)[1:2] = c("plus", "minus")
rbs$id = rownames(rbs)
rbs$size = rbs$plus + rbs$minus
kk = str_split_fixed(rbs$id, "[:,]", 7)
rbs$chr = kk[, 3]
rbs$start = as.numeric(kk[, 4])
rbs$end = as.numeric(kk[, 5])
rbs$insert = rbs$end - rbs$start

# Remove RBs with end breakpoints beyond chr/contig size
chr_coords = as.data.frame(scanFaIndex(genomeFile))
rownames(chr_coords) = chr_coords$seqnames
rbs$chr_end = chr_coords[rbs$chr, "end"]
rbs = rbs[which(rbs$end < rbs$chr_end),]


res = as.data.frame(matrix(nrow = 7, ncol = 3))
colnames(res) = c("size", "GC", "insert_len")
rownames(res) = c(1:7)
res[] = NA
for (size in c(1:7)) {
  rbs_ = rbs[which(rbs$size == size),]
  if (nrow(rbs_) == 0) {
    next
  }
  res[size, "size"] = size
  res[size, "insert_len"] = mean(rbs_$end - rbs_$start)
  rbs_ = rbs_[sample(1:nrow(rbs_), min(10000, nrow(rbs_))),] # 10000 random

  seqs_ = as.vector(scanFa(genomeFile, GRanges(rbs_$chr, IRanges(start = rbs_$start, end = rbs_$end))))

  seqs_collapsed = paste(seqs_, collapse = "")

  tri_ = trinucleotideFrequency(DNAString(seqs_collapsed))

  tri_freqs = tri_ / sum(tri_)

  gc_ = GC(s2c(seqs_collapsed))
  res[size, "GC"] = gc_

}
write.table(res, file = paste(rb_file, ".GC_inserts.tsv", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)


