#!/usr/bin/env Rscript

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

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# builds the following objects
#   burdens
#   callsvsqpos
#   coverage
#   pyrvsmask
#   readbundles
#   variants
#   unique_variants
#   unique_strand_context_vs_qpos
#   unique_pyr_var_vs_mask
#   metrics


# directory name from args
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  cat("Directory name argument must be supplied\n")
  quit(save="no",status=0)
}

dirname <- args[1]

if (! dir.exists(dirname)){
  stop("Input directory not found : ", dirname, call.=FALSE)
}

# remove trailing forward slash from dirname
dirname <- sub("/$", "", dirname)

# read in csv files
cat(getwd())
cat("\n")
cat( paste(dirname, 'burdens.csv', sep="/") )
cat("\n")
burdens     <- fread(paste(dirname, 'burdens.csv', sep="/"))
callsvsqpos <- fread(paste(dirname, 'callvsqpos.csv', sep="/"))
coverage    <- fread(paste(dirname, 'coverage.csv', sep="/"))
pyrvsmask   <- fread(paste(dirname, 'pyrvsmask.csv', sep="/"))
readbundles <- fread(paste(dirname, 'readbundles.csv', sep="/"))
variants    <- fread(paste(dirname, 'variants.csv', sep="/"))
mismatches  <- fread(paste(dirname, 'mismatches.csv', sep="/"))


# combine burdens from merged bed files
burdens <- burdens[, .(count=sum(count)), by=.(ismasked, isvariant)]

# combine call vs qpos from merged bed files
callsvsqpos <- callsvsqpos[, .(count=sum(count)), by=.(base, qpos, ismasked)]

# combine coverage from merged bed files
coverage <- sum(as.double(coverage$count))

# combine pyrimidine vs mask from merged bed files
pyrvsmask <- pyrvsmask[, .(count=sum(count)), by=.(pyrcontext, ismasked)]

# combine read bundles from merged bed files
readbundles <- readbundles[, .(count=sum(count)), by=.(fwd, rev, ismasked,
isvariant)]


# unique variants
#   unique defined by (chromosome, coordinate and substitution)
variants[, pyrkey:=paste(chrom, chromStart, pyrsub, sep=","),
         by=seq_len(nrow(variants))]
setkey(variants, pyrkey)
unique_variants <- unique(variants)

# print(dim(variants[ismasked == 0]))
# print(dim(unique_variants[ismasked == 0]))

# strand context versus query position
#   used to study asymmetries
#   uses unique variants
unique_strand_context_vs_qpos <- unique_variants[, .(.N), by = .(stdcontext,
  qpos, ismasked)]


# pyrimidine variant versus mask
#   used for trinucleotide profiles and cosine similarities
#   uses unique variants
unique_pyr_var_vs_mask <- unique_variants[, .(.N), by = .(pyrcontext,
  ismasked)]


# mismatches
# mismatches <- mismatches[, .(.N), by = .(mismatch, left, right)]]
#mismatches$mismatch2 <- paste(substring(mismatches$mismatch, 2, 2), substring(mismatches$mismatch, 4, 5), sep="")
#mismatches <- mismatches[mismatch != "UNKNOWN"][left > 2][right > 2][, .(.N), by = .(mismatch2)]
#print (mismatches, n = 1000)


# summary metrics
#   unmasked
n_variants  <- burdens[ismasked == 0][isvariant == 1]$count
n_reference <- burdens[ismasked == 0][isvariant == 0]$count
n_unique    <- nrow(unique_variants[ismasked == 0])

metrics <- data.frame("Metric" = c("total variants", "unique variants", "reference",
  "total variant + reference", "uncorrected burden", "physical coverage"),
  "Value" = c(format(n_variants, scientific = FALSE),
              format(n_unique, scientific = FALSE),
              format(n_reference, scientific = FALSE),
              format((n_variants + n_reference), scientific = FALSE),
              (n_variants/(n_variants + n_reference)), coverage))

print(metrics)


# strand 
# unique_variants$stdsub2 <- paste(substring(unique_variants$stdsub, 2, 2), substring(unique_variants$stdsub, 4, 5), sep="")
# unique_variants$qpos2 <- cut(unique_variants$qpos, breaks=seq(lower=0, upper=200, by = 25))
# stdsub_by_qpos_unmasked <- unique_variants[ismasked == 0][, .(.N), by = .(stdsub2, qpos2)][order(qpos2,stdsub2)]

# print(stdsub_by_qpos_unmasked[stdsub2 == "A>C"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "T>G"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "A>G"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "T>C"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "A>T"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "T>A"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "C>A"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "G>T"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "C>G"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "G>C"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "C>T"])
# print(stdsub_by_qpos_unmasked[stdsub2 == "G>A"])


# trinucleotide plot
unique_variants$pyrsubtype <- paste(substring(unique_variants$pyrsub, 2, 2),
  substring(unique_variants$pyrsub, 4, 5), sep="")

unique_variants$pyrcontexttype <- paste(substring(unique_variants$pyrcontext,
  1, 1), substring(unique_variants$pyrcontext, 3, 3), sep=".")

vars <- unique_variants[ismasked == 0
                      ][, .(.N), by = .(pyrcontexttype, pyrsubtype)]

COLORS6 = c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")

p1 <- ggplot(vars, aes(x=pyrcontexttype, y=N)) + 
  geom_bar(stat='identity', aes(fill=factor(pyrsubtype))) +
  facet_grid(. ~ pyrsubtype) +
  theme(panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "none") + 
  theme(axis.title=element_blank()) +
  scale_fill_manual(values=COLORS6)

ggsave(paste(dirname, "/triNucleotide.pdf", sep=""))




