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
  library(data.table)
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
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Directory name argument must be supplied\n")
  quit(save = "no", status = 0)
}

if (length(args) != 1) {
  cat("Directory name argument must be supplied\n")
  quit(save = "no", status = 1)
}

dirname <- args[1]

if (!dir.exists(dirname)) {
  stop("Input directory not found : ", dirname, call. = FALSE)
}

# remove trailing forward slash from dirname
dirname <- sub("/$", "", dirname)

# read in csv files
burdens <- fread(paste(dirname, 'burdens.csv', sep = "/"))
callsvsqpos <- fread(paste(dirname, 'callvsqpos.csv', sep = "/"))
coverage <- fread(paste(dirname, 'coverage.csv', sep = "/"))
pyrvsmask <- fread(paste(dirname, 'pyrvsmask.csv', sep = "/"))
readbundles <- fread(paste(dirname, 'readbundles.csv', sep = "/"))
variants <- fread(paste(dirname, 'variants.csv', sep = "/"))
mismatches <- fread(paste(dirname, 'mismatches.csv', sep = "/"))


# combine burdens from merged bed files
burdens <- burdens[, .(count = sum(count)), by = .(ismasked, isvariant)]

# combine call vs qpos from merged bed files
callsvsqpos <- callsvsqpos[, .(count = sum(count)), by = .(base, qpos, ismasked)]

# combine coverage from merged bed files
coverage <- sum(as.double(coverage$count))

# combine pyrimidine vs mask from merged bed files
pyrvsmask <- pyrvsmask[, .(count = sum(count)), by = .(pyrcontext, ismasked)]

# combine read bundles from merged bed files
readbundles <- readbundles[, .(count = sum(count)), by = .(fwd, rev, ismasked,
isvariant)]

# unique variants
#   unique defined by (chromosome, coordinate and substitution)
variants[, pyrkey := paste(chrom, chromStart, pyrsub, sep = ","),
         by = seq_len(nrow(variants))]
setkey(variants, pyrkey)
unique_variants <- unique(variants)

# summary metrics
#   unmasked
n_variants <- burdens[ismasked == 0][isvariant == 1]$count
n_reference <- burdens[ismasked == 0][isvariant == 0]$count
n_unique <- nrow(unique_variants[ismasked == 0])

metrics <- data.frame("Metric" = c("total variants", "unique variants", "reference",
  "total variant + reference", "uncorrected burden", "physical coverage"),
  "Value" = c(format(n_variants, scientific = FALSE),
              format(n_unique, scientific = FALSE),
              format(n_reference, scientific = FALSE),
              format((n_variants + n_reference), scientific = FALSE),
              (n_variants / (n_variants + n_reference)), coverage))

print(metrics)
