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


#########################################################################################
# This scripts receives a SNV vcf, an indel vcf, the sample deduplicated bam and
# the xx.cov.bed.gz file and with all this creates a new VAF by:
# * merging subs into DNVs and MNVs, using RB barcodes to discern if they are independent 
#   events
# * collapse and annotate repeated mutations
# * calculate duplex & mut VAFs for subs and indels
##########################################################################################

suppressPackageStartupMessages({
  library(deepSNV)
})
options(warn = 2) #turn warnings into errors

args = commandArgs(trailingOnly = TRUE)
muts_vcf = args[1]
indel_vcf = args[2]
dedup_bam = args[3] #neat NanoSeq BAM
cov_bed = args[4]
output_file = args[5]

if (length(args) == 0) {
  cat("Wrong number of arguments\n")
  cat("snv_merge_and_vaf_calc.R muts_vcf indel_vcf dedup_bam cov_bed output_file\n")
  quit(save = "no", status = 0)
}

if (length(args) != 5) {
  cat("Wrong number of arguments\n")
  cat("snv_merge_and_vaf_calc.R muts_vcf indel_vcf dedup_bam cov_bed output_file\n")
  quit(save = "no", status = 1)
}

if (!file.exists(muts_vcf)) {
  stop("SNV vcf file not found : ", muts_vcf, call. = FALSE)
}
if (!file.exists(indel_vcf)) {
  if (indel_vcf != '-') {
    stop("Indel vcf file not found : ", indel_vcf, call. = FALSE)
  }
}

if (!file.exists(dedup_bam)) {
  stop("Dedup BAM (neat) file not found : ", dedup_bam, call. = FALSE)
}
if (!file.exists(cov_bed)) {
  stop("Coverage BED file not found : ", cov_bed, call. = FALSE)
}

if (length(grep("\\.gz", muts_vcf)) > 0) {
  system(paste("gzip -t ", muts_vcf), intern = TRUE)
  num_snvs = system(paste("zgrep -cv \"^#\" ", muts_vcf, "|| true"), intern = TRUE)
} else {
  num_snvs = system(paste("grep -cv \"^#\" ", muts_vcf, "|| true"), intern = TRUE)
}
num_snvs = as.integer(num_snvs)

num_indels = 0
if (indel_vcf != '-') {
  if (length(grep("\\.gz", indel_vcf)) > 0) {
    system(paste("gzip -t ", indel_vcf), intern = TRUE)
    num_indels = system(paste("zgrep -cv \"^#\" ", indel_vcf, "|| true"), intern = TRUE)
  } else {
    num_indels = system(paste("grep -cv \"^#\" ", indel_vcf, "|| true"), intern = TRUE)
  }
}
num_indels = as.integer(num_indels)

col_t = c("character","numeric","character","character","character","character","character")
if (num_snvs == 0) {
  cat("No SNVs, skipping snv analysis...\n")
} else {
  if (length(grep("\\.gz", muts_vcf)) > 0) {
    zz = gzfile(muts_vcf, 'rt')
    snvs = read.table(zz, header = F, sep = "\t", colClasses = col_t, stringsAsFactors = F)
    close(zz)
  } else {
    snvs = read.table(muts_vcf, header = F, sep = "\t", colClasses = col_t, stringsAsFactors = F)
  }
  colnames(snvs)[1:8] = c("chr", "pos", "kk", "ref", "mut", "qual", "filter", "info")
}

col_t2 = c("character","numeric","character","character","character","numeric","character","character","character","character")
if (num_indels == 0) {
  cat("No Indels, skipping indel analysis...\n")
} else {
  if (length(grep("\\.gz", indel_vcf)) > 0) {
    zz = gzfile(indel_vcf, 'rt')
    indels = read.table(zz, header = F, sep = "\t", colClasses = col_t2, stringsAsFactors = F)
    close(zz)
  } else {
    indels = read.table(indel_vcf, header = F, sep = "\t", colClasses = col_t, stringsAsFactors = F)
  }
  colnames(indels)[1:8] = c("chr", "pos", "kk", "ref", "mut", "qual", "filter", "info")
}

#cat("Parsing...\n")
snvs_final = read.table(text = "", colClasses =
  col_t,
  col.names = c("chr", "pos", "kk", "ref", "mut", "qual", "filter", "INFO"))

xx = system(paste("gzip -t ", cov_bed), intern = TRUE)
tbx = TabixFile(cov_bed)

if (num_snvs > 0) {
  snvs$info = gsub("=", ";", snvs$info)
  mat = read.table(text = snvs$info, sep = ';', strip.white = TRUE, stringsAsFactors = F)
  for (i in seq(from = 1, to = ncol(mat) - 1, 2)) {
    col_name = mat[1, i]
    values = mat[, i + 1]
    snvs[, col_name] = values
  }
  if (length(intersect("BTAG", colnames(snvs))) == 0) {
    # for compatibility with pipeline v2
    snvs$BTAG = "KKK|KKK"
  }
  snvs$rb_id = paste(snvs$chr, snvs$BBEG, snvs$BEND, snvs$BTAG, sep = ":")
  snvs$rb_id = gsub("\\|", ":", snvs$rb_id)

  ##########################################################################################
  # Merge DNVs and MNVs
  # sort by rb_id and find consecutive muts
  #cat("Merging...\n")
  snvs = snvs[order(snvs$rb_id, snvs$chr, snvs$pos),]

  snvs$consecutive = 0
  counter = 1
  at_least_one = 0
  if (nrow(snvs) >= 2) {
    for (i in c(2:nrow(snvs))) {
      if (snvs[i - 1, "rb_id"] == snvs[i, "rb_id"] & snvs[i - 1, "pos"] + 1 == snvs[i, "pos"]) {
        snvs[i, "consecutive"] = counter
        snvs[i - 1, "consecutive"] = counter
        at_least_one = 1
      } else {
        counter = counter + 1
      }
    }
  }
  snvs_new = snvs[which(snvs$consecutive == 0),]
  snvs_new$TYPE = "snv"
  if (at_least_one == 1) {
    for (i in setdiff(unique(snvs$consecutive), 0)) {
      snvs_tmp = snvs[which(snvs$consecutive == i),]
      # merge snvs_tmp and add to snvs_new:
      new_row = nrow(snvs_new) + 1
      snvs_new[new_row, "chr"] = snvs_tmp[1, "chr"]
      snvs_new[new_row, "pos"] = snvs_tmp[1, "pos"]
      snvs_new[new_row, "kk"] = snvs_tmp[1, "kk"]
      snvs_new[new_row, "ref"] = paste(snvs_tmp[, "ref"], collapse = "")
      snvs_new[new_row, "mut"] = paste(snvs_tmp[, "mut"], collapse = "")
      snvs_new[new_row, "qual"] = snvs_tmp[1, "qual"]
      if (length(unique(snvs_tmp$filter)) > 1) {

        snvs_new[new_row, "filter"] = paste(snvs_tmp[, "filter"], collapse = ";") # replaced comman with semicolon (fa8)

      } else {
        snvs_new[new_row, "filter"] = snvs_tmp[1, "filter"]
      }
      snvs_new[new_row, "rb_id"] = snvs_tmp[1, "rb_id"]
      snvs_new[new_row, "BBEG"] = snvs_tmp[1, "BBEG"]
      snvs_new[new_row, "BEND"] = snvs_tmp[1, "BEND"]
      snvs_new[new_row, "TRI"] = NA
      snvs_new[new_row, "QPOS"] = snvs_tmp[1, "QPOS"]
      snvs_new[new_row, "DEPTH_FWD"] = sum(snvs_tmp[, "DEPTH_FWD"])
      snvs_new[new_row, "DEPTH_REV"] = sum(snvs_tmp[, "DEPTH_REV"])
      snvs_new[new_row, "DEPTH_NORM_FWD"] = sum(snvs_tmp[, "DEPTH_NORM_FWD"])
      snvs_new[new_row, "DEPTH_NORM_REV"] = sum(snvs_tmp[, "DEPTH_NORM_REV"])
      snvs_new[new_row, "DPLX_ASXS"] = snvs_tmp[1, "DPLX_ASXS"]
      snvs_new[new_row, "DPLX_CLIP"] = snvs_tmp[1, "DPLX_CLIP"]
      snvs_new[new_row, "DPLX_NM"] = snvs_tmp[1, "DPLX_NM"]
      snvs_new[new_row, "BULK_ASXS"] = snvs_tmp[1, "BULK_ASXS"]
      snvs_new[new_row, "BULK_NM"] = snvs_tmp[1, "BULK_NM"]         
      if (nrow(snvs_tmp) == 2) {
        snvs_new[new_row, "TYPE"] = "dnv"
      }
      if (nrow(snvs_tmp) > 2) {
        snvs_new[new_row, "TYPE"] = "mnv"
      }
    }
  }
  # if any sub is PASS, set as PASS (fa8)
  snvs_new[grep("PASS",snvs_new$filter),"filter"] = "PASS"
  # end
  snvs_new = snvs_new[order(snvs_new$chr, snvs_new$pos, snvs_new$rb_id),]

  ##########################################################################################
  # Create mutation ids and collapse repeated mutation calls
  #cat("Collapsing...\n")
  snvs_new$mut_id = paste(snvs_new$chr, snvs_new$pos, snvs_new$ref, snvs_new$mut, sep = ":")
  counts_per_mut = table(snvs_new$mut_id)
  repeat_muts = counts_per_mut[which(counts_per_mut > 1)]
  unique_muts = counts_per_mut[which(counts_per_mut == 1)]

  snvs_new2 = snvs_new[which(snvs_new$mut_id %in% names(unique_muts)),]
  if (nrow(snvs_new2) > 0) {
    snvs_new2$TIMES_CALLED = 1
  }
  if (length(repeat_muts) > 0) {
    for (mut_id in names(repeat_muts)) {
      new_row = nrow(snvs_new2) + 1
      snvs_tmp = snvs_new[which(snvs_new$mut_id == mut_id),]

      freq = nrow(snvs_tmp)
      snvs_new2[new_row, "chr"] = snvs_tmp[1, "chr"]
      snvs_new2[new_row, "pos"] = snvs_tmp[1, "pos"]
      snvs_new2[new_row, "kk"] = snvs_tmp[1, "kk"]
      snvs_new2[new_row, "ref"] = snvs_tmp[1, "ref"]
      snvs_new2[new_row, "mut"] = snvs_tmp[1, "mut"]
      snvs_new2[new_row, "qual"] = snvs_tmp[1, "qual"]
      snvs_new2[new_row, "filter"] = snvs_tmp[1, "filter"]
      snvs_new2[new_row, "rb_id"] = paste(snvs_tmp[, "rb_id"], collapse = ",")
      snvs_new2[new_row, "BBEG"] = paste(snvs_tmp[, "BBEG"], collapse = ",")
      snvs_new2[new_row, "BEND"] = paste(snvs_tmp[, "BEND"], collapse = ",")
      snvs_new2[new_row, "TRI"] = snvs_tmp[1, "TRI"]
      snvs_new2[new_row, "QPOS"] = paste(snvs_tmp[, "QPOS"], collapse = ",")
      snvs_new2[new_row, "DEPTH_FWD"] = median(as.numeric(snvs_tmp[, "DEPTH_FWD"]))
      snvs_new2[new_row, "DEPTH_REV"] = median(as.numeric(snvs_tmp[, "DEPTH_REV"]))
      snvs_new2[new_row, "DEPTH_NORM_FWD"] = snvs_tmp[1, "DEPTH_NORM_FWD"]
      snvs_new2[new_row, "DEPTH_NORM_REV"] = snvs_tmp[1, "DEPTH_NORM_REV"]
      snvs_new2[new_row, "TIMES_CALLED"] = freq
      snvs_new2[new_row, "DPLX_ASXS"] = paste(snvs_tmp[, "DPLX_ASXS"], collapse = ",")
      snvs_new2[new_row, "DPLX_CLIP"] = paste(snvs_tmp[, "DPLX_CLIP"], collapse = ",")
      snvs_new2[new_row, "DPLX_NM"] = paste(snvs_tmp[, "DPLX_NM"], collapse = ",")
      snvs_new2[new_row, "BULK_ASXS"] = paste(snvs_tmp[, "BULK_ASXS"], collapse = ",")
      snvs_new2[new_row, "BULK_NM"] = paste(snvs_tmp[, "BULK_NM"], collapse = ",")

      snvs_new2[new_row, "TYPE"] = snvs_tmp[1, "TYPE"]
    }
  }
  snvs_new2 = snvs_new2[order(snvs_new2$chr, snvs_new2$pos),]

  # drop some columns:
  snvs_new2 = snvs_new2[, c("chr", "pos", "kk", "ref", "mut", "qual", "filter", "TRI", "rb_id", "QPOS", "DEPTH_FWD",
                         "DEPTH_REV", "DEPTH_NORM_FWD", "DEPTH_NORM_REV", "TYPE", "TIMES_CALLED", "DPLX_ASXS", 
                         "DPLX_CLIP", "DPLX_NM", "BULK_ASXS", "BULK_NM")]


  ##########################################################################################
  # Calculate VAFs for subs (SNVs, DNVs, MNVs)
  # Duplex VAFs, use tabix to get the duplex coverage for each mutation:

  #cat("VAFs(1)...\n")
  # SUBS:
  snvs_new2_filt = snvs_new2[grep("PASS", snvs_new2$filter, invert = T),]
  snvs_new2_ok = snvs_new2[grep("PASS", snvs_new2$filter),]
  if (nrow(snvs_new2_filt) > 0) {
    snvs_new2_filt$DUPLEX_VAF = NA
    snvs_new2_filt$DUPLEX_COV = NA
  }
  if (nrow(snvs_new2_ok) > 0) {
    ranges = GRanges(snvs_new2_ok$chr, IRanges(start = snvs_new2_ok$pos, end = snvs_new2_ok$pos + nchar(snvs_new2_ok$ref)))
    res <- scanTabix(tbx, param = ranges)
    duplex_covs = as.numeric(unlist((sapply(res, function(x) unlist(strsplit(x[1], ";"))[3]))))
    snvs_new2_ok$DUPLEX_VAF = snvs_new2_ok$TIMES / duplex_covs
    snvs_new2_ok$DUPLEX_COV = duplex_covs
  }
  snvs_final = rbind(snvs_new2_ok, snvs_new2_filt)
  snvs_final = snvs_final[order(snvs_final$chr, snvs_final$pos),]

  #cat("VAFs(2)...\n")
  # BAM VAFs / bam2R
  for (i in c(1:nrow(snvs_final))) {
    kk = bam2R(dedup_bam, snvs_final[i, "chr"], snvs_final[i, "pos"], snvs_final[i, "pos"] + nchar(snvs_final[i, "ref"]) - 1, q = 30, mask = 3844, mq = 30)
    muts = unlist(strsplit(snvs_final[i, "mut"], ""))
    total_cov = 0
    total_mut = 0
    for (j in c(1:length(muts))) {
      total_mut = total_mut + kk[j, muts[j]] + kk[j, tolower(muts[j])]
      # total_cov = total_cov + sum(kk[j, c("A", "C", "G", "T", "a", "c", "g", "t","DEL","INS","del","ins")], na.rm = T)
      total_cov = total_cov + sum(kk[j, c("A", "C", "G", "T", "a", "c", "g", "t","DEL","del")], na.rm = T) # shouldn't count INS / ins
    }
    snvs_final[i, "BAM_MUT"] = total_mut
    snvs_final[i, "BAM_COV"] = total_cov
    snvs_final[i, "BAM_VAF"] = total_mut / total_cov
  }

  #cat("VAFs(3)...\n")
  # BAM VAFs / bam2R with relaxed BQ filters
  for (i in c(1:nrow(snvs_final))) {
    kk = bam2R(dedup_bam, snvs_final[i, "chr"], snvs_final[i, "pos"], snvs_final[i, "pos"] + nchar(snvs_final[i, "ref"]) - 1, q = 10, mask = 3844, mq = 30)
    muts = unlist(strsplit(snvs_final[i, "mut"], ""))
    total_cov = 0
    total_mut = 0
    for (j in c(1:length(muts))) {
      total_mut = total_mut + kk[j, muts[j]] + kk[j, tolower(muts[j])]
      # total_cov = total_cov + sum(kk[j, c("A", "C", "G", "T", "a", "c", "g", "t","DEL","INS","del","ins")], na.rm = T)
      total_cov = total_cov + sum(kk[j, c("A", "C", "G", "T", "a", "c", "g", "t","DEL","del")], na.rm = T) # shouldn't count INS / ins
    }
    snvs_final[i, "BAM_MUT_BQ10"] = total_mut
    snvs_final[i, "BAM_COV_BQ10"] = total_cov
    snvs_final[i, "BAM_VAF_BQ10"] = total_mut / total_cov
  }
}
# Finished with SNVs!
##########################################################################################
# Create INDEL ids and collapse repeated indel calls
#cat("Indels...\n")
if (num_indels > 0) {
  indels$info = gsub("=", ";", indels$info)
  for (i in c(1:nrow(indels))) {
    kk = unlist(strsplit(indels[i, "info"], ";"))
    for (j in seq(from = 2, to = length(kk), 2)) {
      col_name = kk[j]
      value = kk[j + 1]
      indels[i, col_name] = value
    }
  }
  indels$rb_id = gsub(",", ":", indels$RB)
  indels$TYPE = "del"
  indels[which(nchar(indels$ref) < nchar(indels$mut)), "TYPE"] = "ins"
  indels = indels[order(indels$chr, indels$pos, indels$rb_id),]

  # Collapse repeated indels:
  indels$mut_id = paste(indels$chr, indels$pos, indels$ref, indels$mut, sep = ":")

  counts_per_mut = table(indels$mut_id)
  repeat_muts = counts_per_mut[which(counts_per_mut > 1)]
  unique_muts = counts_per_mut[which(counts_per_mut == 1)]

  indels_new = indels[which(indels$mut_id %in% names(unique_muts)),]
  if (nrow(indels_new) > 0) {
    indels_new$TIMES_CALLED = 1
  }
  if (length(repeat_muts) > 0) {
    for (mut_id in names(repeat_muts)) {
      new_row = nrow(indels_new) + 1
      indels_tmp = indels[which(indels$mut_id == mut_id),]
      freq = nrow(indels_tmp)
      indels_new[new_row, "chr"] = indels_tmp[1, "chr"]
      indels_new[new_row, "pos"] = indels_tmp[1, "pos"]
      indels_new[new_row, "kk"] = indels_tmp[1, "kk"]
      indels_new[new_row, "ref"] = indels_tmp[1, "ref"]
      indels_new[new_row, "mut"] = indels_tmp[1, "mut"]
      indels_new[new_row, "qual"] = mean(indels_tmp[, "qual"])
      indels_new[new_row, "filter"] = indels_tmp[1, "filter"]
      indels_new[new_row, "rb_id"] = paste(indels_tmp[, "rb_id"], collapse = ",")
      indels_new[new_row, "TIMES_CALLED"] = freq
      indels_new[new_row, "TYPE"] = indels_tmp[1, "TYPE"]
      indels_new[new_row, "SEQ"] = indels_tmp[1, "SEQ"]
      indels_new[new_row, "MQ"] = indels_tmp[1, "MQ"]
      indels_new[new_row, "DP4"] = indels_tmp[1, "DP4"]
      indels_new[new_row, "mut_id"] = indels_tmp[1, "mut_id"]

      indels_new[new_row, "BBEG"] = paste(indels_tmp[, "BBEG"], collapse = ",")
      indels_new[new_row, "BEND"] = paste(indels_tmp[, "BEND"], collapse = ",")
      indels_new[new_row, "QPOS"] = paste(indels_tmp[, "QPOS"], collapse = ",")
      indels_new[new_row, "DEPTH_FWD"] = median(as.numeric(indels_tmp[, "DEPTH_FWD"]))
      indels_new[new_row, "DEPTH_REV"] = median(as.numeric(indels_tmp[, "DEPTH_REV"]))
      indels_new[new_row, "DEPTH_NORM_FWD"] = indels_tmp[1, "DEPTH_NORM_FWD"]
      indels_new[new_row, "DEPTH_NORM_REV"] = indels_tmp[1, "DEPTH_NORM_REV"]
      indels_new[new_row, "DPLX_ASXS"] = paste(indels_tmp[, "DPLX_ASXS"], collapse = ",")
      indels_new[new_row, "DPLX_CLIP"] = paste(indels_tmp[, "DPLX_CLIP"], collapse = ",")
      indels_new[new_row, "DPLX_NM"] = paste(indels_tmp[, "DPLX_NM"], collapse = ",")
      indels_new[new_row, "BULK_ASXS"] = paste(indels_tmp[, "BULK_ASXS"], collapse = ",")
      indels_new[new_row, "BULK_NM"] = paste(indels_tmp[, "BULK_NM"], collapse = ",")
    }
  }

  indels_new = indels_new[order(indels_new$chr, indels_new$pos),]

  # drop some columns:
  indels_new = indels_new[, c("chr", "pos", "kk", "ref", "mut", "qual", "filter", "rb_id", "TYPE", "TIMES_CALLED", "SEQ","BBEG","BEND","QPOS","DEPTH_FWD","DEPTH_REV","DEPTH_NORM_FWD","DEPTH_NORM_REV","DPLX_ASXS","DPLX_CLIP","DPLX_NM","BULK_ASXS","BULK_NM")]

  ##########################################################################################
  # Calculate VAFs for indels 
  # Duplex VAFs
  # INDELS:
  # Do it separately. Because indels don't pass filters, we need to sum them (TIMES_SEEN) to 
  # the duplex coverage seen in cov_bed
  #cat("VAFs(1)...\n")
  indels_new_filt = indels_new[grep("PASS", indels_new$filter, invert = T),]
  indels_new_ok = indels_new[grep("PASS", indels_new$filter),]
  if (nrow(indels_new_filt) > 0) {
    indels_new_filt$DUPLEX_VAF = NA
    indels_new_filt$DUPLEX_COV = NA
  }
  if (nrow(indels_new_ok) > 0) {
    ranges = GRanges(indels_new_ok$chr, IRanges(start = indels_new_ok$pos, end = indels_new_ok$pos + nchar(indels_new_ok$ref)))
    res <- scanTabix(tbx, param = ranges)
    duplex_covs = as.numeric(unlist((sapply(res, function(x) unlist(strsplit(x[1], ";"))[3]))))
    indels_new_ok$DUPLEX_VAF = indels_new_ok$TIMES / (duplex_covs + indels_new_ok$TIMES) # difference with subs
    indels_new_ok$DUPLEX_COV = duplex_covs + indels_new_ok$TIMES
  }
  indels_final = rbind(indels_new_ok, indels_new_filt)
  indels_final = indels_final[order(indels_final$chr, indels_final$pos),]

  #cat("VAFs(2)...\n")
  # BAM VAFs / bam2R
  for (i in c(1:nrow(indels_final))) {
    kk = bam2R(dedup_bam, indels_final[i, "chr"], indels_final[i, "pos"], indels_final[i, "pos"], q = 20, mask = 3844, mq = 30)
    # total_cov = sum(kk[1, c("A", "C", "G", "T", "a", "c", "g", "t","DEL","INS","del","ins")], na.rm = T)
    total_cov = total_cov + sum(kk[1, c("A", "C", "G", "T", "a", "c", "g", "t","DEL","del")], na.rm = T) # shouldn't count INS / ins
    if (indels_final[i, "TYPE"] == "del") {
      total_mut = sum(kk[, c("DEL", "del")])
    } else {
      #ins
      total_mut = sum(kk[, c("INS", "ins")])
    }
    indels_final[i, "BAM_MUT"] = total_mut
    indels_final[i, "BAM_COV"] = total_cov
    indels_final[i, "BAM_VAF"] = total_mut / total_cov
  }

  # BAM VAFs / bam2R with relaxed BQ filters
  #cat("VAFs(3)...\n")
  for (i in c(1:nrow(indels_final))) {
    kk = bam2R(dedup_bam, indels_final[i, "chr"], indels_final[i, "pos"], indels_final[i, "pos"], q = 10, mask = 3844, mq = 30)
    # total_cov = sum(kk[1, c("A", "C", "G", "T", "a", "c", "g", "t","DEL","INS","del","ins")], na.rm = T)
    total_cov = total_cov + sum(kk[1, c("A", "C", "G", "T", "a", "c", "g", "t","DEL","del")], na.rm = T) # shouldn't count INS / ins
    if (indels_final[i, "TYPE"] == "del") {
      total_mut = sum(kk[, c("DEL", "del")])
    } else {
      #ins
      total_mut = sum(kk[, c("INS", "ins")])
    }
    indels_final[i, "BAM_MUT_BQ10"] = total_mut
    indels_final[i, "BAM_COV_BQ10"] = total_cov
    indels_final[i, "BAM_VAF_BQ10"] = total_mut / total_cov
  }
}
#cat("Preparing final matrix SNVs...\n")
if (num_snvs > 0) {
  snvs_final$INFO = paste(rep("TRI=", nrow(snvs_final)), snvs_final$TRI, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("TIMES_CALLED=", nrow(snvs_final)), snvs_final$TIMES_CALLED, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("TYPE=", nrow(snvs_final)), snvs_final$TYPE, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("DUPLEX_VAF=", nrow(snvs_final)), snvs_final$DUPLEX_VAF, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("BAM_VAF=", nrow(snvs_final)), snvs_final$BAM_VAF, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("BAM_VAF_BQ10=", nrow(snvs_final)), snvs_final$BAM_VAF_BQ10, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("DEPTH_NORM_FWD=", nrow(snvs_final)), snvs_final$DEPTH_NORM_FWD, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("DEPTH_NORM_REV=", nrow(snvs_final)), snvs_final$DEPTH_NORM_REV, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("DEPTH_FWD=", nrow(snvs_final)), snvs_final$DEPTH_FWD, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("DEPTH_REV=", nrow(snvs_final)), snvs_final$DEPTH_REV, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("DUPLEX_COV=", nrow(snvs_final)), snvs_final$DUPLEX_COV, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("BAM_MUT=", nrow(snvs_final)), snvs_final$BAM_MUT, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("BAM_COV=", nrow(snvs_final)), snvs_final$BAM_COV, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("BAM_MUT_BQ10=", nrow(snvs_final)), snvs_final$BAM_MUT_BQ10, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("BAM_COV_BQ10=", nrow(snvs_final)), snvs_final$BAM_COV_BQ10, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("RB=", nrow(snvs_final)), snvs_final$rb_id, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("QPOS=", nrow(snvs_final)), snvs_final$QPOS, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("DPLX_ASXS=", nrow(snvs_final)), snvs_final$DPLX_ASXS, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("DPLX_CLIP=", nrow(snvs_final)), snvs_final$DPLX_CLIP, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("DPLX_NM=", nrow(snvs_final)), snvs_final$DPLX_NM, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("BULK_ASXS=", nrow(snvs_final)), snvs_final$BULK_ASXS, ";", sep = "")
  snvs_final$INFO = paste(snvs_final$INFO, rep("BULK_NM=", nrow(snvs_final)), snvs_final$BULK_NM, "", sep = "")

  snvs_final = snvs_final[, c("chr", "pos", "kk", "ref", "mut", "qual", "filter", "INFO")]
}

if (num_indels > 0) {
  #cat("Preparing final matrix Indels...\n")
  indels_final$INFO = paste(rep("TYPE=", nrow(indels_final)), indels_final$TYPE, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("TIMES_CALLED=", nrow(indels_final)), indels_final$TIMES_CALLED, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("DUPLEX_VAF=", nrow(indels_final)), indels_final$DUPLEX_VAF, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("BAM_VAF=", nrow(indels_final)), indels_final$BAM_VAF, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("BAM_VAF_BQ10=", nrow(indels_final)), indels_final$BAM_VAF_BQ10, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("DEPTH_NORM_FWD=", nrow(indels_final)), indels_final$DEPTH_NORM_FWD, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("DEPTH_NORM_REV=", nrow(indels_final)), indels_final$DEPTH_NORM_REV, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("DEPTH_FWD=", nrow(indels_final)), indels_final$DEPTH_FWD, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("DEPTH_REV=", nrow(indels_final)), indels_final$DEPTH_REV, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("SEQ=", nrow(indels_final)), indels_final$SEQ, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("DUPLEX_COV=", nrow(indels_final)), indels_final$DUPLEX_COV, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("BAM_MUT=", nrow(indels_final)), indels_final$BAM_MUT, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("BAM_COV=", nrow(indels_final)), indels_final$BAM_COV, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("BAM_MUT_BQ10=", nrow(indels_final)), indels_final$BAM_MUT_BQ10, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("BAM_COV_BQ10=", nrow(indels_final)), indels_final$BAM_COV_BQ10, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("RB=", nrow(indels_final)), indels_final$rb_id, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("QPOS=", nrow(indels_final)), indels_final$QPOS, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("DPLX_ASXS=", nrow(indels_final)), indels_final$DPLX_ASXS, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("DPLX_CLIP=", nrow(indels_final)), indels_final$DPLX_CLIP, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("DPLX_NM=", nrow(indels_final)), indels_final$DPLX_NM, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("BULK_ASXS=", nrow(indels_final)), indels_final$BULK_ASXS, ";", sep = "")
  indels_final$INFO = paste(indels_final$INFO, rep("BULK_NM=", nrow(indels_final)), indels_final$BULK_NM, "", sep = "")
}

# Describe in header:
if (length(grep("\\.gz", muts_vcf)) > 0) {
  ref_and_contigs = system(paste("zgrep '^##reference\\|^##contig' ", muts_vcf, "|| true"), intern = TRUE)
} else {
  ref_and_contigs = system(paste("grep '^##reference\\|^##contig' ", muts_vcf, "|| true"), intern = TRUE)
}
header = vector()
header[length(header) + 1] = "##fileformat=VCFv4.2"
header[length(header) + 1] = "##FILTER=<ID=PASS,Description=\"All filters passed\">"
header = c(header, ref_and_contigs)
header[length(header) + 1] = "##FILTER=<ID=NEI_IND,Description=\"Indel was found in an indel rich region of the matched normal\">"
header[length(header) + 1] = "##FILTER=<ID=MISSINGBULK,Description=\"Indel was not found in the matched normal\">"
header[length(header) + 1] = "##FILTER=<ID=MASKED,Description=\"Indel overlaps with SW or SNP site\">"
header[length(header) + 1] = "##FILTER=<ID=dbsnp,Description=\"SNV/MNV in common SNP site\">"
header[length(header) + 1] = "##FILTER=<ID=shearwater,Description=\"SNV/MNV in Noisy site\">"
header[length(header) + 1] = "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Type of mutation: snv, dnv, mnv, ins or del\">"
header[length(header) + 1] = "##INFO=<ID=TRI,Number=1,Type=String,Description=\"Pyrimidine trinucleotide context SNV\">"
header[length(header) + 1] = "##INFO=<ID=TIMES_CALLED,Number=1,Type=String,Description=\"Number of times mutation has been called\">"
header[length(header) + 1] = "##INFO=<ID=DUPLEX_VAF,Number=1,Type=Float,Description=\"Duplex VAF\">"
header[length(header) + 1] = "##INFO=<ID=DUPLEX_COV,Number=1,Type=Integer,Description=\"Duplex coverage at site\">"
header[length(header) + 1] = "##INFO=<ID=BAM_VAF,Number=1,Type=Float,Description=\"Deduplicated BAM VAF (bam2R: q=30,mask=3844,mq=30; q=20 for indels)\">"
header[length(header) + 1] = "##INFO=<ID=BAM_COV,Number=1,Type=Integer,Description=\"Deduplicated BAM coverage\">"
header[length(header) + 1] = "##INFO=<ID=BAM_VAF_BQ10,Number=1,Type=Float,Description=\"Deduplicated BAM VAF (bam2R: q=10,mask=3844,mq=30)\">"
header[length(header) + 1] = "##INFO=<ID=BAM_COV_BQ10,Number=1,Type=Integer,Description=\"Deduplicated BAM coverage\">"
header[length(header) + 1] = "##INFO=<ID=BAM_MUT,Number=1,Type=Integer,Description=\"Number of mut alleles in deduplicated BAM\">"
header[length(header) + 1] = "##INFO=<ID=BAM_MUT_BQ10,Number=1,Type=Integer,Description=\"Number of mut alleles in deduplicated BAM\">"
header[length(header) + 1] = "##INFO=<ID=QPOS,Number=.,Type=Integer,Description=\"Read position(s) closest to 5-prime end\">"
header[length(header) + 1] = "##INFO=<ID=DEPTH_FWD,Number=1,Type=Float,Description=\"Read bundle forward reads depth\">"
header[length(header) + 1] = "##INFO=<ID=DEPTH_REV,Number=1,Type=Float,Description=\"Read bundle reverse reads depth\">"
header[length(header) + 1] = "##INFO=<ID=DEPTH_NORM_FWD,Number=1,Type=Float,Description=\"Matched normal forward reads depth\">"
header[length(header) + 1] = "##INFO=<ID=DEPTH_NORM_REV,Number=1,Type=Float,Description=\"Matched normal reverse reads depth\">"
header[length(header) + 1] = "##INFO=<ID=RB,Number=.,Type=String,Description=\"Read bundle id(s): chr:breakpoints:barcodes\">"
header[length(header) + 1] = "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Sequence context for indels\">"
header[length(header) + 1] = "##INFO=<ID=DPLX_ASXS,Number=.,Type=Integer,Description=\"AS-XS for duplex\">"
header[length(header) + 1] = "##INFO=<ID=DPLX_CLIP,Number=.,Type=Integer,Description=\"Clipping for duplex\">"
header[length(header) + 1] = "##INFO=<ID=DPLX_NM,Number=.,Type=Integer,Description=\"Mismatches in duplex\">"
header[length(header) + 1] = "##INFO=<ID=BULK_ASXS,Number=.,Type=Integer,Description=\"AS-XS for bulk\">"
header[length(header) + 1] = "##INFO=<ID=BULK_NM,Number=.,Type=Integer,Description=\"Mismatches in bulk\">"
header[length(header) + 1] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

muts_final = snvs_final[, c("chr", "pos", "kk", "ref", "mut", "qual", "filter", "INFO")]
if (num_indels > 0) {
  muts_final = rbind(muts_final, indels_final[, c("chr", "pos", "kk", "ref", "mut", "qual", "filter", "INFO")])
}
write.table(header, file = output_file, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(muts_final, file = output_file, sep = "\t", quote = F, row.names = F, col.names = F, append = T)
