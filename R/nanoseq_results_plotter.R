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

suppressPackageStartupMessages({
  library(data.table)
  library(epitools)
  library(gridExtra)
  library(grid)
  library(Biostrings)
})

options(stringsAsFactors = FALSE)
READ_LENGTH = 151

args = commandArgs(TRUE)

if (length(args) == 0 || length(args) < 2 || length(args) > 3) {
  message("nanoseq_results_plotter.R  directory  output_prefix [trinucleotide_frequencies_file]\n\n")
  message("Must specify a directory with the various CSV files and a prefix for the output files.\n")
  message("Optionally a file containing the background genomic (pyrimidine-based) trinucleotide absolute counts (for normalization purposes). If not provided, human frequencies will be assumed.\n\n")
  message("This script generates a series of plots and tables summarizing the results of the variantcaller.\n\n")
  quit(save = "no", status = 1)
}

dir_res = args[1]
out_name = args[2]
if (length(args) == 3) {
  message(sprintf("Reading genomic trinucleotide absolute counts from file: %s\n", args[3]))
  order = c("ACA", "ACC", "ACG", "ACT", "ATA", "ATC", "ATG", "ATT", "CCA", "CCC", "CCG", "CCT", "CTA", "CTC", "CTG", "CTT",
            "GCA", "GCC", "GCG", "GCT", "GTA", "GTC", "GTG", "GTT", "TCA", "TCC", "TCG", "TCT", "TTA", "TTC", "TTG", "TTT")
  genome_counts_tmp = read.table(args[3], sep = "\t", stringsAsFactors = F, header = F, row.names = 1)
  if (nrow(genome_counts_tmp) != length(order)) {
    message(sprintf("%s not properly formatted. Expected %d rows, one for each pyrimidine trinucleotide\n", args[3], length(order)))
    quit(save = "no", status = 0)
  }
  if (sum(genome_counts_tmp[, 1]) < 1e6) {
    message(sprintf("%s does not seem to contain absolute counts.\n", args[3]))
    quit(save = "no", status = 0)
  }
  genome_counts = vector()
  genome_counts[rownames(genome_counts_tmp)] = genome_counts_tmp[, 1]
  genome_counts = genome_counts[order]
} else {
  # human genome trinucleotide frequencies assumed
  genome_counts = vector()
  genome_counts[c("ACA", "ACC", "ACG", "ACT", "ATA", "ATC", "ATG", "ATT", "CCA", "CCC", "CCG", "CCT", "CTA", "CTC", "CTG", "CTT",
                  "GCA", "GCC", "GCG", "GCT", "GTA", "GTC", "GTG", "GTT", "TCA", "TCC", "TCG", "TCT", "TTA", "TTC", "TTG", "TTT")] =
                  c(115415924, 66550070, 14381094, 92058521, 117976329, 76401029, 105094288, 142651503, 105547494, 75238490,
                    15801067, 101628641, 73791042, 96335416, 115950255, 114180747, 82414099, 68090507, 13621251, 80004082,
                    64915540, 54055728, 86012414, 83421918, 112085858, 88336615, 12630597, 126566213, 119020255, 112827451,
                    108406418, 219915599);
}

if (!dir.exists(dir_res)) {
  stop("Directory : ", dir_res, " not found", call. = FALSE)
}

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

dirname <- dir_res

# read in csv files
burdens <- fread(paste(dirname, 'burdens.csv', sep = "/"))
callsvsqpos <- fread(paste(dirname, 'callvsqpos.csv', sep = "/"))
coverage <- fread(paste(dirname, 'coverage.csv', sep = "/"))
pyrvsmask <- fread(paste(dirname, 'pyrvsmask.csv', sep = "/"))
readbundles <- fread(paste(dirname, 'readbundles.csv', sep = "/"))
variants <- fread(paste(dirname, 'variants.csv', sep = "/"))

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
unique_variants <- variants[!duplicated(variants[, c("chrom", "chromStart", "pyrsub")]),]

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

# summary metrics
#   unmasked
n_variants <- burdens[ismasked == 0][isvariant == 1]$count
n_reference <- burdens[ismasked == 0][isvariant == 0]$count
n_unique <- nrow(unique_variants[ismasked == 0])

if(length(n_reference) == 0) {
  message("0 reference calls. Exiting…\n")
  quit(save = "no", status = 1)
}

if (length(n_variants) == 0) { n_variants = 0; }
metrics <- data.frame("Metric" = c("total variants", "unique variants",
  "total variant + reference", "uncorrected burden", "physical coverage"),
  "Value" = c(n_variants, n_unique, (n_variants + n_reference),
  (n_variants / (n_variants + n_reference)), coverage))

print(metrics)


##########################################################################################
# First PDF: asymmetries, 6 subst types (12 really)
pdf(width = 5, height = 5, file = paste(out_name, ".subst_asym.pdf", sep = ""))
#colors <- rep(c("deepskyblue","black","firebrick2","gray","darkolivegreen3","rosybrown2"),each=2)
subs = c("C>A", "G>T", "C>G", "G>C", "C>T", "G>A", "T>A", "A>T", "T>C", "A>G", "T>G", "A>C")
colors = rep(c(rgb(34 / 255, 159 / 255, 198 / 255), rgb(26 / 255, 26 / 255, 26 / 255),
          rgb(201 / 255, 93 / 255, 94 / 255), rgb(178 / 255, 182 / 255, 180 / 255),
          rgb(153 / 255, 208 / 255, 62 / 255), rgb(217 / 255, 190 / 255, 217 / 255)), each = 2)
par(mar = c(5, 5, 5, 5))

# Masking:
unique_variants = unique_variants[which(unique_variants$ismasked == 0),]

unique_variants$from = substr(unique_variants$stdsub, 2, 2)
unique_variants$to = substr(unique_variants$stdsub, 5, 5)
unique_variants$sub = paste(unique_variants$from, ">", unique_variants$to, sep = "")
total = nrow(unique_variants)
counts = table(unique_variants$sub)
tmp_ = names(counts)
counts = as.vector(counts)
names(counts) = tmp_
counts[setdiff(subs, names(counts))] = 0
counts = counts[intersect(subs, names(counts))]
counts = counts[subs]
write.table(counts, file = paste(out_name, ".subst_asym.tsv", sep = ""), quote = F, sep = "\t", row.names = F, col.names = F)

tick_10 = 10 * total / 100
bar = barplot(counts, las = 2, col = colors,
              xlab = "Substitution (read strand)", ylab = "Number of mutations", border = "NA",
              ylim = c(0, max(counts) + 0.2 * max(counts)))
axis(side = 4, at = c(0, tick_10, tick_10 * 2, tick_10 * 3), labels = c("0", "10", "20", "30"))
mtext("% of mutations", side = 4, line = 2)
# Add significance indicators
pvalues = vector()
for (i in c(1, 3, 5, 7, 9, 11)) {
  j = i + 1
  if (counts[i] + counts[j] > 0) {
    pvalues[subs[i]] = binom.test(counts[i], counts[i] + counts[j])$p.value
  } else {
    pvalues[subs[i]] = 1
  }
}
write.table(pvalues, file = paste(out_name, ".subst_asym.pvals", sep = ""), quote = F, sep = "\t", row.names = T, col.names = F)
qvalues = p.adjust(pvalues, method = "BH")
for (i in c(1, 3, 5, 7, 9, 11)) {
  indicator = ""
  if (qvalues[subs[i]] < 0.01) {
    indicator = "**"
  } else if (qvalues[subs[i]] < 0.1) {
    indicator = "*"
  }
  text((bar[i] + bar[i + 1]) / 2, max(counts[i], counts[i + 1]), indicator, pos = 3)
}
dev.off()


##########################################################################################
# Second: the subst type asymmetries by bins
if (n_variants > 0) {
  message("       Calculating conf int for sub types...\n")
  pdf(width = 9, height = 7, file = paste(out_name, ".subst_asym_and_rates_binned.pdf", sep = ""))
  BIN_SIZE = 10


  # Now I need to do it by position:
  # callsvsqpos
  #    base qpos ismasked   count
  # 1:    A    0        0    2861
  # 2:    A    1        0 8376857
  # 3:    A    1        1   18505
  # 4:    A    2        0 5057418
  # Convert to:
  #    nucleotide qpos   count   BIN
  # 1           A   11 4824183 11-20
  # 2           A   12 4842390 11-20
  # 3           A   13 4831192 11-20
  # 4           A   14 4828102 11-20

  # Masking:
  nt_counts_ = callsvsqpos[which(callsvsqpos$ismasked == 0),]

  nt_counts_$qpos = nt_counts_$qpos + 1
  nt_counts_ = as.data.frame(nt_counts_)
  qposs = unique(callsvsqpos$qpos)
  nt_counts = data.frame()
  nt_counts[paste(c(rep("A", length(qposs)), rep("C", length(qposs)), rep("G", length(qposs)), rep("T", length(qposs))), qposs, sep = ":"), "qpos"] = as.numeric(rep(qposs, 4))
  nt_counts[paste(c(rep("A", length(qposs)), rep("C", length(qposs)), rep("G", length(qposs)), rep("T", length(qposs))), qposs, sep = ":"), "nucleotide"] = c(rep("A", length(qposs)), rep("C", length(qposs)), rep("G", length(qposs)), rep("T", length(qposs)))
  for (qpos in as.numeric(qposs)) {
    for (base in c("A", "C", "G", "T")) {
      row = paste(base, qpos, sep = ":")
      nt_counts[row, "count"] = sum(nt_counts_[which(nt_counts_$base == base & nt_counts_$qpos == qpos), "count"])
    }
  }
  nt_counts$BIN <- paste(floor(nt_counts$qpos / BIN_SIZE - 0.001) * BIN_SIZE + 1, "-", floor(nt_counts$qpos / BIN_SIZE - 0.001) * BIN_SIZE + BIN_SIZE, sep = "")
  mm = as.matrix(table(unique_variants[, c("sub", "qpos")]))
  nuc_specific = data.frame(matrix(nrow = 10, ncol = 3))
  colnames(nuc_specific) = c("sub", "qpos", "count")
  count = 1
  for (sub in rownames(mm)) {
    for (qpos in colnames(mm)) {
      nuc_specific[count, c("sub", "qpos", "count")] = c(sub, qpos, mm[sub, qpos])
      count = 1 + count
    }
  }
  nuc_specific$count = as.numeric(nuc_specific$count)
  nuc_specific$qpos = as.numeric(nuc_specific$qpos)
  nuc_specific = nuc_specific[which(nuc_specific$count > 0),]

  nt_muts <- nuc_specific
  nt_muts$qpos <- nt_muts$qpos + 1
  nt_muts$BIN <- paste(floor(nt_muts$qpos / BIN_SIZE - 0.001) * BIN_SIZE + 1, "-", floor(nt_muts$qpos / BIN_SIZE - 0.001) * BIN_SIZE + BIN_SIZE, sep = "")
  max_qpos = max(nt_counts$qpos, na.rm = T)
  bins <- paste(seq(from = 1, to = max_qpos, by = BIN_SIZE), "-", seq(from = BIN_SIZE, to = max_qpos + BIN_SIZE - 1, by = BIN_SIZE), sep = "")
  bins <- intersect(unique(nt_counts$BIN), bins)
  cosmic_colours <- c("deepskyblue", "black", "firebrick2", "gray", "darkolivegreen3", "rosybrown2")
  subs_pyr <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
  subs_pur <- c("G>T", "G>C", "G>A", "A>T", "A>G", "A>C");
  bin_subs_mat <- matrix(nrow = length(bins), ncol = 6 * 2 * 5)
  rownames(bin_subs_mat) <- bins
  colnames(bin_subs_mat) <- sort(c(paste(subs_pyr, rep(c("_counts", "_total", "_rate", "_rate_lci", "_rate_uci"), 6), sep = ""), paste(subs_pur, rep(c("_counts", "_total", "_rate", "_rate_lci", "_rate_uci"), 6), sep = "")))
  # Now fill the matrix:
  for (bin in bins) {
    for (sub in c(subs_pyr, subs_pur)) {
      nt <- unlist(strsplit(sub, ">"))[1]
      total <- paste(sub, "_total", sep = "")
      counts <- paste(sub, "_counts", sep = "")
      rate <- paste(sub, "_rate", sep = "")
      rate_lci <- paste(sub, "_rate_lci", sep = "")
      rate_uci <- paste(sub, "_rate_uci", sep = "")
      bin_subs_mat[bin, total] <- sum(nt_counts[which(nt_counts$BIN == bin & nt_counts$nucleotide == nt), "count"], na.rm = T)
      bin_subs_mat[bin, counts] <- sum(nt_muts[which(nt_muts$BIN == bin & nt_muts$sub == sub), "count"], na.rm = T) # somatic calls
      bin_subs_mat[bin, rate] <- bin_subs_mat[bin, counts] / bin_subs_mat[bin, total]
      if (bin_subs_mat[bin, total] != 0) {
        #conf.ints <- as.vector(as.matrix(binom.wilson(bin_subs_mat[bin,counts],bin_subs_mat[bin,total])[c("lower","upper")]))
        conf.ints <- poisson.test(bin_subs_mat[bin, counts])$conf.int / bin_subs_mat[bin, total]
        #conf.ints <- binom.test(bin_subs_mat[bin,counts],bin_subs_mat[bin,total])$conf.int
        bin_subs_mat[bin, rate_lci] <- conf.ints[1]
        bin_subs_mat[bin, rate_uci] <- conf.ints[2]
      } else {
        bin_subs_mat[bin, rate_lci] <- NA
        bin_subs_mat[bin, rate_uci] <- NA
      }
    }
  }

  bin_mat <- matrix(nrow = length(bins), ncol = 5)
  rownames(bin_mat) <- bins
  colnames(bin_mat) <- c("counts", "total", "rate", "rate_lci", "rate_uci")
  for (bin in bins) {
    bin_mat[bin, "counts"] <- sum(nt_muts[which(nt_muts$BIN == bin), "count"])
    bin_mat[bin, "total"] <- sum(nt_counts[which(nt_counts$BIN == bin), "count"])
    bin_mat[bin, "rate"] <- bin_mat[bin, "counts"] / bin_mat[bin, "total"]
    #conf.ints <- binom.test(bin_mat[bin,"counts"],bin_mat[bin,"total"])$conf.int
    conf.ints <- poisson.test(bin_mat[bin, "counts"])$conf.int / bin_mat[bin, "total"]
    #conf.ints <- as.vector(as.matrix(binom.wilson(bin_mat[bin,"counts"],bin_mat[bin,"total"])[c("lower","upper")]))
    bin_mat[bin, "rate_lci"] <- conf.ints[1]
    bin_mat[bin, "rate_uci"] <- conf.ints[2]
  }

  message("       Plotting results...\n")
  plotting_order <- c("C>A", "G>T", "C>G", "G>C", "C>T", "G>A", "T>A", "A>T", "T>C", "A>G", "T>G", "A>C")
  colors <- rep(c("deepskyblue", "black", "firebrick2", "gray", "darkolivegreen3", "rosybrown2"), each = 2)

  write.table(bin_subs_mat, file = paste(out_name, ".subst_asym_binned.tsv", sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)

  # IMPORTANT NOTE: mutation rates cannot be compared for the overall and by channel 
  # (by channel the territory is counted multiple times)

  par(mfrow = c(2, 3))
  #par(mar=c(3,2,2,2))
  par(mar = c(5, 3, 2, 1))
  tipos <- list()
  tipos[[1]] <- c("C>A", "G>T", "deepskyblue")
  tipos[[2]] <- c("C>G", "G>C", "black")
  tipos[[3]] <- c("C>T", "G>A", "firebrick2")
  tipos[[4]] <- c("T>A", "A>T", "gray")
  tipos[[5]] <- c("T>C", "A>G", "darkolivegreen3")
  tipos[[6]] <- c("T>G", "A>C", "rosybrown2")
  #The following three lines were just for testing something
  #bin_mat <- bin_mat[-which(rownames(bin_mat)=="1-10"),]
  #bin_subs_mat <- bin_subs_mat[-which(rownames(bin_subs_mat)=="1-10"),]
  #bins <- bins[c(2:length(bins))]
  for (j in c(1:length(tipos))) {
    tipo <- tipos[[j]]
    colorcete <- tipo[3]
    plot(bin_mat[, "rate"], type = "l", col = "gray", lwd = 1, xaxt = 'n', ylab = "Call rate", ylim = c(0, max(bin_subs_mat[, grep("_rate$", colnames(bin_subs_mat), perl = T)], na.rm = T)), xlab = "")
    segments(c(1:nrow(bin_mat)), bin_mat[, "rate_lci"], c(1:nrow(bin_mat)), bin_mat[, "rate_uci"], col = "gray")
    #plot(bin_mat[,"rate"],type="l",col="white",lwd=1,xaxt='n',ylab="Call rate",ylim=c(0,max(bin_subs_mat[,grep("_rate$",colnames(bin_subs_mat),perl=T)],na.rm=T)),xlab="")
    #segments(c(1:nrow(bin_mat)),bin_mat[,"rate_lci"],c(1:nrow(bin_mat)),bin_mat[,"rate_uci"],col="white")
    points(bin_subs_mat[, paste(tipo[1], "_rate", sep = "")], type = "l", col = colorcete, lty = 1, lwd = 2)
    points(bin_subs_mat[, paste(tipo[2], "_rate", sep = "")], type = "l", col = colorcete, lty = 2, lwd = 2)
    #segments(c(1:nrow(bin_mat)),bin_subs_mat[,paste(sub,"_rate_lci",sep="")],c(1:nrow(bin_mat)),bin_subs_mat[,paste(sub,"_rate_uci",sep="")],col="gray")
    axis(1, at = c(1:length(bins)), labels = bins, las = 2)
    legend(2, max(bin_subs_mat[, grep("_rate$", colnames(bin_subs_mat), perl = T)], na.rm = T), legend = c(tipo[1], tipo[2]), col = colorcete, lty = c(1, 2))
  }
  dev.off()
}

##########################################################################################
# Third PDF: trinuc profiles - four plots:
#   1. Observed counts
#   2. Diff to genome
#   3. Genome-corrected
#   4. In rates


# Masking:
pyrvsmask_unmasked = pyrvsmask[which(pyrvsmask$ismasked == 0),]

tri_bg = pyrvsmask_unmasked$count
names(tri_bg) = pyrvsmask_unmasked$pyrcontext
tri_bg[setdiff(names(genome_counts), names(tri_bg))] = 0

if (n_variants > 0) {

  colours = rep(c("deepskyblue", "black", "firebrick2", "gray", "darkolivegreen3", "rosybrown2"), each = 16)
  #colours        = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  sub_vec = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ctx_vec = paste(rep(c("A", "C", "G", "T"), each = 4), rep(c("A", "C", "G", "T"), times = 4), sep = "-")
  full_vec = paste(rep(sub_vec, each = 16), rep(ctx_vec, times = 6), sep = ",")
  xstr = paste(substr(full_vec, 5, 5), substr(full_vec, 1, 1), substr(full_vec, 7, 7), sep = "")
  ordered_names = paste(xstr, ">", rep(c("A", "G", "T", "A", "C", "G"), each = 16), sep = "")
  tmp_ = table(unique_variants[which(unique_variants$ismasked == 0), "pyrsub"])
  tri_obs = as.vector(tmp_)
  names(tri_obs) = names(tmp_)
  tri_obs[setdiff(ordered_names, names(tri_obs))] = 0
  tri_obs = tri_obs[ordered_names]

  pdf(width = 9, height = 13.5, file = paste(out_name, ".trinuc-profiles.pdf", sep = ""))
  par(mfrow = c(5, 1))
  par(mar = c(4, 6, 2, 2))
  # Observed counts
  y = tri_obs;
  maxy = max(y)
  h = barplot(y, las = 2, col = colours, border = NA, ylim = c(0, maxy * 1.5), space = 1, cex.names = 0.6, names.arg = xstr, ylab = "Observed mutation counts", main = "Observed mutation counts")
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j - 1) * 16 + 1, j * 16)]
    rect(xpos[1] - 0.5, maxy * 1.2, xpos[2] + 0.5, maxy * 1.3, border = NA, col = colours[j * 16])
    text(x = mean(xpos), y = maxy * 1.3, pos = 3, label = sub_vec[j])
  }

  # Difference to genome trinuc freqs
  difference_from_genome_wide <- as.vector((tri_bg / sum(as.numeric(tri_bg))) / (genome_counts / sum(genome_counts)))
  names(difference_from_genome_wide) = names(tri_bg)
  ratio2genome = difference_from_genome_wide
  write.table(cbind(tri_bg, ratio2genome), file = paste(out_name, ".trint_counts_and_ratio2genome.tsv", sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)
  barplot(difference_from_genome_wide, las = 2, main = "Obs/Genome trinuc. freqs", cex.names = .7, cex.lab = .7, col = rgb(0.2, 0.4, .6, .5), ylim = c(0.0, max(difference_from_genome_wide)), xpd = FALSE)
  abline(h = 1, col = "red", lty = 2, lwd = 2)

  tris = sapply(names(tri_obs), function(x) unlist(strsplit(x, ">"))[1])
  y = tri_obs / difference_from_genome_wide[tris];
  maxy = max(y)
  h = barplot(y, las = 2, col = colours, border = NA, ylim = c(0, maxy * 1.5), space = 1, cex.names = 0.6, names.arg = xstr, ylab = "Corrected mutation counts", main = "Corrected-to-genome mutation counts")
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j - 1) * 16 + 1, j * 16)]
    rect(xpos[1] - 0.5, maxy * 1.2, xpos[2] + 0.5, maxy * 1.3, border = NA, col = colours[j * 16])
    text(x = mean(xpos), y = maxy * 1.3, pos = 3, label = sub_vec[j])
  }
  trint_subst_obs = tri_obs
  trint_onto_genome = y
  trint_bg = tri_bg

  write.table(cbind(trint_subst_obs, trint_onto_genome), file = paste(out_name, ".trint_subs_obs_corrected.tsv", sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)


  barplot(tri_bg / sum(tri_bg), las = 2, main = "Obs trinuc. freqs", cex.names = .7, cex.lab = .7, col = rgb(0.2, 0.4, .6, .5), xpd = FALSE)

  tris = sapply(names(tri_obs), function(x) unlist(strsplit(x, ">"))[1])
  y = tri_obs / tri_bg[tris];
  maxy = max(y)
  h = barplot(y, las = 2, col = colours, border = NA, ylim = c(0, maxy * 1.5), space = 1, cex.names = 0.6, names.arg = xstr, ylab = "Mutation rates", main = "Mutation rates")
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j - 1) * 16 + 1, j * 16)]
    rect(xpos[1] - 0.5, maxy * 1.2, xpos[2] + 0.5, maxy * 1.3, border = NA, col = colours[j * 16])
    text(x = mean(xpos), y = maxy * 1.3, pos = 3, label = sub_vec[j])
  }
  dev.off()

  #############
  # Write variants:
  write.table(variants[which(variants$ismasked == 0),], file = paste(out_name, ".muts.tsv", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)

}


###########
# redo calculations including repeated mutations (important for clonal samples!)

# Masking:
pyrvsmask_unmasked = pyrvsmask[which(pyrvsmask$ismasked == 0),]

tri_bg = pyrvsmask_unmasked$count
names(tri_bg) = pyrvsmask_unmasked$pyrcontext
tri_bg[setdiff(names(genome_counts), names(tri_bg))] = 0

colours = rep(c("deepskyblue", "black", "firebrick2", "gray", "darkolivegreen3", "rosybrown2"), each = 16)
#colours        = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
sub_vec = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
ctx_vec = paste(rep(c("A", "C", "G", "T"), each = 4), rep(c("A", "C", "G", "T"), times = 4), sep = "-")
full_vec = paste(rep(sub_vec, each = 16), rep(ctx_vec, times = 6), sep = ",")
xstr = paste(substr(full_vec, 5, 5), substr(full_vec, 1, 1), substr(full_vec, 7, 7), sep = "")
ordered_names = paste(xstr, ">", rep(c("A", "G", "T", "A", "C", "G"), each = 16), sep = "")
tmp_ = table(variants[which(variants$ismasked == 0), "pyrsub"])
tri_obs = as.vector(tmp_)
names(tri_obs) = names(tmp_)
tri_obs[setdiff(ordered_names, names(tri_obs))] = 0
tri_obs = tri_obs[ordered_names]


# Difference to genome trinuc freqs
difference_from_genome_wide <- as.vector((tri_bg / sum(as.numeric(tri_bg))) / (genome_counts / sum(genome_counts)))
names(difference_from_genome_wide) = names(tri_bg)
ratio2genome = difference_from_genome_wide

tris = sapply(names(tri_obs), function(x) unlist(strsplit(x, ">"))[1])
y = tri_obs / difference_from_genome_wide[tris];
maxy = max(y)
trint_subst_obs = tri_obs
trint_onto_genome = y
trint_bg = tri_bg

###################
# Fourth: correct and save mutation burden
mut_burden = sum(trint_subst_obs) / sum(trint_bg)
mut_burden_corrected = sum(trint_onto_genome) / sum(trint_bg)
correction_factor = mut_burden_corrected / mut_burden

burdens = as.data.frame(matrix(nrow = 2, ncol = 5))
rownames(burdens) = c("observed", "corrected")
colnames(burdens) = c("muts", "total", "burden", "burden_lci", "burden_uci")

burdens["observed", "muts"] = sum(trint_subst_obs)
burdens["observed", "total"] = sum(trint_bg)
burdens["observed", "burden"] = sum(trint_subst_obs) / sum(trint_bg)
burdens["corrected", "muts"] = sum(trint_onto_genome)
burdens["corrected", "total"] = sum(trint_bg)
burdens["corrected", "burden"] = sum(trint_onto_genome) / sum(trint_bg)

# Add confidence intervals:
ci = poisson.test(sum(trint_subst_obs))$conf.int[1:2]
#ci = binom.wilson(sum(trint_subst_obs),sum(trint_bg))[,c("lower","upper")]
burdens["observed", "burden_lci"] = ci[1] / sum(trint_bg)
burdens["observed", "burden_uci"] = ci[2] / sum(trint_bg)
ci_corrected = ci * correction_factor
#ci = binom.wilson(burdens["corrected","muts"],burdens["corrected","total"])[,c("lower","upper")]
burdens["corrected", "burden_lci"] = ci_corrected[1] / sum(trint_bg)
burdens["corrected", "burden_uci"] = ci_corrected[2] / sum(trint_bg)

# Save to table:
write.table(burdens, file = paste(out_name, ".mut_burden.tsv", sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)


##########################################################################################
# Fifth: save mutations in VCF format 
# now done directly in the python wrapper


##########################################################################################
# Plot burden w/o masking common SNPs
# Qualitative check for contamination
burdens <- fread(paste(dirname, 'burdens.csv', sep = "/"))
burdens <- burdens[, .(count = sum(count)), by = .(ismasked, isvariant)]

burdens = as.data.frame(burdens)

##################
# fa8 / ao7 bug fix
if (nrow(burdens[which(burdens$ismasked == 0 & burdens$isvariant == 1), ]) > 0) {
  n_variants_masked <- burdens[which(burdens$ismasked == 0 & burdens$isvariant == 1), "count"]
} else {
  n_variants_masked = 0
}
n_reference_masked <- burdens[which(burdens$ismasked == 0 & burdens$isvariant == 0), "count"]
if (nrow(burdens[which(burdens$isvariant == 1), ]) > 0) {
  n_variants_unmasked <- sum(burdens[which(burdens$isvariant == 1), "count"])
} else {
  n_variants_unmasked = 0
}
n_reference_unmasked <- sum(burdens[which(burdens$isvariant == 0), "count"])
##################
burden_masked = n_variants_masked / (n_reference_masked + n_variants_masked)
burden_unmasked = n_variants_unmasked / (n_reference_unmasked + n_variants_unmasked)
ci_masked = poisson.test(n_variants_masked)$conf.int / (n_reference_masked + n_variants_masked)
ci_unmasked = poisson.test(n_variants_unmasked)$conf.int / (n_reference_unmasked + n_variants_unmasked)
if (n_variants_unmasked > 0 & n_variants_masked > 0 ) {
  # fa8
  pdf(width = 5, height = 5, file = paste(out_name, ".burden.masked-vs-unmasked.pdf", sep = ""))
  bar = barplot(c(burden_masked, burden_unmasked), names = c("burden (masked)", "burden (unmasked)"), ylim = c(0, max(burden_masked, burden_unmasked) + 0.1 * max(burden_masked, burden_unmasked)), main = "Qualitative contamination check")
  segments(bar[1], ci_masked[1], bar[1], ci_masked[2], col = "black", lwd = 2)
  segments(bar[2], ci_unmasked[1], bar[2], ci_unmasked[2], col = "black", lwd = 2)
  dev.off()
}
# fa8: save burdens pre / post masking in numeric format
tosave = as.data.frame(matrix(nrow=2,ncol=5))
colnames(tosave) = c("subs","total","burden","burden_lci","burden_uci")
rownames(tosave) = c("masked","unmasked")
tosave["masked",] = c(n_variants_masked, n_reference_masked+ n_variants_masked,
                      n_variants_masked / (n_reference_masked + n_variants_masked),
                      ci_masked[1],
                      ci_masked[2])
tosave["unmasked",] = c(n_variants_unmasked, n_reference_unmasked + n_variants_unmasked,
                      n_variants_unmasked / (n_reference_unmasked + n_variants_unmasked),
                      ci_unmasked[1],
                      ci_unmasked[2])
write.table(tosave,file=paste(out_name, ".obs_burdens.pre_vs_post_masking.tsv", sep = ""),sep="\t",row.names=T,col.names=T,quote=F)
# end.

##########################################################################################
# Mismatches
# To understand error profile - damage during library preparation
mismatches = read.table(paste(dirname, 'mismatches.csv', sep = "/"), sep = ",", header = T)

if(nrow(mismatches) > 0) {
  mismatches$sub = paste(unlist(sapply(mismatches$mismatch, function(x) unlist(strsplit(x, ""))[2])), ">", unlist(sapply(mismatches$mismatch, function(x) unlist(strsplit(x, ""))[5])), sep = "")
  mismatches = mismatches[which(mismatches$ismasked == 0),]
  mismatches = mismatches[which(mismatches$mismatch != "UNKNOWN"),]
  total = nrow(mismatches)

  pdf(width = 5, height = 5, file = paste(out_name, ".mismatches.subst_asym.pdf", sep = ""))

  subs = c("C>A", "G>T", "C>G", "G>C", "C>T", "G>A", "T>A", "A>T", "T>C", "A>G", "T>G", "A>C")
  colors = rep(c(rgb(34 / 255, 159 / 255, 198 / 255), rgb(26 / 255, 26 / 255, 26 / 255),
            rgb(201 / 255, 93 / 255, 94 / 255), rgb(178 / 255, 182 / 255, 180 / 255),
            rgb(153 / 255, 208 / 255, 62 / 255), rgb(217 / 255, 190 / 255, 217 / 255)), each = 2)
  par(mar = c(5, 5, 5, 5))
  counts = table(mismatches$sub)
  tmp_ = names(counts)
  counts = as.vector(counts)
  names(counts) = tmp_
  counts[setdiff(subs, names(counts))] = 0
  counts = counts[intersect(subs, names(counts))]
  counts = counts[subs]
  write.table(counts, file = paste(out_name, ".mismatches.subst_asym.tsv", sep = ""), quote = F, sep = "\t", row.names = F, col.names = F)

  tick_10 = 10 * total / 100
  bar = barplot(counts, las = 2, col = colors,
                xlab = "Substitution (read strand)", ylab = "Number of mutations", border = "NA",
                ylim = c(0, max(counts) + 0.2 * max(counts)))
  axis(side = 4, at = c(0, tick_10, tick_10 * 2, tick_10 * 3), labels = c("0", "10", "20", "30"))
  mtext("% of mutations", side = 4, line = 2)
  # Add significance indicators
  pvalues = vector()
  for (i in c(1, 3, 5, 7, 9, 11)) {
    j = i + 1
    if (counts[i] + counts[j] > 0) {
      pvalues[subs[i]] = binom.test(counts[i], counts[i] + counts[j])$p.value
    } else {
      pvalues[subs[i]] = 1
    }
  }
  write.table(pvalues, file = paste(out_name, ".mismatches.subst_asym.pvals", sep = ""), quote = F, sep = "\t", row.names = T, col.names = F)
  qvalues = p.adjust(pvalues, method = "BH")
  for (i in c(1, 3, 5, 7, 9, 11)) {
    indicator = ""
    if (qvalues[subs[i]] < 0.01) {
      indicator = "**"
    } else if (qvalues[subs[i]] < 0.1) {
      indicator = "*"
    }
    text((bar[i] + bar[i + 1]) / 2, max(counts[i], counts[i + 1]), indicator, pos = 3)
  }
  dev.off()

  #########
  # Now the trinuc profiles
  pdf(width = 15, height = 8, file = paste(out_name, ".mismatches.trinuc-profile.pdf", sep = ""))
  par(mfrow = c(2, 1))
  par(mar = c(4, 6, 2, 2))

  complement_subs2pyr = function(sub) {
    first = unlist(strsplit(sub, ""))[1]
    from = unlist(strsplit(sub, ""))[2]
    third = unlist(strsplit(sub, ""))[3]
    to = unlist(strsplit(sub, ""))[5]
    complement = vector()
    complement[c("A", "C", "G", "T")] = c("T", "G", "C", "A")
    pyr = c("C", "T")
    if (from %in% pyr) {
      return(sub)
    } else {
      sub = paste(complement[third], complement[from], complement[first], ">", complement[to], sep = "")
      return(sub)
    }
  }

  mismatches$mismatch_pyr = sapply(mismatches$mismatch, function(x) complement_subs2pyr(x))

  colours = rep(c("deepskyblue", "black", "firebrick2", "gray", "darkolivegreen3", "rosybrown2"), each = 16)
  sub_vec = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ctx_vec = paste(rep(c("A", "C", "G", "T"), each = 4), rep(c("A", "C", "G", "T"), times = 4), sep = "-")
  full_vec = paste(rep(sub_vec, each = 16), rep(ctx_vec, times = 6), sep = ",")
  xstr = paste(substr(full_vec, 5, 5), substr(full_vec, 1, 1), substr(full_vec, 7, 7), sep = "")
  ordered_names = paste(xstr, ">", rep(c("A", "G", "T", "A", "C", "G"), each = 16), sep = "")
  tmp_ = table(mismatches[which(mismatches$mismatch == mismatches$mismatch_pyr), "mismatch"])
  tri_obs = as.vector(tmp_)
  names(tri_obs) = names(tmp_)
  tri_obs[setdiff(ordered_names, names(tri_obs))] = 0
  tri_obs = tri_obs[ordered_names]

  # Observed counts
  y = tri_obs;
  maxy = max(y)
  h = barplot(y, las = 2, col = colours, border = NA, ylim = c(0, maxy * 1.5), space = 1, cex.names = 0.6, names.arg = xstr, ylab = "Observed Pyrimidine mismatches counts", main = "Observed Pyirimidine mismatches counts")
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j - 1) * 16 + 1, j * 16)]
    rect(xpos[1] - 0.5, maxy * 1.2, xpos[2] + 0.5, maxy * 1.3, border = NA, col = colours[j * 16])
    text(x = mean(xpos), y = maxy * 1.3, pos = 3, label = sub_vec[j])
  }
  write.table(tri_obs, file = paste(out_name, ".SSC-mismatches-Pyrimidine.triprofiles.tsv", sep = ""), quote = F, sep = "\t", row.names = T, col.names = F)
  tri_obs_pyr = tri_obs

  tmp_ = table(mismatches[which(mismatches$mismatch != mismatches$mismatch_pyr), "mismatch_pyr"])
  tri_obs = as.vector(tmp_)
  names(tri_obs) = names(tmp_)
  tri_obs[setdiff(ordered_names, names(tri_obs))] = 0
  tri_obs = tri_obs[ordered_names]

  # Observed counts
  y = tri_obs;
  maxy = max(y)
  h = barplot(y, las = 2, col = colours, border = NA, ylim = c(0, maxy * 1.5), space = 1, cex.names = 0.6, names.arg = xstr, ylab = "Observed Purine mismatches counts", main = "Observed Purine mismatches counts")
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j - 1) * 16 + 1, j * 16)]
    rect(xpos[1] - 0.5, maxy * 1.2, xpos[2] + 0.5, maxy * 1.3, border = NA, col = colours[j * 16])
    text(x = mean(xpos), y = maxy * 1.3, pos = 3, label = sub_vec[j])
  }
  dev.off()
  write.table(tri_obs, file = paste(out_name, ".SSC-mismatches-Purine.triprofiles.tsv", sep = ""), quote = F, sep = "\t", row.names = T, col.names = F)
  tri_obs = tri_obs_pyr + tri_obs
  write.table(tri_obs, file = paste(out_name, ".SSC-mismatches-Both.triprofiles.tsv", sep = ""), quote = F, sep = "\t", row.names = T, col.names = F)


  ##########################################################################################
  # Now error rates based on SSC, by each of 96 trinuc-sub channels and global
  tmp_ = table(mismatches[which(mismatches$mismatch == mismatches$mismatch_pyr), "mismatch"])
  tri_obs = as.vector(tmp_)
  names(tri_obs) = names(tmp_)
  tri_obs[setdiff(ordered_names, names(tri_obs))] = 0
  tri_obs = tri_obs[ordered_names]
  tri_obs_pyr = tri_obs
  tmp_ = table(mismatches[which(mismatches$mismatch != mismatches$mismatch_pyr), "mismatch_pyr"])
  tri_obs = as.vector(tmp_)
  names(tri_obs) = names(tmp_)
  tri_obs[setdiff(ordered_names, names(tri_obs))] = 0
  tri_obs = tri_obs[ordered_names]
  tri_obs_pur = tri_obs

  pdf(width = 15, height = 4, file = paste(out_name, ".DSC_estimated_error_rates.pdf", sep = ""))
  par(mar = c(4, 6, 2, 2))
  tris = sapply(names(tri_obs), function(x) unlist(strsplit(x, ">"))[1])
  y = (tri_obs_pyr / (tri_bg[tris] / 2)) * (tri_obs_pur / (tri_bg[tris] / 2));
  maxy = max(y)
  ##y = (tri_obs_pyr/(tri_bg[tris])) * (tri_obs_pur/(tri_bg[tris])); maxy = max(y)
  h = barplot(y, las = 2, col = colours, border = NA, ylim = c(0, maxy * 1.5), space = 1, cex.names = 0.6, names.arg = xstr, ylab = "Error rates", main = "Error rates")
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j - 1) * 16 + 1, j * 16)]
    rect(xpos[1] - 0.5, maxy * 1.2, xpos[2] + 0.5, maxy * 1.3, border = NA, col = colours[j * 16])
    text(x = mean(xpos), y = maxy * 1.3, pos = 3, label = sub_vec[j])
  }
  dev.off()


  # Global:
  # Number of errors per channel:
  errors_per_channel = y * tri_bg[tris]
  total_error_rate = sum(errors_per_channel) / sum(tri_bg)

  pdf(width = 15, height = 4, file = paste(out_name, ".DSC_errors_per_channel.pdf", sep = ""))
  par(mar = c(4, 6, 2, 2))
  tris = sapply(names(tri_obs), function(x) unlist(strsplit(x, ">"))[1])
  y = errors_per_channel;
  maxy = max(y)
  h = barplot(y, las = 2, col = colours, border = NA, ylim = c(0, maxy * 1.5), space = 1, cex.names = 0.6, names.arg = xstr, ylab = "Predicted errors", main = "Predicted errors")
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j - 1) * 16 + 1, j * 16)]
    rect(xpos[1] - 0.5, maxy * 1.2, xpos[2] + 0.5, maxy * 1.3, border = NA, col = colours[j * 16])
    text(x = mean(xpos), y = maxy * 1.3, pos = 3, label = sub_vec[j])
  }
  dev.off()

  output = list()
  output[["errors_per_channel"]] = errors_per_channel
  output[["error_rates_per_channel"]] = y
  output[["total_error_rate"]] = total_error_rate
  output[["total_errors"]] = sum(errors_per_channel)
  kkk = vector()
  kkk["total_error_rate"] = total_error_rate
  kkk["total_errors"] = output[["total_errors"]]
  write.table(kkk, file = paste(out_name, ".estimated_error_rates.tsv", sep = ""), quote = F, sep = "\t", row.names = T, col.names = F)
}
