#!/usr/bin/env Rscript

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

#install R packages
#install R packages
args = commandArgs(T)
instLib = args[1] 
r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, ask=FALSE, lib=instLib, lib.loc=instLib)
  sapply(pkg, library, character.only = TRUE)
}

ipak_bioc <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg, ask=FALSE, lib=instLib, lib.loc=instLib)
  sapply(pkg, library, character.only = TRUE)
}

if( (version$major == 3 && version$minor >=5) || version$major > 3) {
  # biocmanager versions of R
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(ask=FALSE, lib=instLib, lib.loc=instLib)
  ipak_bioc(c("deepSNV"))
  ipak_bioc(c("vcfR"))
  ipak_bioc(c("VGAM"))
  #ipak_bioc(c("data.table"))
  #ipak_bioc(c("epitools")) 
  #ipak_bioc(c("gridExtra"))
  #ipak_bioc(c("seqinr"))
} else {
  # OLD versions of R
  stop("Must update to R v3.5.0 or greater")
}
