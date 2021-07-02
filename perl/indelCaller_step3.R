#!/usr/bin/env Rscript
#    Rscript to check whether each of the potential indel calls
#    happen at loci rich in indels in the matched normal. Those would be flagged as 'NEI_IND'.
#    Reliable indels are flagged with 'PASS'.
#    Three parameters are expected:
#      1. Reference genome
#      2. The name of the output file from step 2 of indelCaller
#      3. The BAM file containing the matched normal

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("indelCaller_step3.R  reference  vcf  bam\n\n Check identified indels against matched normal.\n\n Must specify: reference, VCF from indelCaller_step2, BAM file for the matched normal\n\n", call.=FALSE)
}

genomeFile = args[1]
vcf_file = args[2]
bam_file = args[3]

if (! file.exists(genomeFile)){
  stop("Reference file not found : ",genomeFile, call.=FALSE)
}
if (! file.exists(vcf_file) ){
  stop("VCF file not found : ",vcf_file, call.=FALSE)
}
if (! file.exists(bam_file) ){
  stop("Matched BAM file not found : ",bam_file, call.=FALSE)
}

library(deepSNV)
library(vcfR)
library("GenomicRanges")
library("Rsamtools")
library("MASS")

FLANK = 5

vcf <- read.vcfR( vcf_file, verbose = FALSE )

# Create regions:
#  * For deletions: get pos + length(deletion)   // length(deletion) = length(ref)-1
#  * For insertions: get pos
#  * Then, sum +/-5 to each side to count indels in the vecinity

for(i in c(1:nrow( getFIX(vcf) ))) {
	pos = strtoi( vcf@fix[i,"POS"] )
	chr = vcf@fix[i,"CHROM"]
	len = max(1,length( vcf@fix[i,"REF"])-1)
	start = pos - FLANK
	end   = pos + len + FLANK
	kk = bam2R(bam_file, chr,start,end, q=-100, mask=3844, mq=10) 
	# dont want a filter on BQ because in some bams BQ of indels have -1
    n_bases  = sum(kk[,c("A","C","G","T","a","c","g","t")])
	n_indels = sum(kk[,c("-","INS","_","ins")            ]) # Number of reads with an indel around the mutation
	cat(chr,pos,len,n_bases,n_indels,"\n");
	cat(n_bases,"/",n_indels,"\n")
  cat(i,"\n")
	if(n_bases == 0) {
    vcf@fix[i,"FILTER"] = "MISSINGBULK"
    vcf@fix[i,"INFO"] = paste(vcf@fix[i,"INFO"],";NN=[",n_indels,"/",n_bases,"]",sep="")

	} else if(n_indels/n_bases > 0.01) {
    vcf@fix[i,"FILTER"] = "NEI_IND"
		vcf@fix[i,"INFO"] = paste(vcf@fix[i,"INFO"],";NN=[",n_indels,"/",n_bases,"]",sep="")
	}	
	sequence = as.vector(scanFa(genomeFile, GRanges(chr, IRanges(start-3, end+3))))
  vcf@fix[i,"INFO"] = paste(vcf@fix[i,"INFO"],";SEQ=",sequence,sep="")

}

#Add new fields to the header
newhead = c()
once = TRUE
for( i in vcf@meta ) {
  if( once & grepl("##FILTER=<ID=MASKED", i, fixed=TRUE) ) {
    cat("here\n")
    i = paste("##FILTER=<ID=MISSINGBULK,Description=\"Site was not found in the matched normal\">\n",i,sep="")
    i = paste("##FILTER=<ID=NEI_IND,Description=\"Site was found in an indel rich region of the matched normal\">\n",i,sep="")
    i = paste("##INFO=<ID=NN,Number=1,Type=String,Description=\"n indels / n bases\">\n",i,sep="")
    i = paste("##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Sequence of indel plus flanking sequences\">\n",i,sep="")
    once = FALSE
  }
  newhead = c(newhead, i)
}

vcf@meta = newhead

out_vcf_file = gsub(".vcf",".filtered.vcf",vcf_file)
write.vcf( vcf, file=out_vcf_file)

