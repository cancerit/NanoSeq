# CHANGES

## 3.3.0

* Added new INFO fileds to the final vcf
* Changes to indel calling to reduce false calls
* Added burdens output to nanoseq plotter
* Fixed bugs in the VAF script

## 3.2.1

* Fixed bug in part with empty intervals at end of genome

## 3.2.0

* Modified BAM processing as to keep unmapped reads
* Fixed VCF filter bug in VAF script 

## 3.1.0

* Fixed type error of VAF script

## 3.0.0

* Made changes to accomodate targeted NanoSeq

## 2.3.3

* Added quotes add various places to support ant chromosomes
* Fixed error with 1 bp intervals
* Fixed efficiency calculation problem with coordinates outside of chr ends

## 2.3.2

* Reverted to old efficiency calculation

## 2.3.1

* Added safe creation of tmp dir
* Re-read input file to reduce mem usage in efficiency calculation

## 2.3.0

* Cleaned up/ simplified variantcaller.R
* Added extra check for turncated dsa output

## 2.2.0

* Tested with targeted data
* Fixed issue with coverage done files
* Fixed bamcov so that it handles chromosomes without coverage
* Fixed required/optional arguments from docustring
* Fixed incorrect default argument (-d)
* Added argument for file name of results
* Removed plugin option from htslib that was causing issues with CRAM
* Set seqs_per_slice to 1000 in htslib for improved CRAM streaming 
* Fixed off by one errors in part
* Fixed interval test in part
* Added consistency checks for contigs in BAM files 

## 2.1.0

* Added CRAM support
* Fixed issue with sort in indel pipeline

## 2.0.0

* Fixed 1-off error in dsa code.
* Incorporated Indel scripts.
* Incorporated new code for quick coverage estimation.
* Rewrote Python wrapper.
* Updated licenses.

## 1.5.3

* Fixed bug with last interval for outlier case

## 1.5.1

* Added more error captures for htslib calls

## 1.5.0

* Updated R library installation script
* Added option to skip post processing step
* Simplified excluded chromosome option
* Reduced default coverage window

## 1.4.0

* Initial public release
* Updated to htslib 1.11
* Added licencing information
* Added error checks when reading

