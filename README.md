# NanoSeq

Nanorate sequencing (NanoSeq) is a DNA library preparation and sequencing protocol based on Duplex Sequencing ([Schmitt et al, 2012](https://doi.org/10.1073/pnas.1208715109)) and BotSeqS ([Hoang et al, 2016](https://doi.org/10.1073/pnas.1607794113)). NanoSeq allows calling mutations with single molecule resolution and extremely low error rates ([Abascal et al, 2021](https://doi.org/10.1038/s41586-021-03477-4)). The pipeline and code in this repository cover the preprocessing of NanoSeq sequencing data, the assessment of data quality and efficiency, and the calling of mutations (substitutions and indels) and the estimating mutation burdens and substitution profiles.

The wet-lab protocol is described in the original publication ([Abascal et al, 2021](https://doi.org/10.1038/s41586-021-03477-4)) and on ProtocolExchange ([Lensing et al, 2021](https://protocolexchange.researchsquare.com/article/pex-1298/v1)).

## Dependencies

Execution of the scripts from this repository requires that these dependencies are on PATH :

* samtools
* bcftools
* Rscript
* biobambam
* bwa

## Installation

```
./setup.sh path_to_install                          #install code from this repo
export PATH=$PATH:path_to_install/bin
Rscript ./build/manualInstall.R <R libraries path>  #install all the required R libraries
```

## Preprocessing of the sequencing data

### Steps overview

1) Extract the duplex barcodes from the fastq files and add them to the fastq header of each read (`extract_tags.py`)
2) Map reads to the reference genome using bwa with option -C to add the barcodes as tags in the bam
3) Add rc and mc tags, mark optical duplicates, and filter the bam for unpaired reads, creating a molecule-unique read bundle (RB) tag identifier for each read pair.

### 1/2. Extract barcodes and map reads

Prior to mapping, fastq files must be pre-processed with `extract_tags.py` in order to trim adapter sequences and to add the appropiate tags (rb,mb) to the read headers.
```
#(trim 3 bases, skip 4 bases, add rb & mb tags), for reads of read length 151 bps
python extract-tags.py -a R1.fastq -b R2.fastq -c extrR1.fastq -d extrR2.fastq -m 3 -s 4 -l 151

#(align with bwa appending rb & mb tags)
bwa mem -C reference_genome.fa extrR1.fastq extrR2.fastq > mapped.sam
```

### 3. Add rc and mc tags, mark optical duplicates, create read bundle tags

A read bundle tag must be appended to each read-pair of a BAM to determine which reads are PCR duplicates. The tag consists of: chromosome, read coordinate, mate corrdinate, read rb tag, and mate mb tag (RB:rc,mc,rb,mb).

With `bamsormadup` the rc and mc tags are added. We also recommend running `bammarkduplicatesopt` with optminpixeldif=2500 to flag optical duplicates.

```
#(sort sam with biobambam, append rc & mc tags)
bamsormadup inputformat=sam rcsupport=1 threads=1 < mapped.sam > mapped_od.bam
```

With `bamaddreadbundles`, optical duplicates and unpaired mates are filtered and the RB tag is created:
```
#(append RB tag, filter OD and unpaired read mates)
bamaddreadbundles -I mapped_od.bam -O filtered.bam
```

## Preparation of a matched normal

The NanoSeq analysis requires a matched normal to distinguish somatic mutations from germline SNPs. Sequencing undiluted NanoSeq libraries is the most cost-efficient solution to create a matched normal because all the coverage will concentrate in the fraction of the genome "seen" with the selected restriction enzyme. If the matched normal happens to be an undiluted NanoSeq library it must be processed further as to just keep one read-pair from each read bundle to produce a 'neat' normal (i.e. to remove PCR duplicates).

```
randomreadinbundle -I filtered.bam -O neat.bam
```

Deduplication of bams with `randomreadinbundle` is also required to run VerifyBamId and the efficiency estimates (`efficiency_nanoseq.pl`).

### Note

Correct pre-processing means that duplex BAMs must have @PG tags for `bamsormadup` , `bammarkduplicatesopt` and `bamaddreadbundles`. A neat bulk (NanoSeq library) BAM must have a @PG tag for `bamaddreadbundles` & `randomreadinbundle`. A WGS bulk will NOT have a tag for `bamaddreadbundles`. 

The pipeline checks that all these programs have been run on the bam, exiting with an error otherwise. At users' own risk, the bam header checking can be disabled using option `--no_test` in `dsa`.

## Contamination check

It is highly recommended to carry out a contamination check of the sample pair with [`verifyBAMId`](https://github.com/Griffan/VerifyBamID). This contamination check must be done on a bam generated with `randomreadinbundle`, where only one read per read bundle is kept in the bam (see above).
An alpha < 0.005 would be acceptable for most situations.

## Efficiency

The script efficiency_nanoseq.pl analyses the information in the NanoSeq original bam and its deduplicated version.
The output provides information on duplicate rates, read counts... Theoretically, the optimal duplicate rate in terms of efficiency (duplex bases / sequenced bases) is 81% for read bundles of size >= 2+2, with 65% and 90% yielding ≥80% of the maximum of efficiency. Empirically the optimal duplicate rate is 75-76%.

Apart of the duplicate rate, the following outputs are important to assess the quality of the experiment: F-EFF, EFFICIENCY, GC_BOTH/GC_SINGLE

F-EFF or strand drop out fraction: This shows the fraction of read bundles missing one of the two original strands beyond what would be expected under random sampling (assuming a binomial process). Good values are between 0.10-0.30, and larger values are likely due to DNA damage such as modified bases or internal nicks that prevent amplification of one of the two strands. Larger values do not impact the quality of the results, just reduce the efficiency of the protocol.

EFFICIENCY: This is the number of duplex bases divided by the number of sequenced bases. Efficiency is maximised to ~0.07 when duplicate rates and strand drop outs are optimal

GC_BOTH and GC_SINGLE: the GC content of RBs with both strands and with just one strand. The two values should be similar between them and similar to the genome average. If there are large deviations that is possibly due to biases during PCR amplification. If GC_BOTH is substantially larger than GC_SINGLE, DNA denaturation before dilution may have taken place.


## NanoSeq analysis

For a matched normal and a duplex pair of samples (hereafter referred to as "normal" and "tumour") an analysis requires the following steps:

1) Generation of tables files (`dsa`)
2) Generation of variant files (`variants`)
3) Indel identification (`indelCaller_step1.pl`, `indelCaller_step2.pl` & `indelCaller_step3.R`)
4) Summarizing of results (`variantcaller.R` & `nanoseq_results_plotter.R`)

The wrapper script `runNanoSeq.py` provides a convinient way of running all the steps of a NanoSeq analyisis. It is meant to be run as a job array or in a multithreaded environment.

The wrapper scrpt has subcommands that that are meant to roughly follow the same steps that were outlined before. The steps are: cov, part, dsa, var, indel & post. Output for each step is located in the tmpNanoSeq directory.

### Coverage (cov)

A coverage histogram is computed for the bulk BAM.

```
runNanoSeq.py -t 10 \
  -A normal.bam \
  -B tumour.bam \
  -R genome.fa \
  cov \
  -Q 0 \
  --exclude "MT,GL%,NC_%,hs37d5"
```

### Partition (part)

Divide the coverage so that each job in the NanoSeq analysis gets roughly the same work. The -n argument idicates the number of tasks that will be used in the dsa, var and indel steps.

```
runNanoSeq.py -t 1 \
  -A normal.bam \
  -B tumour.bam \
  -R genome.fa \
  part \
  -n 60 \
```

### dsa beds (dsa)

Compute the dsa bed files. SNP and NOISE BED files contain sites to be marked on the output VCF file.

```
runNanoSeq.py -t 60 \
  -A normal.bam \
  -B tumour.bam \
  -R genome.fa \
  dsa \
  -C SNP.sorted.bed.gz \
  -D NOISE.sorted.bed.gz \
  -d 2 \
  -q 30 \
```

### Variant tables (var)

Compute variants tables.

```
runNanoSeq.py -t 60 \
  -A normal.bam \
  -B tumour.bam \
  -R genome.fa \
  var \
  -a 50 \
  -b 5 \
  -c 0 \
  -f 0.9 \
  -i 1 \
  -m 8 \
  -n 3 \
  -p 0 \
  -q 60 \
  -r 144 \
  -v 0.01 \
  -x 8 \
  -z 12
```

### Indel vcfs (indel)

Compute vcf files for indels.

```
runNanoSeq.py -t 60 \
  -A normal.bam \
  -B tumour.bam \
  -R genome.fa \
  indel \
  -s sample \
  --rb 2 \
  --t3 135 \
  --t5 10 \
  --mc 16
```
### Post processing (post)

Merge final files, produce summaries, compute efficiency. Results can be found in tmpNanoSeq/post.

```
runNanoSeq.py -t 2 \
  -A normal.bam \
  -B tumour.bam \
  -R genome.fa \
  post
```
### Note

A bash script is provided in the LSF directory as a template for execution of all the steps in an computing environment using the LSF job scheduler. Script should be edited to fit individual needs.

### Genomic masks

Genomic masks for common SNP masking and detection of noisy/variable genomic sites. Masks for GRCh37 are available [here](https://drive.google.com/drive/folders/1wqkgpRTuf4EUhqCGSLA4fIg9qEEw3ZcL?usp=sharing). 

When running NanoSeq on a species for which no common SNP data is available, an empty bed file should be used.


## Output

The most relevant summary files include the following.

* `muts.vcf.gz / muts.tsv`: substitutions called in vcf and tsv format. "PASS" substitutions are those not filtered by the common SNP and noisy sites masks (see Genomic masks).
* `indels.vcf.gz`: indel calls
* `burden.masked-vs-unmasked.pdf`: estimated burden before and after filtering common SNPs. Provides a qualitative view on contamination.
* `mut_burden.tsv`: estimated substitution burden with Poisson confidence intervals. The corrected burden shows the burden after normalizing observed trinucleotide frequencies to the genomic trinucleotide frequencies.
* `trinuc-profiles.pdf / trint_subs_obs_corrected.tsv / trint_counts_and_ratio2genome.tsv`: Trinucleotide substitution profiles (observed and corrected), using the trinucleotide substitution counts in `trint_subs_obs_corrected.tsv` and the normalization of trinucleotide frequencies in `trint_counts_and_ratio2genome.tsv`. Normalization is required because NanoSeq results are depleted of trinucleotides overlapping the restriction site and of CpGs due to extensive filtering of common SNPs.
* `cov.bed.gz`: large file containing the effective duplex coverage for each genomic site, also showing the trinucleotide context of each site. This file is required to to calculate burdens and substitution profiles in sets of specific genomic regions (e.g. highly expressed genes, heterochromatin, ...).
* `subst_asym.pdf / subst_asym_and_rates_binned.pdf / subst_asym_binned.tsv / subst_asym.pvals / subst_asym.tsv`: These files are not generally needed for NanoSeq analysis. They were originally used to detect asymmetries in the original DuplexSeq & BotSeqS protocols. 
* `mismatches.trinuc-profile.pdf / mismatches.subst_asym.pdf / mismatches.subst_asym.pvals / mismatches.subst_asym.tsv`: These files show the asymmetries and pyrimidine/purine-based trinucleotide substitution profiles for single-strand consensus calls. These profiles are useful to understand DNA damage during library preparation.
* `DSC_errors_per_channel.pdf / DSC_estimated_error_rates.pdf / estimated_error_rates.tsv / SSC-mismatches-Both.triprofiles.tsv / SSC-mismatches-Purine.triprofiles.tsv / SSC-mismatches-Pyrimidine.triprofiles.tsv`: Based on the independent error rates in the purine and pyrimidine channels (e.g. G>T and C>A), we calculate the probability of having independent errors affecting both strands and resulting in double-strand consensus.


