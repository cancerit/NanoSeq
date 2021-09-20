# NanoSeq

NanoSeq is a protocol based on Duplex Sequencing ([Schmitt et al, 2012](https://doi.org/10.1073/pnas.1208715109)) and BotSeqS ([Hoang et al, 2016](https://doi.org/10.1073/pnas.1607794113)) that allows detection of low occurrence mutations with a high degree of confience. ([Abascal et al, 2021](https://doi.org/10.1038/s41586-021-03477-4)).

## Dependencies

Execution of the scripts from this repo requires that these dependencies be on PATH :

* samtools
* bcftools
* Rscript

## Installation

```
./setup.sh path_to_install                       #install code from this repo
export PATH=$PATH:path_to_install/bin
Rscript ./build/libInstall.R <R libraries path>  #install all the required R libraries
```

## FASTQ and BAM processing of NanoSeq libraries

1) Mapping of trimmed fastq files for each sample (`extract_tags.py`)
2) Appending RB: tag (read bundle) to reads (`bamaddreadbundles`)

### Mapping of NanoSeq reads

Prior to mapping fastq files must be pre-processed with `extract_tags.py` in order to trim adapter sequences and to add the appropiate tags (rb,mb) to the read headers.
```
#(trim 3 bases, skip 4 bases, add rb & mb tags)
python extract-tags.py -a 1_R1.fastq -b 1_R2.fastq -c 70#1R1.fastq -d 70#1R2.fastq -m 3 -s 4 -l 151

#(align with bwa appending rb & mb tags)
bwa mem -C hs37d5.fa 70#1R1.fastq.gz 70#1R2.fastq > 70#1.sam
```

### Appending read bundle tag

The read bundle tag ( RB: ) must be appended to each read-pair of a BAM. The tag consists of : read coordinate, mate corrdinate, read rb tag, mate mb tag (RG:rc,mc,rb,mb).

This can be accomplished doing the following:

```
#(sort sam with biobambam, append rc & mc tags)
bamsormadup inputformat=sam rcsupport=1 threads=1 < 70#1.sam > 70#1.bam

#(append RB tag)
bamaddreadbundles -I 70#1.bam -O 70#1.tag.bam
```

### Processing of normal

The NanoSeq analysis uses a normal/tumour sample pair. If the normal happens to be a NanoSeq experiment it must be processed further as to just keep one read-pair from each read bundle to produce a 'neat' normal.

```
randomreadinbundle -I 70#1.tag.bam -O 70#1.neat.bam
```
### Note

Correct pre-processing means that duplex BAMs must have @PG tags for bamsormadup , bammarkduplicatesopt & bamaddreadbundles. A neat bulk (NanoSeq library) BAM must have a @PG tag for bamaddreadbundles & randomreadinbundle. A WGS bulk will NOT have a tag for bamaddreadbundles

### Contamination check

It is highly recommended to carry out a contamination check of the sample pair with [`verifyBAMId`](https://github.com/Griffan/VerifyBamID). This contamination check must be done on a bam generated with `randomreadinbundle`, where only one read per read bundle is kept in the bam (see above).
An alpha < 0.005 would be acceptable for most situations.


## NanoSeq analysis

For a bulk/duplex (normal/tumour) pair of samples an analysis requires the following steps:

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
  -a 2 \
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

Merge final files, produce summaries, compute efficiency. Results can be found in here.

```
runNanoSeq.py -t 2 \
  -A normal.bam \
  -B tumour.bam \
  -R genome.fa \
  post
```
### Note

A bash script is provided in the LSF directory to facilitate generation and running of all the steps with LSF job arrays. Script should be edited for each analyisis.

## Output

Intermediate files are placed in the tables and variants. The resulting files and plots are placed in the summary folder.

The most relevant summary files include the following.

- **muts.vcf**

    This file provides the standard variant information with the following flags for each variant:
    
    - PASS these are SNVs of interest
    - SNVs occurring at common SNP sites (DBSNP mask)
    - SHEARWATER marks error-prone sites defined from a large panel of normals

    Note that contaminated samples will usually show many mutations that overlap with DBSNP sites.

- **burden.masked-vs-unmasked.pdf**

    Provides a qualitative look at contamination. Burdens should be similar for both bars.

- **mut_burden.tsv**

    Summarizes the mutation burden estimates, observed and corrected. Burdens are expressed as number of mutations per base. The correction is done to account for differences between the trinucleotide frequencies of the genome and of NanoSeq data.

- **trint_counts_and_ratio2genome.tsv**

    File providing information on the frequency of each 32 pyrimidine trinucleotide in the NanoSeq data and the ratio to the genome frequency.

    tri_bg is the trinucleotide frequency. ratio2genome is the trinucleotide frequency vs that seen in the genome

- **subs_obs_corrected**

    Shows the correction for each of the 96 trinucleotides substitutions, using the results in **trint_counts_and_ratio2genome.tsv**.

- **trinuc-profiles.pdf**

    Shows the observed substitution profile, the corrected substitution profile, and the actual mutation rates


- **subst_asym**

    These files are not generally needed for NanoSeq analysis. They have been employed to detect asymetries in the end repair of the original BotSeq protocol.
