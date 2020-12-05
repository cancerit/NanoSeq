# NanoSeq

NanoSeq is a protocol based on Duplex Sequencing ([Schmitt et al, 2012](https://doi.org/10.1073/pnas.1208715109)) and BotSeqS ([Hoang et al, 2016](https://doi.org/10.1073/pnas.1607794113)) that allows detection of low occurrence mutations with a high degree of confience. (add citation NanoSeq).

Note that some of the code might use the term BotSeq or NanoSeq interchangebly due to legacy reasons but all the code in this repository is intended for use in the analysis of NanoSeq experiments.

## Quick installation

```
./setup.sh path_to_install_to
export PATH=$PATH:path_to_install_to/bin
Rscript ./build/libInstall.R <default R library path>
```
Samtools must be in execution path for use of the wrapper script.

## General workflow

For a bulk/duplex (normal/tumour) pair of samples an analysis would require:

1) Mapping of trimmed fastq files for each sample (`extract_tags.py`)
2) Appending RB: tag (read bundle) to reads (`bamaddreadbundles`)
3) Contamination check of each sample (`verifyBAMId`)
4) Generation of tables files (`dsa`)
5) Generation of variant files (`variants`)
6) Summarizing of results (`variantcaller.R` & `botseq_results_plotter.R`)

A wrapper script is provided to facilitate parallel execution of steps 4,5 & 6.

## Mapping of NanoSeq reads

Prior to mapping fastq files must be pre-processed with `extract_tags.py` in order to trim adapter sequences and to add the appropiate tags (rb,mb) to the read headers.
```
#(trim 3 bases, skip 4 bases, add rb & mb tags)
python extract-tags.py -a 1_R1.fastq -b 1_R2.fastq -c 70#1R1.fastq -d 70#1R2.fastq -m 3 -s 4 -l 151

#(align with bwa appending rb & mb tags)
bwa mem -C hs37d5.fa 70#1R1.fastq.gz 70#1R2.fastq > 70#1.sam
```

## Appending read bundle tag

The read bundle tag ( RB: ) must be appended to each read-pair of a BAM. The tag consists of : read coordinate, mate corrdinate, read rb tag, mate mb tag (RG:rc,mc,rb,mb).

This can be accomplished doing the following:

```
#(sort sam with biobambam, append rc & mc tags)
bamsormadup inputformat=sam rcsupport=1 threads=1 < 70#1.sam > 70#1.bam

#(append RB tag)
bamaddreadbundles -I 70#1.bam -O 70#1.tag.bam
```

## Processing of normal

The NanoSeq analysis uses a normal/tumour sample pair. If the normal happens to be a NanoSeq experiment it must be processed further as to just keep one read-pair from each read bundle to produce a 'neat' normal.

```
randomreadinbundle -I 70#1.tag.bam -O 70#1.neat.bam
```

## Contamination check

It is highly recommended to carry out a contamination check of the sample pair with [`verifyBAMId`](https://github.com/Griffan/VerifyBamID). This contamination check must be done on a bam generated with `randomreadinbundle`, where only one read per read bundle is kept in the bam (see above).
An alpha < 0.005 would be acceptable for most situations.

## NanoSeq analysis

The wrapper script `runBotSeq.py` provides a convinient way of running the NanoSeq analyisis. It wraps calls to `dsa` (generation of tables), `variants` (generation of variant files), `variantcaller.R` & `botseq_results_plotter.R` (creation of summary csv files, plots and vcf file). All of these programs should be present in the execution path.



A typical execution with 8 threads would be:

```
runBotSeq.py -t 8 -A 70#1.neat.bam -B 80#1.bam -C SNP.sorted.bed.gz -D NOISE.sorted.bed.gz -R hs37d5.fa -b 0 -o ./output
```
Where 
-   -A is the normal (bulk) BAM
-   -B is the tumour (duplex) BAM
-   -C BED SNP mask. Common SNPs (need to add link to file)
-   -D BED Noise mask. Sites that are unreliable or artefactual (need to add link to file)
-   -R reference fasta
-   -b optimal parameter for variantcaller for a neat normal ( -b 5 for WGS normal)
-   -o output directory 

The wrapper can pass several advanced parameters to the variantcaller program. The default values have been found to be optimal for most cases.

** Note to assure correct pre-processing, duplex BAMs must have @PG tags for bamsormadup , bammarkduplicatesopt & bamaddreadbundles. A neat bulk (NanoSeq library) BAM must have a @PG tag for bamaddreadbundles & randomreadinbundle. A WGS bulk must NOT have a tag for bamaddreadbundles **

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
