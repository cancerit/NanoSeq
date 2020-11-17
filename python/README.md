## Helper python script for BotSeq

# check_lsf_batch.py

Check completion of a batch of LSF jobs

# extract_tags.py 

# gridsearch_merge.py

# intervals.py

Generate a BED of genomic intervals for the GRC37 reference

# variantcaller_merge.py

Merge the results of variantcaller and produce output csv files

# runBotSeq.py

Automate a BotSeq execution. It can be used to launch calculations in a LSF job array using -j and -t options. (see run.bjob script) or tentatively it can be used to run calculations in a multithreaded machine using the -t option.

Sample LSF job array

```
#!/bin/bash
##run a 300 CPU BotSeq job using an LSF array

#BSUB -q normal
#BSUB -J BotSeq[01-300]
#BSUB -R "span[hosts=1] select[mem>200] rusage[mem=200]" -M200
#BSUB -e %I.err
#BSUB -o %I.out
 
#add samtools to the path
module add samtools
 
#add the location of BotSeq executables to the path
export PATH=/lustre/scratch119/realdata/mdt1/team78/ra11/BotSeqS/code/botseq/opt/bin/:$PATH
 
runBotSeq.py -j $LSB_JOBINDEX -t 300 -A ../BotSeqS/BOTSEQ_BAMS/10.final.neat.bam -B ../BotSeqS/BOTSEQ_BAMS/filtered-10.bam -C SNP.sorted.bed.gz -R /lustre/scratch117/casm/team78/ro4/hs37d5/hs37d5.fa -D NOISE.sorted.bed.gz -b 0
```

The script requires that the BotSeq executables and samtools are in the user's path.

