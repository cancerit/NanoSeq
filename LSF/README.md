Sample bash script for use in with an LSF job sheduler. This launches an entire NanoSeq analysis. 

The following should be edited with the desired values:

```
NANOSEQ_PATH="~/NanoSeq/opt/bin"  # location of NanoSeq scripts and executables
TUMOUR="filtered-3.bam"           # Duplex (tumour) BAM with index
NORMAL="3.final.bam"              # Reference (normal) BAM with index
REF="genome.fa"                   # Genome reference with index
SNP_BED="SNP.sorted.bed.gz"       # SNP BED file
NOISE_BED="NOISE.sorted.bed.gz"   # NOISE BED file
CPU=80                            # Number of CPUs used to run dsa,var and indel steps
NAME="MY_NAME"                    # Name given to all the jobs in the LSF queue
```
