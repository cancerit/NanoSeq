#!/bin/bash

NANOSEQ_PATH="/nfs/users/nfs_r/ra11/lustre2/NanoSeq/NanoSeq/opt/bin"

TUMOUR="filtered-3.bam"
NORMAL="3.final.bam"
REF="genome.fa"
SNP_BED="SNP.sorted.bed.gz"
NOISE_BED="NOISE.sorted.bed.gz"
CPU=60
NAME="NANO3"

if [ ! -f $TUMOUR ] || [ ! -f $TUMOUR.bai ] ; then
  echo "TUMOUR file or index not found!"
  exit 1
fi

if [ ! -f $NORMAL ] || [ ! -f $NORMAL.bai ]; then
  echo "NORMAL files or index not found!"
  exit 1
fi

if [ ! -f $REF ] || [ ! -f $REF.fai ]; then
  echo "REF file or index not found!"
  exit 1
fi

if [ ! -f $SNP_BED ] || [ ! -f $SNP_BED.tbi ]; then
  echo "SNP_BED file or index not found!"
  exit 1
fi

if [ ! -f $NOISE_BED ] || [ ! -f $NOISE_BED.tbi ]; then
  echo "NOISE_BED file or index not found!"
  exit 1
fi



DEPENDS="module add samtools; module add bcftools; module add ISG/R/4.1.0"

#cov
cat > run_cov.bsub <<- EOF
#!/bin/bash
#BSUB -q normal
#BSUB -J $NAME.cov[1-10]
#BSUB -n 1
#BSUB -R "select[mem>900] rusage[mem=900]" -M900
#BSUB -e %I.cov.err
#BSUB -o %I.cov.out

$DEPENDS

export PATH=/nfs/users/nfs_r/ra11/lustre2/NanoSeq/NanoSeq/opt/bin:\$PATH

runNanoSeq.py -k 10 -j \$LSB_JOBINDEX  \\
  -A $NORMAL \\
  -B $TUMOUR \\
  -R $REF \\
  cov \\
  -Q 0 \\
  --exclude "MT,GL%,NC_%,hs37d5"
EOF

#part
cat > run_part.bsub <<- EOF
#!/bin/bash
#BSUB -q normal
#BSUB -J $NAME.part
#BSUB -n 1
#BSUB -R "select[mem>8000] rusage[mem=8000]" -M8000
#BSUB -e part.err
#BSUB -o part.out

$DEPENDS

export PATH=/nfs/users/nfs_r/ra11/lustre2/NanoSeq/NanoSeq/opt/bin:\$PATH

runNanoSeq.py -t 1  \\
  -A $NORMAL \\
  -B $TUMOUR \\
  -R $REF \\
  part \\
  -n $CPU \\
# --excludedBED BED 
EOF

#dsa
cat > run_dsa.bsub <<- EOF
#!/bin/bash
#BSUB -q normal
#BSUB -J $NAME.dsa[1-$CPU]
#BSUB -n 1
#BSUB -R "select[mem>900] rusage[mem=900]" -M900
#BSUB -e %I.dsa.err
#BSUB -o %I.dsa.out

$DEPENDS

export PATH=/nfs/users/nfs_r/ra11/lustre2/NanoSeq/NanoSeq/opt/bin:\$PATH

runNanoSeq.py -k $CPU -j \$LSB_JOBINDEX  \\
  -A $NORMAL \\
  -B $TUMOUR \\
  -R $REF \\
  dsa \\
  -C $SNP_BED \\
  -D $NOISE_BED \\
  -d 2 \\
  -q 30
EOF

#var
cat > run_var.bsub <<- EOF
#!/bin/bash
#BSUB -q normal
#BSUB -J $NAME.var[1-$CPU]
#BSUB -n 1
#BSUB -R "select[mem>900] rusage[mem=900]" -M900
#BSUB -e %I.var.err
#BSUB -o %I.var.out

$DEPENDS

export PATH=/nfs/users/nfs_r/ra11/lustre2/NanoSeq/NanoSeq/opt/bin:\$PATH

runNanoSeq.py -k $CPU -j \$LSB_JOBINDEX  \\
  -A $NORMAL \\
  -B $TUMOUR \\
  -R $REF \\
  var \\
  -a 2 \\
  -b 5 \\
  -c 0 \\
  -d 2 \\
  -f 0.9 \\
  -i 1 \\
  -m 8 \\
  -n 3 \\
  -p 0 \\
  -q 60 \\
  -r 144 \\
  -v 0.01 \\
  -x 8 \\
  -z 12
EOF

#indel
cat > run_indel.bsub <<- EOF
#!/bin/bash
#BSUB -q normal
#BSUB -J $NAME.indel[1-$CPU]
#BSUB -n 1
#BSUB -R "select[mem>900] rusage[mem=900]" -M900
#BSUB -e %I.indel.err
#BSUB -o %I.indel.out

$DEPENDS

export PATH=/nfs/users/nfs_r/ra11/lustre2/NanoSeq/NanoSeq/opt/bin:\$PATH

runNanoSeq.py -k $CPU -j \$LSB_JOBINDEX  \\
  -A $NORMAL \\
  -B $TUMOUR \\
  -R $REF \\
  indel \\
  -s sample \\
  --rb 2 \\
  --t3 135 \\
  --t5 10 \\
  --mc 16
EOF

#post
cat > run_post.bsub <<- EOF
#!/bin/bash
#BSUB -q normal
#BSUB -J $NAME.post
#BSUB -n 2
#BSUB -R "select[mem>5000] rusage[mem=5000]" -M5000
#BSUB -e post.err
#BSUB -o post.out

$DEPENDS

export PATH=/nfs/users/nfs_r/ra11/lustre2/NanoSeq/NanoSeq/opt/bin:\$PATH

runNanoSeq.py -t 2  \\
  -A $NORMAL \\
  -B $TUMOUR \\
  -R $REF \\
  post
EOF

function get_jobid {
  OUT=$*
  ID=$(echo $OUT | head -n1 | cut -d'<' -f2 | cut -d'>' -f1)
}

echo "Submit jobs to queue? (y/n)"
read ANSWER


if [ $ANSWER = 'y' ]; then
  get_jobid `bsub < run_cov.bsub`
  COV_ID=$ID
  echo "submitted cov as $COV_ID"

  get_jobid `bsub -w "done($COV_ID)" < run_part.bsub`
  PART_ID=$ID
  echo "submitted part as $PART_ID"

  get_jobid `bsub -w "done($PART_ID)" < run_dsa.bsub`
  DSA_ID=$ID
  echo "submitted dsa as $DSA_ID"

  get_jobid `bsub -w "done($DSA_ID)" < run_var.bsub`
  VAR_ID=$ID
  echo "submitted var as $VAR_ID"

  get_jobid `bsub -w "done($DSA_ID)" < run_indel.bsub`
  INDEL_ID=$ID
  echo "submitted indel as $INDEL_ID"

  get_jobid `bsub -w "done($INDEL_ID) && done($VAR_ID)" < run_post.bsub`
  POST_ID=$ID
  echo "submitted post as $ID"
fi

