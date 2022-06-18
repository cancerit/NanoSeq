nextflow.enable.dsl=2

MAX_IN_PARALLEL = 30
QUEUE = "long"

process BWAMEM2_INDEX {

    container params.bwa_container

    input:
        path fasta
        val species
        val assembly

    output:
        path "bwamem2"      , emit: index
        path "versions.yml" , emit: versions


    maxForks MAX_IN_PARALLEL
    queue QUEUE

    maxRetries 2
    cpus 1
    memory { task.exitStatus == 130  ? 80.GB * task.attempt : 80.GB }

    script:
        """
        touch ${task.process}
        mkdir bwamem2
        bwa-mem2 index $fasta -p bwamem2/genome.fa
        samtools dict -a $assembly -s $species $fasta > bwamem2/genome.fa.dict
        ln -s ../${fasta} ./bwamem2 
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        """
        mkdir bwamem2
        touch bwamem2/
        touch bwamem2/genome.fa.0123
        touch bwamem2/genome.fa.ann
        touch bwamem2/genome.fa.pac
        touch bwamem2/genome.fa.amb
        touch bwamem2/genome.fa.bwt.2bit.64
        touch bwamem2/${fasta}.dict
        ln -s ../${fasta} ./bwamem2/genome.fa
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process BWAMEM2_MAP {

    tag "$meta.name"

    container params.bwa_container

    input:
        tuple val(meta), path(reads)
        path index_dir

    output:
        tuple val(meta), path("out/${meta.name}.bam"),path("out/${meta.name}.bam.bai"), emit: bam
        path  "versions.yml", emit: versions

    maxForks MAX_IN_PARALLEL
    queue QUEUE
    clusterOptions "-R avx512 " //ask for latest vector instructions

    cpus 6
    maxRetries 4
    memory { task.exitStatus == 130 ? 40.GB * task.attempt : 40.GB }

    script:
        def args = task.ext.args ?: '-C'
        def args2 = task.ext.args2 ?: '-n'
        """
        touch ${task.process}_${meta.name}
        mkdir out
        bwa-mem2 mem $args -t $task.cpus ${index_dir}/genome.fa $reads \\
            | samtools sort -@ $task.cpus $args2 -o out/${meta.name}.bam -
        touch out/${meta.name}.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        """
        mkdir out
        touch out/${meta.name}.bam
        touch out/${meta.name}.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process BWAMEM2_REMAP {

    tag "$meta.name"

    container params.bwa_container

    input:
        tuple val(meta), path(bam)
        path index_dir

    output:
        tuple val(meta), path("sort/${meta.name}.bam"), path("sort/${meta.name}.bam.bai"), emit: bam
        path  "versions.yml", emit: versions

    maxForks MAX_IN_PARALLEL
    queue QUEUE
    clusterOptions "-R avx512 " //ask for latest vector instructions
    
    maxRetries 4
    cpus 6
    memory { task.exitStatus == 130  ? 30.GB * task.attempt : 30.GB }

    script:
        def args = task.ext.args ?: '-n -tags BC,QT,mb,rb -b \'-T 30 -Y\''
        def args2 = task.ext.args ?: '-n'
        """
        touch ${task.process}_${meta.name}
        mkdir bwamem
        mkdir sort
        bwa_mem.pl $args -bwamem2 -t $task.cpus -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $bam
        samtools sort -@ $task.cpus $args2 ./bwamem/${meta.name}.bam -o ./sort/${meta.name}.bam
        rm -rf ./bwamem/*
        touch sort/${meta.name}.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
            bwa_mem.pl: \$(bwa_mem.pl -v )
        END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        """
        mkdir bwamem
        mkdir sort
        touch ./sort/${meta.name}.bam
        touch ./sort/${meta.name}.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
            bwa_mem.pl: \$(bwa_mem.pl -v )
        END_VERSIONS
        """
}

