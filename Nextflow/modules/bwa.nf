nextflow.enable.dsl=2

MAX_IN_PARALLEL = 25
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

    maxRetries 2
    cpus 1
    memory { task.exitStatus == 130  ? 80.GB * task.attempt : 80.GB }

    script:
        """
        touch ${task.process}
        mkdir -p bwamem2
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
        mkdir -p bwamem2
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
    clusterOptions "-R avx512" //ask for a node with latest vector instructions

    cpus 10
    maxRetries 4
    memory { ( task.exitStatus == 130 || task.exitStatus == 140) ? 50.GB * task.attempt : 50.GB }
    queue { task.exitStatus == 140 ? "basement" : QUEUE }

    script:
        def args = task.ext.args ?: '-C'
        def args2 = task.ext.args2 ?: '-n'
        """
        touch ${task.process}_${meta.name}
        mkdir -p out
        bwa-mem2 mem $args -t $task.cpus ${index_dir}/genome.fa $reads \\
            | samtools sort -@ $task.cpus $args2 -m 2G -o out/${meta.name}.bam -
        touch out/${meta.name}.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        """
        mkdir -p out
        touch out/${meta.name}.bam
        touch out/${meta.name}.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

//divided bwa_mem.pl (remap) into two tasks to assign resources more efficiently
process REMAP_SPLIT {

    tag "$meta.name"

    container params.bwa_container

    input:
        tuple val(meta), path(bam)
        path index_dir

    output:
        tuple val(meta), path("bwamem"), path(bam), emit: split_out
    
    maxRetries 4
    cpus 3
    memory { task.exitStatus == 130  ? 6.GB * task.attempt : 6.GB }

    script:
        def args = task.ext.args ?: '-n -tags BC,QT,mb,rb -b \'-T 30 -Y\''
        """
        touch ${task.process}_${meta.name}
        mkdir -p bwamem
        bwa_mem.pl -p setup $args -bwamem2 -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $bam
        bwa_mem.pl -p split -t $task.cpus $args -bwamem2 -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $bam
        """

    stub:
        def args = task.ext.args ?: ''
        """
        mkdir -p bwamem
        touch ./bwamem/${meta.name}.bam
        touch ./bwamem/${meta.name}.bam.bai
        """
}

process REAMAP_BWAMEM2 {

    tag "$meta.name"

    container params.bwa_container

    input:
        tuple val(meta), path(bwamem), path(bam)
        path index_dir

    output:
        tuple val(meta), path("sort/${meta.name}.bam"), path("sort/${meta.name}.bam.bai"), emit: bam
        path  "versions.yml", emit: versions

    maxForks MAX_IN_PARALLEL
    clusterOptions "-R avx512 " //ask for a node with latest vector instructions
    
    maxRetries 4
    cpus 10
    memory { ( task.exitStatus == 130 || task.exitStatus == 140) ? 50.GB * task.attempt : 50.GB }
    queue { task.exitStatus == 140 ? "basement" : QUEUE }

    script:
        def args = task.ext.args ?: '-n -tags BC,QT,mb,rb -b \'-T 30 -Y\''
        def args2 = task.ext.args ?: '-n'
        """
        touch ${task.process}_${meta.name}
        mkdir -p sort
        bwa_mem.pl -p bwamem $args -bwamem2 -t $task.cpus -mt $task.cpus -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $bam
        #need the mark process so the final bam file gets placed in the correct location
        bwa_mem.pl -p mark -n $args -bwamem2 -t $task.cpus -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $bam
        samtools sort -@ $task.cpus $args2 -m 2G -o ./sort/${meta.name}.bam ./bwamem/${meta.name}.bam
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
        def args2 = task.ext.args ?: ''
        """
        mkdir -p bwamem
        mkdir -p sort
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

workflow BWAMEM2_REMAP {
    take :
        bam_in
        index_dir

    main :
        REMAP_SPLIT(bam_in, index_dir)
        REAMAP_BWAMEM2(REMAP_SPLIT.out.split_out, index_dir)

    emit :
        versions = REAMAP_BWAMEM2.out.versions
        bam = REAMAP_BWAMEM2.out.bam
}
