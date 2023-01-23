nextflow.enable.dsl=2

MAX_IN_PARALLEL = 25
QUEUE = "long"

process BWAMEM2_INDEX {

    container params.bwa_image

    input:
        path fasta
        val species
        val assembly

    output:
        path "bwamem2"      , emit: index
        path "versions.yml" , emit: versions


    maxForks MAX_IN_PARALLEL

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

    container params.bwa_image

    input:
        tuple val(meta), path(reads)
        path index_dir
        val min_mapQ

    output:
        tuple val(meta), path("out/${meta.name}.cram"),path("out/${meta.name}.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    maxForks MAX_IN_PARALLEL
    clusterOptions "-R avx512" //ask for a node with latest vector instructions

    cpus 10
    errorStrategy { sleep(Math.pow(2, task.attempt) * 10000 as long); return 'retry' } //delayed retry
    memory { ( task.exitStatus == 130 || task.exitStatus == 140) ? 60.GB * task.attempt : 60.GB }
    queue { task.exitStatus == 140 ? "basement" : QUEUE }

    script:
        def args = task.ext.args ?: '-C'
        def args2 = task.ext.args2 ?: '-n'
        """
        set -o pipefail
        touch ${task.process}_${meta.name}
        mkdir -p out
        bwa-mem2 mem $args -t $task.cpus ${index_dir}/genome.fa $reads \\
            | samtools sort -@ $task.cpus $args2 -O bam -l 0 -m 2G - \\
            | samtools view -@ $task.cpus -q $min_mapQ  -T ${index_dir}/genome.fa -o out/${meta.name}.cram
        touch out/${meta.name}.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        """
        sleep \$[ ( \$RANDOM % 20 )  + 1 ]s
        mkdir -p out
        touch out/${meta.name}.cram
        touch out/${meta.name}.cram.crai
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

    container params.bwa_image

    input:
        tuple val(meta), path(cram), path(crai)
        path index_dir

    output:
        tuple val(meta), path("bwamem"), path(cram), emit: split_out
    
    cpus 3
    memory { task.exitStatus == 130  ? 15.GB * task.attempt : 15.GB }

    script:
        def args = task.ext.args ?: '-n -tags BC,QT,mb,rb -b \'-T 30 -Y\''
        """
        touch ${task.process}_${meta.name}
        mkdir -p bwamem
        bwa_mem.pl -p setup $args -bwamem2 -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $cram
        bwa_mem.pl -p split -t $task.cpus $args -bwamem2 -cram -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $cram
        """

    stub:
        def args = task.ext.args ?: ''
        """
        sleep \$[ ( \$RANDOM % 20 )  + 1 ]s
        mkdir -p bwamem
        touch ./bwamem/${meta.name}.cram
        touch ./bwamem/${meta.name}.cram.crai
        """
}

process REAMAP_BWAMEM2 {

    tag "$meta.name"

    container params.bwa_image

    input:
        tuple val(meta), path(bwamem), path(cram)
        path index_dir
        val min_mapQ

    output:
        tuple val(meta), path("sort/${meta.name}.cram"), path("sort/${meta.name}.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    maxForks MAX_IN_PARALLEL
    clusterOptions "-R avx512 " //ask for a node with latest vector instructions
    
    errorStrategy { sleep(Math.pow(2, task.attempt) * 10000 as long); return 'retry' } //delayed retry
    cpus 10
    memory { ( task.exitStatus == 130 || task.exitStatus == 140) ? 50.GB * task.attempt : 50.GB }
    queue { task.exitStatus == 140 ? "basement" : QUEUE }

    script:
        def args = task.ext.args ?: "-n -tags BC,QT,mb,rb -b \'-T $min_mapQ -Y\'"
        def args2 = task.ext.args ?: '-n'
        """
        touch ${task.process}_${meta.name}
        mkdir -p sort
        bwa_mem.pl -p bwamem $args -bwamem2 -cram -t $task.cpus -mt $task.cpus -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $cram
        #need the mark process so the final cram file gets placed in the correct location
        bwa_mem.pl -p mark -n $args -bwamem2  -cram -t $task.cpus -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $cram
        bwa_mem.pl -p stats -n $args -bwamem2  -cram -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.name} $cram
        samtools sort -@ $task.cpus $args2 -O cram -m 2G -o ./sort/${meta.name}.cram ./bwamem/${meta.name}.cram
        touch sort/${meta.name}.cram.crai
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
        sleep \$[ ( \$RANDOM % 20 )  + 1 ]s
        mkdir -p bwamem
        mkdir -p sort
        touch ./sort/${meta.name}.cram
        touch ./sort/${meta.name}.cram.crai
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
        cram_in
        index_dir
        min_mapQ

    main :
        REMAP_SPLIT(cram_in, index_dir)
        REAMAP_BWAMEM2(REMAP_SPLIT.out.split_out, index_dir, min_mapQ)

    emit :
        versions = REAMAP_BWAMEM2.out.versions
        cram = REAMAP_BWAMEM2.out.cram
}
