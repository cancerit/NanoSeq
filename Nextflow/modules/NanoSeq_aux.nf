nextflow.enable.dsl=2

//Run this for fastq files prior to mapping
process ADD_NANOSEQ_FASTQ_TAGS {
    
    tag "${meta.id}_${meta.type}"

    container params.nanoseq_image

    input:
        tuple val(meta), path(reads)
        val m
        val s

    output:
        tuple val(meta), path("out/*.fastq.gz"), emit: fastqs
        path  "versions.yml", emit: versions

    maxRetries 4
    cpus 1
    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { ( task.exitStatus == 130 ) ? 500.MB * task.attempt : 500.MB }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 500.MB * task.attempt : 500.MB }

    script:
        def read1 = reads[0]
        def read2 = reads[1]
        """
        touch ${task.process}_${meta.id}_${meta.type}
        mkdir -p out
        L=`zcat $read1 | head -2 | tail -1 | awk '{ print length }'`
        extract_tags.py -a $read1 -b $read2 -c ./out/${meta.name}_R1.fastq.gz -d ./out/${meta.name}_R2.fastq.gz -m $m -s $s -l \$L
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            extract-tags.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """

    stub:
        def read1 = reads[0]
        def read2 = reads[1]
        """
        mkdir -p out
        touch ./out/${meta.name}_R1.fastq.gz
        touch ./out/${meta.name}_R2.fastq.gz
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            extract-tags.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """
}

process MARKDUP {

    tag "${meta.id}_${meta.type}"

    container params.bwa_image

    input:
        tuple val(meta), path(cram), path(index)
        path(ref_path)

    output:
        tuple val(meta), path("optdup/${meta.name}.cram"), path("optdup/${meta.name}.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    maxRetries 4
    cpus 5
    // *MODIFIED* (ao7): increased memory; added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    // memory { ( task.exitStatus == 130 || task.exitStatus == 140) ? 25.GB * task.attempt : 25.GB }
    // queue { task.exitStatus == 140 ? "long" : "normal" }
    memory { ( task.exitStatus == 130 || task.exitStatus == 140 || task.exitStatus == Integer.MAX_VALUE ) ? 40.GB * task.attempt : 40.GB }
    queue { ( task.exitStatus == 140 || task.exitStatus == Integer.MAX_VALUE ) ? "long" : "normal" }
    //

    script:
        """
        touch ${task.process}_${meta.id}_${meta.type}
        mkdir -p nsorted
        mkdir -p optdup
        ln -s ../$cram ./nsorted/${meta.name}.cram
        samtools view -H $cram | grep SO:queryname > /dev/null || \\
            ( rm ./nsorted/${meta.name}.cram; samtools sort -@ $task.cpus -O cram -m 2G -n -o ./nsorted/${meta.name}.cram $cram )
        samtools view --no-PG -@ $task.cpus -h ./nsorted/${meta.name}.cram | bamsormadup inputformat=sam level=0 blocksortverbose=0 rcsupport=1 threads=$task.cpus fragmergepar=$task.cpus optminpixeldif=10 | \\
            bammarkduplicatesopt verbose=0 level=0 index=0 optminpixeldif=2500 | samtools view --no-PG -@ $task.cpus -C -o ./optdup/${meta.name}.cram
        samtools index -@ $task.cpus ./optdup/${meta.name}.cram
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bamsormadup: \$(bamsormadup -v 2>&1 | head -1 | sed 's/.*version //' | sed 's/\\.\$// ')
            bammarkduplicatesopt: \$(bammarkduplicatesopt -v 2>&1 | head -1 | sed 's/.*version //' | sed 's/\\.\$// ')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        """
        mkdir -p nsorted
        mkdir -p optdup
        ln -s ../$cram ./nsorted/${meta.name}.cram
        touch ./optdup/${meta.name}.cram
        touch ./optdup/${meta.name}.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bamsormadup: \$(bamsormadup -v 2>&1 | head -1 | sed 's/.*version //' | sed 's/\\.\$// ')
            bammarkduplicatesopt: \$(bammarkduplicatesopt -v 2>&1 | head -1 | sed 's/.*version //' | sed 's/\\.\$// ')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process NANOSEQ_ADD_RB {

    tag "${meta.id}_${meta.type}"
    
    container params.nanoseq_image

    input:
        tuple val(meta), path(cram), path(index)
        path( ref_path )

    output:
        tuple val(meta), path("out/${meta.name}.cram"), path("out/${meta.name}.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    maxRetries 4
    cpus 1
    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { task.exitStatus == 130 ? 2.GB * task.attempt : 2.GB }
    //queue { task.exitStatus == 140 ? "long" : "normal" }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 2.GB * task.attempt : 2.GB }
    queue { ( task.exitStatus == 140 || task.exitStatus == Integer.MAX_VALUE ) ? "long" : "normal" }
    //

    script:
        """
        touch ${task.process}_${meta.id}_${meta.type}
        mkdir -p out
        NLINES=`samtools view $cram | head -100 | grep rb: | grep rc: | grep mb: | grep mc: | wc -l` || true
        ## *MODIFIED* (ao7): fixed numeric comparison operator
        ## if [ \$NLINES != 0 ]; then
        if [ \$NLINES -ne 0 ]; then
        ##
            bamaddreadbundles -I $cram -O ./out/${meta.name}.cram
            samtools index ./out/${meta.name}.cram
        else
          ln -s ../$cram ./out/${meta.name}.cram
          ln -s ../$index ./out/${meta.name}.cram.crai
        fi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bamaddreadbundles: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    stub:
        """
        mkdir -p out
        touch ./out/${meta.name}.cram
        touch ./out/${meta.name}.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bamaddreadbundles: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process NANOSEQ_DEDUP {

    tag "${meta.id}_${meta.type}"
    
    container params.nanoseq_image

    input:
        tuple val(meta), path(cram), path(index)
        path (ref_path)
        val m

    output:
        tuple val(meta), path("out/${meta.name}.neat.cram"), path("out/${meta.name}.neat.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    maxRetries 4
    cpus 1
    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { task.exitStatus == 130 ? 10.GB * task.attempt : 10.GB }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 10.GB * task.attempt : 10.GB }
    //

    script:
        """
        touch ${task.process}_${meta.id}_${meta.type}
        mkdir -p out
        NLINES1=`samtools view -H $cram | grep ^@PG | grep ID:bamaddreadbundles | wc -l` || true
        NLINES2=`samtools view -H $cram | grep ^@PG | grep ID:randomreadinbundle | wc -l` || true
        ## *MODIFIED* (ao7): fixed numeric comparison operators
        ## if [ \$NLINES1 > 0 ] && [ \$NLINES2 == 0 ]; then
        if [ \$NLINES1 -gt 0 ] && [ \$NLINES2 -eq 0 ]; then
        ##
            ## *MODIFIED* (ao7): added -m parameter
            ##randomreadinbundle -I $cram -O ./dedup/${meta.name}.neat.cram
            randomreadinbundle -I $cram -O ./dedup/${meta.name}.neat.cram -m $m
            ##
            samtools index ./dedup/${meta.name}.neat.cram
        else
          ln -s ../$cram ./out/${meta.name}.neat.cram
          ln -s ../$index ./out/${meta.name}.neat.cram.crai
        fi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            randomreadinbundle: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    stub:
        """
        mkdir -p out
        touch ./out/${meta.name}.neat.cram
        touch ./out/${meta.name}.neat.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            randomreadinbundle: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process VERIFY_BAMID {

    tag "${meta.id}_${meta.type}"
    
    container params.nanoseq_image

    input:
        tuple val(meta), path(cram), path(index)
        val epsilon
        path ref_path
        path vb_ud
        path vb_bed
        path vb_mu

    publishDir "$params.outDir/outNextflow/$meta.id", mode: 'copy', pattern: "verifyBAMid/*", overwrite: true

    output:
        path "verifyBAMid/${meta.name}.verifyBAMid.txt", emit: verifybamid
        path  "versions.yml", emit: versions

    maxRetries 1
    cpus 2
    memory 6.GB

    script:
        """
        touch ${task.process}_${meta.id}_${meta.type}
        mkdir -p verifyBAMid
        #processing of CRAM is super slow so need to convert to BAM
        samtools view -b -@ $task.cpus  $cram > temp.bam 
        samtools index -@ $task.cpus temp.bam
        ( VerifyBamID --Epsilon $epsilon --UDPath $vb_ud \\
            --BedPath $vb_bed --MeanPath $vb_mu \\
            --Reference ${ref_path}/genome.fa --BamFile temp.bam > verifyBAMid/${meta.name}.verifyBAMid.txt 2> verifyBAMid/error ) || ( cp verifyBAMid/error verifyBAMid/${meta.name}.verifyBAMid.txt )
        rm temp.bam
        rm temp.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            VerifyBamID: \$(VerifyBamID -v 2>&1 | grep Version | sed 's/.*://')
        END_VERSIONS
        """
    stub:
        """
        mkdir -p verifyBAMid
        touch verifyBAMid/${meta.name}.verifyBAMid.txt
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            VerifyBamID: \$(VerifyBamID -v 2>&1 | grep Version | sed 's/.*://')
        END_VERSIONS
        """
}

MAXN = 2
process NANOSEQ_EFFI {
    
    tag "${meta.id}_${meta.type}"

    container params.nanoseq_image

    input:
        tuple val(meta), path(cram), path(index), path(cram_neat), path(index_neat)
        path ref_path

    publishDir "$params.outDir/outNextflow/$meta.id", mode: 'copy', pattern : "effi/*" , overwrite: true

    output:
        path "effi/${meta.name}.effi.tsv", emit: effi
        path "versions.yml", emit: versions

    maxRetries MAXN
    cpus 4
    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { task.exitStatus == 130  ? 25.GB * task.attempt : 25.GB }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 25.GB * task.attempt : 25.GB }
    //
    errorStrategy { task.attempt == MAXN ? 'ignore' : 'retry' }

    script:
        def cpus = task.cpus - 1
        """
        touch ${task.process}_${meta.id}_${meta.type}
        mkdir -p effi
        efficiency_nanoseq.pl -t $cpus -d $cram_neat -x $cram -o effi/${meta.name}.effi -r ${ref_path}/genome.fa
 
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            efficiency_nanoseq.pl : \$(runNanoSeq.py -v)
        END_VERSIONS
        """
    stub:
        """
        mkdir -p effi
        touch effi/${meta.name}.effi.tsv
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            efficiency_nanoseq.pl : \$(runNanoSeq.py -v)
        END_VERSIONS
        """
}

process NANOSEQ_VAF {

    tag "${meta.id}_${meta.type}"
    
    container params.nanoseq_image

    input:
        tuple val(meta), path(vcf_muts), path(index_muts), path(vcf_indel), path(index_indel),
             path(bed_cov), path(index_cov), path(cram_neat), path(index_neat)

    publishDir "$params.outDir/outNextflow/$meta.id", mode: 'copy', pattern: "*.vcf.*", overwrite: true

    output:
        tuple val(meta),path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi"), emit: vcf
        path  "versions.yml", emit: versions

    maxRetries 4
    cpus 2
    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { task.exitStatus == 130 ? 10.GB * task.attempt : 10.GB }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 10.GB * task.attempt : 10.GB }
    //

    script:
        """
        touch ${task.process}_${meta.id}_${meta.type}
        mkdir -p out
        snv_merge_and_vaf_calc.R $vcf_muts $vcf_indel $cram_neat $bed_cov out/${meta.id}.vcf
        bcftools sort -Oz out/${meta.id}.vcf -o ${meta.id}.vcf.gz
        bcftools index -t ${meta.id}.vcf.gz
        rm out/${meta.id}.vcf
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snv_merge_and_vaf_calc.R : \$(runNanoSeq.py -v)
        END_VERSIONS
        """
    stub:
        """
        mkdir out
        mkdir sorted
        touch ${meta.id}.vcf.gz
        touch ${meta.id}.vcf.gz.tbi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snv_merge_and_vaf_calc.R : \$(runNanoSeq.py -v)
        END_VERSIONS
        """
}
