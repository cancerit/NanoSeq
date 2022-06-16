nextflow.enable.dsl=2

/*
 * coverage calculation (cov)
 */
process COV {

    tag "$meta.id"

    container params.nanoseq_container

    input :
        path ref
        tuple val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal)
        val q
        val exclude

    publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern:"cov/*", overwrite: true

    output :
        tuple val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal) , emit : done
        path 'cov/*.done'
        path 'cov/*.cov.bed.gz'
        path 'cov/gIntervals.dat'
        path 'cov/args.json'
        path 'cov/nfiles'

    cpus 12 // can be modified
    memory { 200.MB  * task.cpus }

    script :
        """
        touch ${task.process}_${meta.id}
        mkdir cov
        mkdir -p $baseDir/work/temps/$meta.id/
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/cov/*; #allow clean resume for hard links
        

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex -t $task.cpus --out $baseDir/work/temps/$meta.id/ cov -Q $q --exclude \"$exclude\";
        
        mv $baseDir/work/temps/$meta.id/tmpNanoSeq/cov/* ./cov
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """

    stub:
        """
        mkdir -p $baseDir/work/temps/$meta.id/tmpNanoSeq/cov
        mkdir -p $baseDir/work/temps/$meta.id/tmpNanoSeq/part
        mkdir -p $baseDir/work/temps/$meta.id/tmpNanoSeq/var
        mkdir -p $baseDir/work/temps/$meta.id/tmpNanoSeq/indel
        mkdir -p $baseDir/work/temps/$meta.id/tmpNanoSeq/post
        mkdir cov
        touch ./cov/{1..$task.cpus}.done
        touch ./cov/{1..$task.cpus}.cov.bed.gz
        touch ./cov/gIntervals.dat
        touch ./cov/args.json
        touch ./cov/nfiles
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """
}

/*
 * partition calculation (part)
 */
process PART {

    tag "$meta.id"

    container params.nanoseq_container

    input :
        path ref
        tuple  val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal)
        val np

    publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern: "part/*", overwrite: true

    output :
        tuple  val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal), emit : done
        path 'part/1.done'
        path 'part/intervalsPerCPU.dat'
        path 'part/args.json'


    cpus 1
    memory '9 GB'

    script :
        """
        touch ${task.process}_${meta.id}
        mkdir part
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/part/*; #allow clean resume for hard links

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex --out $baseDir/work/temps/$meta.id part -n $np;
        
        mv $baseDir/work/temps/$meta.id/tmpNanoSeq/part/* ./part/
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """

    stub :
        """
        mkdir part
        touch part/1.done
        touch part/intervalsPerCPU.dat
        touch part/args.json
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """
}

/*
 * bed table calculation (dsa)
 */

process DSA {

    tag "${meta.id}_${ii}"

    container params.nanoseq_container

    input :
        path ref
        tuple  val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal)
        val np
        val snp_bed
        val s_b_index
        val noise_bed
        val n_b_index
        val d
        val q
        each ii

    publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern: "dsa/*", overwrite: true

    output :
        tuple val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal), emit : done
        path "dsa/${ii}.done"
        path "dsa/${ii}.dsa.bed.gz"
        path "dsa/nfiles" optional true
        path "dsa/args.json" optional true

    cpus 1
    memory '800. MB'
   
    script :
        """
        touch ${task.process}_${meta.id}_${ii}
        mkdir dsa
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/${ii}.done; #allow clean resume for hard links
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/${ii}.dsa.bed.gz;
        if [ $ii -eq 1 ]; then rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/nfiles; rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/args.json; fi;

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex --out $baseDir/work/temps/$meta.id -j $ii -k $np dsa -C $snp_bed -D $noise_bed -d $d -q $q --no_test;
        
        mv $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/${ii}.* ./dsa/ ;
        if [ $ii -eq 1 ]; then mv $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/nfiles ./dsa/; mv $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/args.json ./dsa/; fi
        cat <<-END_VERSIONS > dsa/versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """

    stub :
        """
        mkdir dsa
        touch ./dsa/${ii}.done
        touch ./dsa/${ii}.dsa.bed.gz
        if [ $ii -eq 1 ]; then touch ./dsa/nfiles; touch ./dsa/args.json; fi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """ 
}


/*
 * variant calculation (var)
 */
process VAR {

    tag "${meta.id}_${ii}"
    
    container params.nanoseq_container

    input :
        path ref
        tuple val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal)
        each ii
        val np
        val a
        val b 
        val c
        val d
        val f
        val i
        val m
        val n
        val p
        val q
        val r
        val v
        val x
        val z

    publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern: "var/*", overwrite: true

    output :
        tuple  val(meta), path(duplex), path(index_duplex), path(normal), path( index_normal), emit : done
        path "var/${ii}.done"
        path "var/${ii}.var"
        path "var/${ii}.cov.bed.gz"
        path 'var/nfiles' optional true
        path 'var/args.json' optional true

    cpus 1
    memory '800. MB'

    script :
        """
        touch ${task.process}_${meta.id}_${ii}
        mkdir var
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/var/${ii}.done; #allow clean resume for hard links
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/var/${ii}.var;
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/var/${ii}.cov.bed.gz;
        if [ $ii -eq 1 ]; then rm -f $baseDir/work//temps/$meta.id/tmpNanoSeq/var/nfiles; rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/var/args.json; fi;

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex --out $baseDir/work/temps/$meta.id/ -j $ii -k $np var -a $a -b $b -c $c -d $d -f $f -i $i \
            -m $m -n $n -p $p -q $q -r $r -v $v -x $x -z $z;
        
        mv $baseDir/work/temps/$meta.id/tmpNanoSeq/var/${ii}.* ./var/ ;
        if [ $ii -eq 1 ]; then mv $baseDir/work/temps/$meta.id/tmpNanoSeq/var/nfiles ./var/; mv $baseDir/work/temps/$meta.id/tmpNanoSeq/var/args.json ./var/; fi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """

    stub :
        """
        mkdir var
        touch ./var/${ii}.done
        touch ./var/${ii}.var
        touch ./var/${ii}.cov.bed.gz
        if [ $ii -eq 1 ]; then touch ./var/nfiles; touch ./var/args.json; fi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """ 
}

/*
 * indel calculation (indel)
 */
process INDEL {

    tag "${meta.id}_${ii}"

    container params.nanoseq_container

    input :
        path ref
        tuple  val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal)
        each ii
        val np
        val rb
        val t3
        val t5
        val mc

    publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern: "indel/*", overwrite: true

    output :
        tuple val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal), emit : done
        path "indel/${ii}.done"
        path "indel/${ii}.indel.filtered.vcf.gz"
        path "indel/${ii}.indel.filtered.vcf.gz.tbi"
        path "indel/nfiles" optional true
        path "indel/args.json" optional true

    cpus 1
    memory '2. GB'

    script :
        """
        touch ${task.process}_${meta.id}_${ii}
        mkdir indel
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/${ii}.done; #allow clean resume for hard links
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/${ii}.indel.filtered.vcf.gz;
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/${ii}.indel.filtered.vcf.gz.tbi;
        if [ $ii -eq 1 ]; then rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/nfiles; rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/args.json; fi;

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex --out $baseDir/work/temps/$meta.id -j $ii -k $np indel --rb $rb --t3 $t3 --t5 $t5 --mc $mc;
        
        rm $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/${ii}.indel.bed.gz; #not required for final calculations
        rm $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/${ii}.indel.vcf.gz; #not required for final calculations
        mv $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/${ii}.* ./indel/ ;
        if [ $ii -eq 1 ]; then mv $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/nfiles ./indel/; mv $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/args.json ./indel/; fi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """

    stub :
        """
        mkdir indel
        touch indel/${ii}.done
        touch indel/${ii}.indel.filtered.vcf.gz
        touch indel/${ii}.indel.filtered.vcf.gz.tbi
        if [ $ii -eq 1 ]; then touch indel/nfiles; touch indel/args.json; fi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """ 
}

/*
 * post calculation (post)
 */
process POST {

    tag "${meta.id}"

    container params.nanoseq_container

    input :
        path ref
        tuple  val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal)
        file triNuc

    publishDir "$params.outDir/outNextflow/$meta.id", mode: 'link', pattern : "post/*", overwrite: true

    output :
        path("versions.yml"), emit: versions
        tuple val(meta), path("post/${meta.id}.muts.vcf.gz"), path( "post/${meta.id}.muts.vcf.gz.tbi"), 
            path("post/${meta.id}.indel.vcf.gz"), path( "post/${meta.id}.indel.vcf.gz.tbi"), path("post/${meta.id}.cov.bed.gz"), 
            path( "post/${meta.id}.cov.bed.gz.tbi"), emit: results
        path("post/*.csv"), emit: csv
        path("post/*.tsv"), emit: tsv optional true
        path("post/*.pdf"), emit: pdf optional true


    cpus 2
    memory '1.GB'

    script :
        def triNuc_arg = triNuc.name != 'NO_FILE' ? "--triNuc $triNuc" : ''
        """
        touch ${task.process}_${meta.id}
        mkdir post
        rm -f $baseDir/work/temps/${meta.id}/tmpNanoSeq/post/*;

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex --out $baseDir/work/temps/${meta.id}/ -t 2 post --name $meta.id $triNuc_arg;
        
        mv $baseDir/work/temps/${meta.id}/tmpNanoSeq/post/* ./post/ ;
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
            bcftools: \$(echo \$(bcftools --version 2>&1) | head -1 | sed 's/^.*bcftools //; s/Using.*\$//')
            R : \$(R --version | head -1 | sed -r 's/.*version ([[0-9.]*) .*/\1/')
        END_VERSIONS
        """
    
    stub :
        """
        mkdir post
        touch post/${meta.id}.cov.bed.gz
        touch post/${meta.id}.cov.bed.gz.tbi
        touch post/${meta.id}.muts.vcf.gz
        touch post/${meta.id}.muts.vcf.gz.tbi
        touch post/${meta.id}.indel.vcf.gz
        touch post/${meta.id}.indel.vcf.gz.tbi
        touch post/file.tsv
        touch post/file.pdf
        touch post/file.csv
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
            bcftools: \$(echo \$(bcftools --version 2>&1) | head -1 | sed 's/^.*bcftools //; s/Using.*\$//')
            R : \$(R --version | head -1 | sed -r 's/.*version ([[0-9.]*) .*/\\1/')
        END_VERSIONS
        """
}

workflow NANOSEQ {
    take :
        bams
        reference
        jobs
        // cov parameters
        cov_q
        cov_exclude
        // dsa parameters
        snp_bed
        noise_bed
        dsa_d
        dsa_q
        // variantcaller parameters
        var_a
        var_b
        var_c
        var_d
        var_f
        var_i
        var_m
        var_n
        var_p
        var_q
        var_r
        var_v
        var_x
        var_z
        // indel parameters
        indel_rb
        indel_t3
        indel_t5
        indel_mc
        // post paramaters
        post_triNuc

    main :
        // indexes
        snp_bed_index = snp_bed + ".tbi"
        noise_bed_index = noise_bed + ".tbi"
        // optional triNucleotide file for plots
        if ( post_triNuc == "" ) {
            triNuc_fh = file('NO_FILE')
        } else {
            triNuc_fh = file(post_triNuc, checkIfExists: true)
        }
        
        COV(reference, bams, cov_q, cov_exclude)
        PART(reference, COV.out.done, jobs)
        jobIndexes = Channel.of(1..jobs)
        
        DSA(reference, PART.out.done,
            jobs, snp_bed, snp_bed_index, noise_bed, noise_bed_index, dsa_d, dsa_q, jobIndexes)
        //this allows proper grouping when processing batches of samples so they can proceed with analysis witout holdups
        outDSA = DSA.out.done.map{[it[0].id,it]}.groupTuple(size:jobs).map{it[1][0].flatten()}
        VAR(reference, outDSA, jobIndexes,
            jobs, var_a, var_b, var_c, var_d, var_f, var_i, var_m, var_n, 
            var_p, var_q, var_r, var_v, var_x, var_z)
        INDEL(reference, outDSA, jobIndexes,
                jobs, indel_rb, indel_t3, indel_t5, indel_mc)
        //this allows proper grouping when processing batches of samples so they can proceed with analysis witout holdups
        outVAR_INDEL = VAR.out.done.mix(INDEL.out.done).map{[it[0].id,it]}.groupTuple(size: 2*jobs).map{it[1][0].flatten()}
        POST(reference, outVAR_INDEL, triNuc_fh )

    emit :
        versions = POST.out.versions
        results = POST.out.results
        tsv = POST.out.tsv
        pdf = POST.out.pdf
        csv = POST.out.csv

}