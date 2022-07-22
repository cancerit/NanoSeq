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
        val include
        val larger

    publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern:"cov/*", overwrite: true

    output :
        tuple val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal) , emit : done
        path 'cov/1.done'
        path 'cov/1.cov.bed.gz'
        path 'cov/gIntervals.dat'
        path 'cov/args.json'
        path 'cov/nfiles'

    cpus 12 // can be modified
    memory { 200.MB  * task.cpus }

    script :
        """
        mkdir -p cov
        touch ${task.process}_${meta.id}
        mkdir -p $baseDir/work/temps/$meta.id/
        rm -rf $baseDir/work/temps/$meta.id/tmpNanoSeq/cov; #allow clean resume for hard links
        
        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex -t $task.cpus --out $baseDir/work/temps/$meta.id/ cov -Q $q --exclude \"$exclude\"  \\
            --include \"$include\" --larger $larger;
        
        mv $baseDir/work/temps/$meta.id/tmpNanoSeq/cov ./cov_tmp
        mv ./cov_tmp/args.json ./cov
        mv ./cov_tmp/gIntervals.dat ./cov
        # this step can produce many small coverage files (one for each contig)
        # consolidate all of them into one file ( better for Lustre performance )
        find ./cov_tmp -name "*.cov.bed.gz" | sort -V | xargs -i cat {} > cov/1.cov.bed.gz
        rm -rf ./cov_tmp
        touch cov/1.done
        echo -n 1 > cov/nfiles
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
        mkdir -p cov
        touch ./cov/1.done
        touch ./cov/1.cov.bed.gz
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
        val excludeCov
        file excludeBED

    publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern: "part/*", overwrite: true

    output :
        tuple  val(meta), path(duplex), path(index_duplex), path(normal), path(index_normal), emit : done
        path 'part/1.done'
        path 'part/intervalsPerCPU.dat'
        path 'part/args.json'


    cpus 1
    memory { task.exitStatus == 130  ? 9.GB * task.attempt : 9.GB }

    script :
        //optional arguments
        def excludeCov_arg = excludeCov != 0 ? "--excludeCov $excludeCov" : ''
        def excludeBED_arg = excludeBED.name != 'NO_FILE_excludeBED' ? "--excludeBED $excludeBED" : ''
        """
        touch ${task.process}_${meta.id}
        mkdir -p part
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/part/*; #allow clean resume for hard links

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex --out $baseDir/work/temps/$meta.id part -n $np $excludeCov_arg $excludeBED_arg;
        
        mv $baseDir/work/temps/$meta.id/tmpNanoSeq/part/* ./part/
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """

    stub :
        """
        mkdir -p part
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
        file snp_bed
        file s_b_index
        file noise_bed
        file n_b_index
        val d
        val q
        each ii

    publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern: "dsa/*", overwrite: true

    output :
        tuple val(metaOut), path(duplex), path(index_duplex), path(normal), path(index_normal), emit : done
        path "dsa/${ii}.done"
        path "dsa/${ii}.dsa.bed.gz"
        path "dsa/nfiles" optional true
        path "dsa/args.json" optional true

    cpus 1
    memory '1.GB'
   
    script :
        //add the job index to the output metadata (to sort processes)
        metaOut = meta.clone()
        metaOut["ii"] = ii
        //these files are optional inputs
        def snp_bed_arg = snp_bed.name != 'NO_FILE_snp_bed' ? "-C $snp_bed" : ''
        def noise_bed_arg = noise_bed.name != 'NO_FILE_noise_bed' ? "-D $noise_bed" : ''
        """
        touch ${task.process}_${meta.id}_${ii}
        mkdir -p dsa
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/${ii}.done; #allow clean resume for hard links
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/${ii}.dsa.bed.gz;
        if [ $ii -eq 1 ]; then rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/nfiles; rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/args.json; fi;

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex --out $baseDir/work/temps/$meta.id -j $ii -k $np dsa -d $d -q $q $snp_bed_arg $noise_bed_arg;
        
        mv $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/${ii}.* ./dsa/ ;
        if [ $ii -eq 1 ]; then mv $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/nfiles ./dsa/; mv $baseDir/work/temps/$meta.id/tmpNanoSeq/dsa/args.json ./dsa/; fi
        cat <<-END_VERSIONS > dsa/versions.yml
        "${task.process}":
            runNanoSeq.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """

    stub :
        metaOut = meta.clone()
        metaOut["ii"] = ii
        """
        mkdir -p dsa
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
        tuple  val(metaOut), path(duplex), path(index_duplex), path(normal), path( index_normal), emit : done
        path "var/${ii}.done"
        path "var/${ii}.var"
        path "var/${ii}.cov.bed.gz"
        path 'var/nfiles' optional true
        path 'var/args.json' optional true

    cpus 1
    memory '800. MB'

    script :
        //add the job index to the output metadata (to sort processes)
        metaOut = meta.clone()
        metaOut["ii"] = ii
        """
        touch ${task.process}_${meta.id}_${ii}
        mkdir -p var
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
        metaOut = meta.clone()
        metaOut["ii"] = ii
        """
        mkdir -p var
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
        val z
        val vaf
        val asxs
        val clip

    publishDir "$baseDir/work/temps/$meta.id/tmpNanoSeq/", mode: 'link', pattern: "indel/*", overwrite: true

    output :
        tuple val(metaOut), path(duplex), path(index_duplex), path(normal), path(index_normal), emit : done
        path "indel/${ii}.done"
        path "indel/${ii}.indel.filtered.vcf.gz"
        path "indel/${ii}.indel.filtered.vcf.gz.tbi"
        path "indel/nfiles" optional true
        path "indel/args.json" optional true

    cpus 1
    memory '2.GB'

    script :
        //apend the index to the output metadata (to sort processes)
        metaOut = meta.clone()
        metaOut["ii"] = ii + np //add np so VAR and INDEL processes sort correctly
        """
        touch ${task.process}_${meta.id}_${ii}
        mkdir -p indel
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/${ii}.done; #allow clean resume for hard links
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/${ii}.indel.filtered.vcf.gz;
        rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/${ii}.indel.filtered.vcf.gz.tbi;
        if [ $ii -eq 1 ]; then rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/nfiles; rm -f $baseDir/work/temps/$meta.id/tmpNanoSeq/indel/args.json; fi;

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex --out $baseDir/work/temps/$meta.id -j $ii -k $np indel --rb $rb --t3 $t3 --t5 $t5 -z $z -v $vaf -a $asxs -c $clip;
        
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
        metaOut = meta.clone()
        metaOut["ii"] = ii + np
        """
        mkdir -p indel
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
    memory {  task.exitStatus == 130  ? 5.GB * task.attempt : 5.GB }

    script :
        def triNuc_arg = triNuc.name != 'NO_FILE_triNuc' ? "--triNuc $triNuc" : ''
        """
        touch ${task.process}_${meta.id}
        mkdir -p post
        rm -f $baseDir/work/temps/${meta.id}/tmpNanoSeq/post/*;

        runNanoSeq.py -R ${ref}/genome.fa -A $normal -B $duplex --out $baseDir/work/temps/${meta.id}/ -t $task.cpus post --name $meta.id $triNuc_arg;
        
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
        mkdir -p post
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

def file_exists(x, name) {
    if ( x == "" ) { return }
    assert file(x).exists() : "\n$name file $x was not found!\n\n"
}

workflow NANOSEQ {
    take :
        bams
        reference
        jobs
        // cov parameters
        cov_q
        cov_exclude
        cov_include
        cov_larger
        // part parameters
        part_excludeBED
        part_excludeCov
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
        indel_z
        indel_v
        indel_a
        indel_c
        // post paramaters
        post_triNuc

    main :
        // optional snp and noise BED files
        if (snp_bed == "" ) {
            snp_bed_fh = file( 'NO_FILE_snp_bed')
            snp_bed_index_fh = file('NO_FILE_snp_bed_index')
        } else {
            snp_bed_fh = file( snp_bed )
            snp_bed_index = snp_bed + ".tbi"
            file_exists(snp_bed_index, "SNP BED index")
            snp_bed_index_fh = file(snp_bed_index)
        }
        if (noise_bed == "") {
            noise_bed_fh = file( 'NO_FILE_noise_bed')
            noise_bed_index_fh = file( 'NO_FILE_noise_bed_index')
        } else {
            noise_bed_fh = file( noise_bed )
            noise_bed_index = noise_bed + ".tbi"
            file_exists(noise_bed_index, "Noise BED index")
            noise_bed_index_fh = file(noise_bed_index)
        }
        // optional file excludeBED
        if (part_excludeBED == "") {
            excludeBED_fh = file( 'NO_FILE_excludeBED')  
        } else {
            excludeBED_fh = file( part_excludeBED )
        }
        // optional triNucleotide file for plots
        if ( post_triNuc == "" ) {
            triNuc_fh = file('NO_FILE_triNuc')
        } else {
            file_exists( post_triNuc, "tri nucleotide ")
            triNuc_fh = file(post_triNuc)
        }
        
        COV(reference, bams, cov_q, cov_exclude, cov_include, cov_larger)
        PART(reference, COV.out.done, jobs, part_excludeCov, excludeBED_fh)
        jobIndexes = Channel.of(1..jobs)
        
        DSA(reference, PART.out.done,
            jobs, snp_bed_fh, snp_bed_index_fh, noise_bed_fh, noise_bed_index_fh, dsa_d, dsa_q, jobIndexes)
        //this allows proper grouping when processing batches of samples so they can proceed with
        // analysis witout holdups and also correctly sorts things so that resume works as expected
        outDSA = DSA.out.done.map{[it[0].id,it]}.groupTuple(size:jobs, sort : { a, b -> a[0]["ii"] <=> b[0]["ii"]}).map{it[1][0].flatten() }
        VAR(reference, outDSA, jobIndexes,
            jobs, var_a, var_b, var_c, var_d, var_f, var_i, var_m, var_n, 
            var_p, var_q, var_r, var_v, var_x, var_z)
        INDEL(reference, outDSA, jobIndexes,
                jobs, indel_rb, indel_t3, indel_t5, indel_z, indel_v, indel_a, indel_c)
        //this allows proper grouping when processing batches of samples so they can proceed with
        // analysis witout holdups and also correctly sorts things so that resume works as expected
        outVAR_INDEL = VAR.out.done.mix(INDEL.out.done).map{[it[0].id,it]}.groupTuple(size: 2*jobs, sort : { a, b -> a[0]["ii"] <=> b[0]["ii"]}).map{it[1][0].flatten()}
        POST(reference, outVAR_INDEL, triNuc_fh )

    emit :
        versions = POST.out.versions
        results = POST.out.results
        tsv = POST.out.tsv
        pdf = POST.out.pdf
        csv = POST.out.csv

}