#!/usr/bin/env nextflow

params.outDir = baseDir
params.duplex = "$baseDir/33526_2#4.bam"
params.normal = "$baseDir/33526_2#10.neat.bam"
params.ref = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/genome.fa"
params.jobs = 80

// indexes
d_i_ext = "." + params.duplex.split("\\.")[-1][0..1] + "i"
duplex_index = params.duplex + d_i_ext
n_i_ext = "." + params.normal.split("\\.")[-1][0..1] + "i"
normal_index = params.normal + n_i_ext
ref_index = params.ref + ".fai"

// cov parameters
params.cov_q = 0
params.cov_exclude = "MT,GL%,NC_%,hs37d5" //for GRCh37 reference

// dsa parameters
params.snp_bed = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/botseq/SNP.sorted.bed.gz"
params.noise_bed = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/botseq/NOISE.sorted.bed.gz"
snp_bed_index = params.snp_bed + ".tbi"
noise_bed_index = params.noise_bed + ".tbi"
params.dsa_d = 2
params.dsa_q = 30

// variantcaller parameters
params.var_a = 2 
params.var_b = 5
params.var_c = 0
params.var_d = 2
params.var_f = 0.9
params.var_i = 1
params.var_m = 8
params.var_n = 3
params.var_p = 0
params.var_q = 60
params.var_r = 144
params.var_v = 0.01
params.var_x = 8
params.var_z = 12

// indel parameters
params.indel_rb = 2
params.indel_t3 = 135
params.indel_t5 = 10
params.indel_mc = 16

// post paramaters
params.post_name = "results"
params.post_triNuc = ""
if ( params.post_triNuc == "" ) {
    triNuc_fh = file('NO_FILE')
} else {
    triNuc_fh = file(params.post_triNuc, checkIfExists: true)
}

println "run arguments: $params"

/*
 * coverage calculation (cov)
 */

cov_cpus = 12 //fixed

process cov {
    input :
    path ref from params.ref
    path ri from ref_index
    path duplex from params.duplex
    path di from duplex_index
    path normal from params.normal
    path ni from normal_index
    val q from params.cov_q
    val exclude from params.cov_exclude

    publishDir "$params.outDir/tmpNanoSeq/cov", mode: 'link', overwrite: true

    output :
    path '*.done' into cov_files_done
    path '*.cov.bed.gz' into cov_files_bed
    path 'gIntervals.dat' into cov_files_dat
    path 'args.json' into cov_files_json
    path 'nfiles' into cov_files_n

    cpus cov_cpus 
    memory { 500.MB * cov_cpus}

    """
    rm -f $params.outDir/tmpNanoSeq/cov/*; #allow clean resume for hard links

    runNanoSeq.py -R $ref -A $normal -B $duplex -t $cov_cpus --out $params.outDir cov -Q $q --exclude \"$exclude\";
    
    mv $params.outDir/tmpNanoSeq/cov/* .
    """
}

/*
 * partition calculation (part)
 */
process part {
    input :
    path ref from params.ref
    path ri from ref_index
    path duplex from params.duplex
    path di from duplex_index
    path normal from params.normal
    path ni from normal_index
    path ifiles from cov_files_done.collect()
    val np from params.jobs

    publishDir "$params.outDir/tmpNanoSeq/part", mode: 'link', overwrite: true

    output :
    path 'intervalsPerCPU.dat' into part_files_dat
    path 'args.json' into part_files_json
    path '1.done' into part_files_done

    cpus 1
    memory '9 GB'

    """
    rm -f $params.outDir/tmpNanoSeq/part/*; #allow clean resume for hard links

    runNanoSeq.py -R $ref -A $normal -B $duplex --out $params.outDir part -n $np;
    
    mv $params.outDir/tmpNanoSeq/part/* .
    """
}

/*
 * bed table calculation (dsa)
 */
jobIndexes = Channel.of(1..params.jobs)
process dsa {
    input :
    path ref from params.ref
    path ri from ref_index
    path duplex from params.duplex
    path di from duplex_index
    path normal from params.normal
    path ni from normal_index
    path ifile from part_files_done
    val np from params.jobs
    val snp_bed from params.snp_bed
    val s_b_index from snp_bed_index
    val noise_bed from params.noise_bed
    val n_b_index from noise_bed_index
    val d from params.dsa_d
    val q from params.dsa_q
    val ii from jobIndexes

    publishDir "$params.outDir/tmpNanoSeq/dsa", mode: 'link', overwrite: true

    output :
    val ii into dsa_out
    path "${ii}.done" into dsa_files_done
    path "${ii}.dsa.bed.gz" into dsa_files_bed
    path "nfiles" optional true into dsa_files_n
    path "args.json" optional true into dsa_files_json

    cpus 1
    memory '1 GB'
   
    """
    rm -f $params.outDir/tmpNanoSeq/dsa/${ii}.done; #allow clean resume for hard links
    rm -f $params.outDir/tmpNanoSeq/dsa/${ii}.dsa.bed.gz;
    if [ $ii -eq 1 ]; then rm -f $params.outDir/tmpNanoSeq/dsa/nfiles; rm -f $params.outDir/tmpNanoSeq/dsa/args.json; fi;

    runNanoSeq.py -R $ref -A $normal -B $duplex --out $params.outDir -j $ii -k $np dsa -C $snp_bed -D $noise_bed -d $d -q $q --no_test;
    
    mv $params.outDir/tmpNanoSeq/dsa/${ii}.* . ;
    if [ $ii -eq 1 ]; then mv $params.outDir/tmpNanoSeq/dsa/nfiles .; mv $params.outDir/tmpNanoSeq/dsa/args.json .; fi
    """
}

dsa_out.into{dsa_out1;dsa_out2}
dsa_files_done.into{dsa_files_done1;dsa_files_done2}
/*
 * variant calculation (var)
 */
process var {
    input :
    path ref from params.ref
    path ri from ref_index
    path duplex from params.duplex
    path di from duplex_index
    path normal from params.normal
    path ni from normal_index
    val ii from dsa_out1
    path ifiles from dsa_files_done1.collect() // wait until all dsa jobs are done
    val np from params.jobs
    val a from params.var_a
    val b from params.var_b
    val c from params.var_c
    val d from params.var_d
    val f from params.var_f
    val i from params.var_i
    val m from params.var_m
    val n from params.var_n
    val p from params.var_p
    val q from params.var_q
    val r from params.var_r
    val v from params.var_v
    val x from params.var_x
    val z from params.var_z

    publishDir "$params.outDir/tmpNanoSeq/var", mode: 'link', overwrite: true

    output :
    path "${ii}.done" into var_files_done
    path "${ii}.var" into var_files_var
    path "${ii}.cov.bed.gz" into var_files_bed
    path "nfiles" optional true into var_files_n
    path "args.json" optional true into var_files_json

    cpus 1
    memory '1 GB'

    """
    rm -f $params.outDir/tmpNanoSeq/var/${ii}.done; #allow clean resume for hard links
    rm -f $params.outDir/tmpNanoSeq/var/${ii}.var;
    rm -f $params.outDir/tmpNanoSeq/var/${ii}.cov.bed.gz;
    if [ $ii -eq 1 ]; then rm -f $params.outDir/tmpNanoSeq/var/nfiles; rm -f $params.outDir/tmpNanoSeq/var/args.json; fi;

    runNanoSeq.py -R $ref -A $normal -B $duplex --out $params.outDir -j $ii -k $np var -a $a -b $b -c $c -d $d -f $f -i $i \
        -m $m -n $n -p $p -q $q -r $r -v $v -x $x -z $z;
    
    mv $params.outDir/tmpNanoSeq/var/${ii}.* . ;
    if [ $ii -eq 1 ]; then mv $params.outDir/tmpNanoSeq/var/nfiles .; mv $params.outDir/tmpNanoSeq/var/args.json .; fi
    """
}

/*
 * indel calculation (indel)
 */
process indel {
    input :
    path ref from params.ref
    path ri from ref_index
    path duplex from params.duplex
    path di from duplex_index
    path normal from params.normal
    path ni from normal_index
    val ii from dsa_out2
    path ifiles from dsa_files_done2.collect() // wait until all dsa jobs are done
    val np from params.jobs
    val rb from params.indel_rb
    val t3 from params.indel_t3
    val t5 from params.indel_t5
    val mc from params.indel_mc

    publishDir "$params.outDir/tmpNanoSeq/indel", mode: 'link', overwrite: true

    output :
    path "${ii}.done" into indel_files_done
    path "${ii}.indel.filtered.vcf.gz" into indel_files_vcf
    path "${ii}.indel.filtered.vcf.gz.tbi" into indel_files_tbi
    path "nfiles" optional true into indel_files_n
    path "args.json" optional true into indel_files_json

    cpus 1
    memory '1 GB'

    """
    rm -f $params.outDir/tmpNanoSeq/indel/${ii}.done; #allow clean resume for hard links
    rm -f $params.outDir/tmpNanoSeq/indel/${ii}.indel.filtered.vcf.gz;
    rm -f $params.outDir/tmpNanoSeq/indel/${ii}.indel.filtered.vcf.gz.tbi;
    if [ $ii -eq 1 ]; then rm -f $params.outDir/tmpNanoSeq/indel/nfiles; rm -f $params.outDir/tmpNanoSeq/indel/args.json; fi;

    runNanoSeq.py -R $ref -A $normal -B $duplex --out $params.outDir -j $ii -k $np indel --rb $rb --t3 $t3 --t5 $t5 --mc $mc;
    
    rm $params.outDir/tmpNanoSeq/indel/${ii}.indel.bed.gz; #not required for final calculations
    rm $params.outDir/tmpNanoSeq/indel/${ii}.indel.vcf.gz; #not required for final calculations
    mv $params.outDir/tmpNanoSeq/indel/${ii}.* . ;
    if [ $ii -eq 1 ]; then mv $params.outDir/tmpNanoSeq/indel/nfiles .; mv $params.outDir/tmpNanoSeq/indel/args.json .; fi
    """
}

/*
 * post calculation (post)
 */
process post {
    input :
    path ref from params.ref
    path ri from ref_index
    path duplex from params.duplex
    path di from duplex_index
    path normal from params.normal
    path ni from normal_index
    val name from params.post_name
    val ii from var_files_done.collect() // wait until all var jobs are done
    val jj from indel_files_done.collect() // wait until all indel jobs are done
    file triNuc from triNuc_fh

    publishDir "$params.outDir/tmpNanoSeq/post", mode: 'link', overwrite: true

    output :
    val true into post_out
    path '*' into post_files_out

    cpus 3
    memory '9 GB'

    script:
    def triNuc_arg = triNuc.name != 'NO_FILE' ? "--triNuc $triNuc" : ''

    """
    rm -f $params.outDir/tmpNanoSeq/post/*;

    runNanoSeq.py -R $ref -A $normal -B $duplex --out $params.outDir -t 3 post --name $name $triNuc_arg;
    
    mv $params.outDir/tmpNanoSeq/post/* . ;
    
    """
}
