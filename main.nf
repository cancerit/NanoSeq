#!/usr/bin/env nextflow

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

println "run arguments: $params"

/*
 * coverage calculation (cov)
 */
njobs = 12 //fixed
jobIndexes = Channel.of(1..njobs) //fixed

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
    val ii from jobIndexes

    output :
    val ii into cov_out
    file "done" into done_cov //extra protection from error

    cpus 1
    memory '1 GB'

    """
    rm -f $baseDir/tmpNanoSeq/cov/*.done;
    runNanoSeq.py -R $ref -A $normal -B $duplex -j $ii -k $njobs --out $baseDir cov -Q $q --exclude \"$exclude\";
    date > done; 
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
    val cov from cov_out.collect()
    val np from params.jobs

    output :
    file "done" into done_part

    cpus 1
    memory '9 GB'

    """
    rm -f $baseDir/tmpNanoSeq/part/1.done;
    runNanoSeq.py -R $ref -A $normal -B $duplex --out $baseDir part -n $np;
    date > done
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
    path done from done_part
    val np from params.jobs
    val snp_bed from params.snp_bed
    val s_b_index from snp_bed_index
    val noise_bed from params.noise_bed
    val n_b_index from noise_bed_index
    val d from params.dsa_d
    val q from params.dsa_q
    val ii from jobIndexes

    output :
    val ii into dsa_out
    file "done" into done_dsa

    cpus 1
    memory '1 GB'
   
    """
    rm -f $baseDir/tmpNanoSeq/dsa/${ii}.done;
    runNanoSeq.py -R $ref -A $normal -B $duplex --out $baseDir -j $ii -k $np dsa -C $snp_bed -D $noise_bed -d $d -q $q --no_test;
    date > done
    """
}

dsa_out.into{dsa_out1;dsa_out2;dsa_out3;dsa_out4}
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
    val jj from dsa_out2.collect() // wait until all dsa jobs are done
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

    output :
    val ii into var_out
    file "done" into done_var

    cpus 1
    memory '1 GB'

    """
    rm -f $baseDir/tmpNanoSeq/var/${ii}.done;
    runNanoSeq.py -R $ref -A $normal -B $duplex --out $baseDir -j $ii -k $np var -a $a -b $b -c $c -d $d -f $f -i $i \
        -m $m -n $n -p $p -q $q -r $r -v $v -x $x -z $z;
    date > done
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
    val ii from dsa_out3
    val jj from dsa_out4.collect() // wait until all dsa jobs are done
    val np from params.jobs
    val rb from params.indel_rb
    val t3 from params.indel_t3
    val t5 from params.indel_t5
    val mc from params.indel_mc

    output :
    val ii into indel_out
    file "done" into done_indel

    cpus 1
    memory '1 GB'

    """
    rm -f $baseDir/tmpNanoSeq/indel/${ii}.done;
    runNanoSeq.py -R $ref -A $normal -B $duplex --out $baseDir -j $ii -k $np indel --rb $rb --t3 $t3 --t5 $t5 --mc $mc;
    date > done
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
    val ii from var_out.collect() // wait until all var jobs are done
    val jj from indel_out.collect() // wait until all indel jobs are done

    output :
    val "done" into post_out
    file "done" into done_post

    cpus 3
    memory '9 GB'

    """
    rm -f $baseDir/tmpNanoSeq/post/1.done;
    runNanoSeq.py -R $ref -A $normal -B $duplex --out $baseDir -t 3 post;
    date > done
    """
}
