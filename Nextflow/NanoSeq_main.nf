nextflow.enable.dsl=2

params.bwa_container = "/software/CASM/singularity/pcap-core/pcap-core_5.7.0.sif"
params.nanoseq_container= "/lustre/scratch117/casm/team78/ra11/NanoSeq/run_Next/nanoseq_dev.sif"

params.outDir = baseDir
params.ref = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/genome.fa"
file_exists(params.ref,"ref")
reference_path = params.ref.split("/")[0..-2].join('/')

// *** Preprocessing and mapping params
params.fastq_tags_m = 3
params.fastq_tags_s = 4
params.nanoseq_dedup_m = 2

// *** NanoSeq parameters
params.jobs = 100
//  cov
params.cov_q = 0
params.cov_exclude = "MT,GL%,NC_%,hs37d5" //for GRCh37 reference
params.snp_bed = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/botseq/SNP.sorted.bed.gz"
file_exists( params.snp_bed, "snp_bed" )
params.noise_bed = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/botseq/NOISE.sorted.bed.gz"
file_exists( params.noise_bed, "noise_bed" )
params.dsa_d = 2
params.dsa_q = 30
// variantcaller parameters
params.var_a = 50
params.var_b = 5
params.var_c = 0
params.var_d = 2
params.var_f = 0.9
params.var_i = 1.0
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
params.post_triNuc = ""
file_exists(params.post_triNuc,"post_triNuc")

// ** VerifyBAMid params
params.vb_epsilon ="1e-12"
params.vb_ud ="/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/verifybamid/ALL_500K.strictmasked.ok.vcf.UD"
file_exists(params.vb_ud,"vb_ud")
params.vb_bed ="/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/verifybamid/ALL_500K.strictmasked.ok.vcf.bed"
file_exists(params.vb_bed, "vb_bed")
params.vb_mu ="/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/verifybamid/ALL_500K.strictmasked.ok.vcf.mu"
file_exists(params.vb_mu, "vb_mu")

// *** remapping BAMs
params.remap = true
assert ( params.remap == true || params.remap == false ) : "realign parameter must be true or false"

// *** process sample sheet
params.sample_sheet = "sample_sheet.tsv"
file_exists(params.sample_sheet,"sample_sheet")

ss = file( params.sample_sheet )
fields = ss.readLines()[0].split(',')

assert ( fields.contains("id")) : "\n Must specify a unique id column in the sample sheet\n\n"
list_ids = []

// three cases with diffeerent starting files
fastqIn = false
reMapIn = false
noMapIn = false
println(ss)
//fastqs (need alignment)
if ( fields.contains("d_fastq1") && fields.contains("d_fastq2") && fields.contains("n_fastq1") &&
       fields.contains("n_fastq2") ) {
    input_ss = Channel
        .from( ss.text )
        .splitCsv(header: true).map(row -> { 
            assert (  ! list_ids.contains(row.id ) ) : "\nids in sample sheet must be unique\n\n"
            list_ids.add(row.id)
            file_exists(row.d_fastq1,"d_fastq1")
            file_exists(row.d_fastq2,"d_fastq2")
            file_exists(row.n_fastq1,"n_fastq1")
            file_exists(row.n_fastq2,"n_fastq2")
            [  [id: row.id],row.d_fastq1,row.d_fastq2,row.n_fastq1,row.n_fastq2 ]}).multiMap{ it ->
                meta_duplex = it[0].clone()
                meta_duplex["type"]="duplex"
                meta_normal = it[0].clone()
                meta_normal["type"]="normal"
                duplex : [ meta_duplex, [it[1], it[2]] ]
                normal : [ meta_normal, [it[3], it[4]] ]}
    fastqIn = true
//bams needing realignment (no need to check bam indexes)
} else if ( fields.contains("d_bam") && fields.contains("n_bam") && params.remap ) {
    input_ss = Channel
        .from( ss.text )
        .splitCsv(header: true).map(row -> { 
            assert (  ! list_ids.contains(row.id ) ) : "\nids in sample sheet must be unique\n\n"
            list_ids.add(row.id)
            file_exists(row.d_bam,"d_bam")
            file_exists(row.n_bam,"n_bam")
            [  [id: row.id],row.d_bam,row.n_bam ]}).multiMap{ it ->
                meta_duplex = it[0].clone()
                meta_duplex["type"]="duplex"
                meta_normal = it[0].clone()
                meta_normal["type"]="normal"
                duplex : [ meta_duplex, it[1] ]
                normal : [ meta_normal, it[2] ] 
            }
    reMapIn = true
//bams that dont need realignment, (need to check bam indexes)
} else if ( fields.contains("d_bam") && fields.contains("n_bam") && ! params.remap ) {
    input_ss = Channel
        .from( ss.text ).splitCsv(header: true).map(row -> { 
            assert (  ! list_ids.contains(row.id ) ) : "\nids in sample sheet must be unique\n\n"
            list_ids.add(row.id)
            file_exists(row.d_bam,"d_bam")
            file_exists(row.d_bam+".bai", "BAM index")
            file_exists(row.n_bam,"n_bam")
            file_exists(row.n_bam+".bai", "BAM index")
            [  [id: row.id],row.d_bam, row.d_bam + ".bai", row.n_bam, row.n_bam + ".bai"]}).multiMap{ it ->
                meta_duplex = it[0].clone()
                meta_duplex["type"]="duplex"
                meta_normal = it[0].clone()
                meta_normal["type"]="normal"
                duplex : [ meta_duplex, it[1],it[2] ]
                normal : [ meta_normal, it[3],it[4] ]
            }
    noMapIn = true
} else {
    throw new Exception("\nCan't recognize the Sample sheet format\n")
}

def file_exists(x, name = "") {
    if ( x == "" ) { return }
    assert file(x).exists() : "\n$name file $x was not found!\n\n"
}

process FINALIZE {

    input:
        val versions
        val id

    output:
        path "versions.yml"

    publishDir "$params.outDir/outNextflow/$id", mode: 'copy', overwrite: true

    executor 'local'
    cpus 1
    memory 500.MB

    script:
        def allversions = versions.join(' ')
        """
        cat $allversions > versions.yml
        """
    stub:
        def allversions = versions.join(' ')
        """
        cat $allversions > versions.yml
        """
}

println("\n\n")
println("run arguments: $params")
println("\n\n")
println("sample sheet contents:\n")
println(ss.text)
println("\n")
if ( fastqIn ) {
    println("running with fastqs as input that will be trimmed,tagged and mapped\n")
}
if (reMapIn ) {
    println("running with BAMs as input that will be remapped\n")
}
if (noMapIn) {
    println("running with BAMs as input that will not be remapped\n")
}

include { NANOSEQ } from './modules/NanoSeq_analysis.nf'
include { BWAMEM2_MAP; BWAMEM2_REMAP } from './modules/bwa.nf'
include { ADD_NANOSEQ_FASTQ_TAGS; MARKDUP; NANOSEQ_ADD_RB;
         NANOSEQ_DEDUP; VERIFY_BAMID; NANOSEQ_EFFI; NANOSEQ_VAF} from './modules/NanoSeq_aux.nf'


workflow MAP_FASTQ {
    take :
        ch_fastq
        reference_path

    main :
        ADD_NANOSEQ_FASTQ_TAGS( ch_fastq, params.fastq_tags_m, params.fastq_tags_s)
        BWAMEM2_MAP( ADD_NANOSEQ_FASTQ_TAGS.out.fastqs, reference_path )
    
    emit :
        bam =  BWAMEM2_MAP.out.bam
        versions = BWAMEM2_MAP.out.versions
}

versions = Channel.empty() //accumulate all version files here

workflow {
    if ( fastqIn ) { //fastq
 
        //Merge normal's and duplex's into one channel
        ch_fastq_ss = input_ss.duplex.mix( input_ss.normal ).map{ meta, it ->
            meta["name"] = meta.id + "_" + meta.type
            [ meta, it ]}

        MAP = MAP_FASTQ( ch_fastq_ss , reference_path )

    } else { //bam input

        //Merge normal's and duplex's into one channel
        ch_bam_ss = input_ss.duplex.mix( input_ss.normal ).map{ 
            it[0]["name"] = it[0].id + "_" + it[0].type
            it  }

        ch_bam = ch_bam_ss
    }
    
    if ( reMapIn ) { //remapping

        MAP = BWAMEM2_REMAP( ch_bam_ss, reference_path)
   
    }

    if ( fastqIn || reMapIn ) {
        
        MARKDUP( MAP.bam )

        versions = versions.concat(MAP.versions.first())
        versions = versions.concat(MARKDUP.out.versions.first())

        ch_bam = MARKDUP.out.bam
    }

    NANOSEQ_ADD_RB( ch_bam )

    versions = versions.concat(NANOSEQ_ADD_RB.out.versions.first())

    NANOSEQ_DEDUP( NANOSEQ_ADD_RB.out.bam, params.nanoseq_dedup_m )

    versions = versions.concat(NANOSEQ_DEDUP.out.versions.first())

    if ( params.vb_ud != "" &&  params.vb_bed != "" && params.vb_mu != "" ) {
        VERIFY_BAMID( NANOSEQ_DEDUP.out.bam, params.vb_epsilon, reference_path, params.vb_ud, params.vb_bed, params.vb_mu )
    
        versions = versions.concat(VERIFY_BAMID.out.versions.first())
    }

    //collate channels to prepare input for efficiency calculation
    //must provide the BAM from NANOSEQ_ADD_RB and NANOSEQ_DEDUP as arguments
    ch_add_rb_normal =NANOSEQ_ADD_RB.out.bam.filter{ it[0].type == "normal"}.map{[it[0].id, it ]}
    ch_add_rb_duplex =NANOSEQ_ADD_RB.out.bam.filter{ it[0].type == "duplex"}.map{[it[0].id, it ]}

    ch_dedup_normal = NANOSEQ_DEDUP.out.bam.filter{ it[0].type == "normal"}.map{[it[0].id, it ]}
    ch_dedup_duplex = NANOSEQ_DEDUP.out.bam.filter{ it[0].type == "duplex"}.map{[it[0].id, it ]}

    ch_normal_effi = ch_add_rb_normal.join( ch_dedup_normal ).map{it[1] + it[2][1..-1] }
    ch_duplex_effi = ch_add_rb_duplex.join( ch_dedup_duplex ).map{it[1] + it[2][1..-1] }

    ch_input_effi = ch_normal_effi.mix(ch_duplex_effi)
    
    NANOSEQ_EFFI( ch_input_effi, reference_path)
    
    versions = versions.concat(NANOSEQ_EFFI.out.versions.first())
    
    //Collate BAMs for NanoSeq call
    //Must provide :  NANOSEQ_ADD_RB.out.bam (duplex) & NANOSEQ_DEDUP.out.bam (normal) ; as input arguments

    ch_input_nanoseq = ch_add_rb_duplex.join(ch_dedup_normal).map{
        meta = it[1][0].clone()
        meta.name = meta.id
        meta.type = "pair"
        [ meta ] + it[1][1..-1] + it[2][1..-1] }
    
    NANOSEQ( ch_input_nanoseq, reference_path, params.jobs, params.cov_q, params.cov_exclude, 
        params.snp_bed, params.noise_bed, params.dsa_d, params.dsa_q, params.var_a, params.var_b, params.var_c, 
        params.var_d, params.var_f, params.var_i, params.var_m, params.var_n, params.var_p, params.var_q, 
        params.var_r, params.var_v, params.var_x, params.var_z, params.indel_rb, params.indel_t3, params.indel_t5, 
        params.indel_mc, params.post_triNuc)

    versions = versions.concat(NANOSEQ.out.versions.first())

    //NANOSEQ.out.results contains: muts.vcf, muts.vcf.tbi, indel.vcf, indel.vcf.tbi, cov.bed & cov.bed.tbi
    //VAF calculation also requires NANOSEQ_DEDUP.out.bam (duplex) in addition to the previous
    ch_input_vaf = NANOSEQ.out.results.map{[it[0].id, it ]}.join( ch_dedup_duplex ).map{it[1] + it[2][1..-1] }

    NANOSEQ_VAF(ch_input_vaf)

    versions = versions.concat(NANOSEQ_VAF.out.versions.first())

    ch_ids = ch_input_vaf.map{it[0].id}

    FINALIZE( versions.collect(), ch_ids )
}
