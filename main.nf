#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/pacvar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/pacvar
    Website: https://nf-co.re/pacvar
    Slack  : https://nfcore.slack.com/channels/pacvar
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { getGenomeAttribute	  } from './subworkflows/local/utils_nfcore_pacvar_pipeline'

params.fasta = getGenomeAttribute('fasta')
params.fasta_fai = getGenomeAttribute('fasta_fai')
params.dbsnp = getGenomeAttribute('dbsnp')
params.dbsnp_tbi = getGenomeAttribute('dbsnp_tbi')
params.dict = getGenomeAttribute('dict')

include { PACVAR  } from './workflows/pacvar'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_pacvar_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_pacvar_pipeline'

// Initialize genomic attibutes with associated meta data maps
fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
fasta_fai = params.fasta_fai ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
dict = params.dict ? Channel.fromPath(params.dict).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
dbsnp = params.dbsnp ? Channel.fromPath(params.dbsnp).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
dbsnp_tbi = params.dbsnp_tbi ? Channel.fromPath(params.dbsnp_tbi).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

intervals = params.intervals ? Channel.fromPath(params.intervals).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
id = params.id ? Channel.fromPath(params.id).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PACVAR {


    take:
    samplesheet // channel: samplesheet read in from --input
    fasta
    fasta_fai
    dict
    dbsnp
    dbsnp_tbi
    intervals
    id


    main:

    // println "pacvar Fasta parameter: ${params.fasta}"
    // println "pacvar genome name: ${params.genome}"


    //
    // WORKFLOW: Run pipeline
    //
    PACVAR (
        samplesheet,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi,
        intervals,
        id
    )


    emit:
    multiqc_report = PACVAR.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow {

    main:

    println "main Fasta parameter: ${params.fasta}"
    println "main genome name: ${params.genome}"

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PACVAR (
        PIPELINE_INITIALISATION.out.samplesheet,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi,
        intervals,
        id
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_PACVAR.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
