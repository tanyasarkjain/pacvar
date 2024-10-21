/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pacvar_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SET_VALUE_CHANNEL as SET_BARCODES_CHANNEL } from '../subworkflows/local/set_value_channel'
include { SET_VALUE_CHANNEL } from '../subworkflows/local/set_value_channel'
include { PBMM2 } from '../modules/local/pbmm2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { LIMA } from '../modules/nf-core/lima/main'
include { DEEPVARIANT_RUNDEEPVARIANT } from '../modules/nf-core/deepvariant/rundeepvariant/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { GATK4_HAPLOTYPECALLER } from '../modules/nf-core/gatk4/haplotypecaller/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PACVAR {
    
    //TODO: from the beginning of the workflow, probably want to flatmap
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta // channel: includes the metadata associated with the fasta as well 
    fasta_fai // channel: includes the metadata associateds the fasta index as well
    dict
    dbsnp
    dbsnp_tbi

    main:

    SET_BARCODES_CHANNEL(params.barcodes) 
    LIMA(ch_samplesheet, SET_BARCODES_CHANNEL.out.data)  

    trans_ch = LIMA.out.bam 
        .flatMap{ tuple ->
        def metadata = tuple[0]
        def sampleBams = tuple[1]
        //seperate samples
            sampleBams.collect { bam ->
                [metadata, bam]
        }
        }
        .map{tuple ->
            def bam = tuple[1]
            //change metadata to reflect demultiplexed barcode
            [[id: bam.baseName], bam]
        }.view()
    
    trans_ch.map{tuple ->
    def val1 = tuple[0], 
    de}

    PBMM2(trans_ch, params.fasta)

    SAMTOOLS_SORT(PBMM2.out.aligned_bam_ch, params.fasta)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    intervals = Channel.from([ [], 0 ])
    intervals.view()


    //For now setting intervals to an empty val
    //probably want to add all this logic to a seperate sub-workflow
    deepvar_input_ch = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai).view()
        .map{tuple ->
        def metadata = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        [metadata, bam, bai, []]
     }.view()


    gatk4_input_ch = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai).view()
        .map{tuple ->
        def metadata = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        [metadata, bam, bai, [], []]
     }.view()




    //SINGLE variant calling (may want to move into a different sub-workflow (also should allow somehow for an optional argument))
    //meta data is kinda unneeded 


    dict.view { item -> 
        println "Dict item: ${item}"
    }

    GATK4_HAPLOTYPECALLER(gatk4_input_ch, 
                        fasta,
                        fasta_fai,
                        dict,
                        dbsnp,
                        dbsnp_tbi
                        )

    // Frustrated by the fact that this isn't working :( 

    // DEEPVARIANT_RUNDEEPVARIANT(deepvar_input_ch,
    //                             params.fasta,
    //                             params.fasta_fai)

    //MULTIQC STUFF - NOT QUITE SURE WHAT THIS DOES     
    // MODULE: MultiQC
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions = ch_versions
    LIMA.out.bam
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/