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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow PACVAR {
    
    //TODO: from the beginning of the workflow, probably want to flatmap
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    SET_BARCODES_CHANNEL(params.barcodes) // - barcodes fasta
    
    LIMA(ch_samplesheet, SET_BARCODES_CHANNEL.out.data)  // demultiplex 

    PBMM2(LIMA.out.bam, file(params.fasta)) // align using pbmm2 - prob want to couple with some meta?

    //

    //TODO: May want to change to use flatmap from the get, go (to avoid the for loop)
    // PBMM2.out.aligned_bam_ch
    //     .flatMap { tuple ->
    //         def metadata = tuple[0]
    //         def sampleBams = tuple[1]
    //             // Transform each sample BAM into a tuple with the associated metadata
    //         sampleBams.collect { bam -> [metadata, bam] }
    //     }.view()


    //should maybe transform the metadata info to contain the sample barcode as well

    // PBMM2.out.aligned_bam_ch
    //     .flatMap { tuple ->
    //         def metadata = tuple[0]
    //         def sampleBams = tuple[1]
    //             // Transform each sample BAM into a tuple with the associated metadata
    //         sampleBams.collect { bam -> [metadata, bam] }
    // }.flatMap{val -> val}.view()


    comb = PBMM2.out.aligned_bam_ch
        .flatMap { tuple ->
            def metadata = tuple[0]
            def sampleBams = tuple[1]
            // Transform each sample BAM into a tuple with the associated metadata
            sampleBams.collect { bam ->
                // Create a tuple with metadata and BAM path in the desired format
                [metadata, bam]
            }
        }

    PBMM2.out.aligned_bam_ch.view()
    comb.view()



    

    SAMTOOLS_SORT(comb,params.fasta)
    
    // Assuming 'inputChannel' is your original channel containing the tuples




    //PBMM2.out.aligned_bam_ch.flatMap{ list_of_bams -> list_of_bams }.view()

    //bam_paths_channel =  PBMM2.out.aligned_bam_ch.map { tuple -> tuple[1] }.flatMap{ list_of_bams -> list_of_bams }

    //combined_channel = bam_paths_channel.combine(Channel.fromList([params.fasta]))
    //ombined_channel.view()  // Debugging output

    //SAMTOOLS_SORT(combined_channel, meta_ch)
    

    //SAMTOOLS_INDEX(PBMM2.out.aligned_bam_ch)

    //wrap up the bam and the index and the meta - all into one new channel
    // aligned_indexed_ch = PBMM2.out.aligned_bam_ch
    //     .combine(SAMTOOLS_INDEX.out.bai)
    //     .map{(bam_tuple, bai_tuple) ->
    //         def (meta1, bam) = bam_tuple
    //         def (meta2, bai) = bai_tuple
    //         return [meta1, bam, bai]
    //     }

    //not quite sure if converting to a cram is really better 
    //BAM_TO_CRAM(aligned_indexed_ch, params.fasta, params.fasta_fai)

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
