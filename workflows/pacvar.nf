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
include { SET_VALUE_CHANNEL as SET_INTERVALS_CHANNEL } from '../subworkflows/local/set_value_channel'
include { BAM_SNP_VARIANT_CALLING as BAM_SNP_VARIANT_CALLING } from '../subworkflows/local/bam_snp_variant_calling'
include { BAM_CNV_VARIANT_CALLING as BAM_CNV_VARIANT_CALLING } from '../subworkflows/local/bam_cnv_variant_calling'
include { BAM_SV_VARIANT_CALLING as BAM_SV_VARIANT_CALLING } from '../subworkflows/local/bam_sv_variant_calling'
include { REPEAT_CHARACTERIZATION as REPEAT_CHARACTERIZATION } from '../subworkflows/local/repeat_characterization'
include { SAMTOOLS_SORT_AND_INDEX as SAMTOOLS_SORT_AND_INDEX } from '../subworkflows/local/samtools_sort_and_index'


include { SET_VALUE_CHANNEL } from '../subworkflows/local/set_value_channel'
include { HIFICNV } from '../modules/local/hificnv'
include { TRGT_GENOTYPE } from '../modules/local/trgt/genotype'
include { TRGT_PLOT } from '../modules/local/trgt/plot'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LIMA } from '../modules/nf-core/lima/main'
include { DEEPVARIANT_RUNDEEPVARIANT } from '../modules/nf-core/deepvariant/rundeepvariant/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { GATK4_HAPLOTYPECALLER } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { PBMM2_ALIGN } from '../modules/nf-core/pbmm2/align/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PACVAR {

    take:
    ch_samplesheet
    fasta
    fasta_fai
    dict
    dbsnp
    dbsnp_tbi

    main:

    SET_BARCODES_CHANNEL(params.barcodes)

    LIMA(ch_samplesheet, SET_BARCODES_CHANNEL.out.data)

    lima_ch = LIMA.out.bam
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

    PBMM2_ALIGN(lima_ch, fasta)

    //ensure the output of this gets paired correctly

    SAMTOOLS_SORT(PBMM2_ALIGN.out.bam, fasta)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    intervals = Channel.from([ [], 0 ])
    intervals.view()

    //join the bam and index based off the meta id (ensure correct order)
    bam_bai_ch = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai).view()
    ordered_bam_ch = bam_bai_ch.map { meta, bam, bai -> [meta, bam] }
    ordered_bai_ch = bam_bai_ch.map { meta, bam, bai -> [meta, bai] }

    ordered_bam_ch.view{ bamfile -> println("BAM: ${bamfile}")}
    ordered_bai_ch.view{ bamfile -> println("BAI: ${bamfile}")}

    //if whole genome sequencing call CNV and SV call the WGS workflow + phase
    if (params.workflow == 'wgs') {
        //gatk or deepvariant snp calling
        BAM_SNP_VARIANT_CALLING(ordered_bam_ch,
                                ordered_bai_ch,
                                fasta,
                                fasta_fai,
                                dict,
                                dbsnp,
                                dbsnp_tbi,
                                params.intervals)

        //pbsv structural variant calling
        BAM_SV_VARIANT_CALLING(ordered_bam_ch,
                                ordered_bai_ch,
                                fasta,
                                fasta_fai)

        // // //hificnv
        // BAM_CNV_VARIANT_CALLING(ordered_bam_ch,
        //                         ordered_bai_ch,
        //                         fasta,
        //                         fasta_fai)
    }

    if (params.workflow == 'repeat') {

        intervals_ch = Channel.fromPath(params.intervals).map { file ->[file.baseName, file] }
        id_ch = Channel.fromPath(params.id).map { file ->[file.baseName, file] }

        REPEAT_CHARACTERIZATION(ordered_bam_ch,
                                    ordered_bai_ch,
                                    fasta,
                                    fasta_fai,
                                    intervals_ch,
                                    id_ch)
    }

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