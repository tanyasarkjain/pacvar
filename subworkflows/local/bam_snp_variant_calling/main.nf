include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { DEEPVARIANT_RUNDEEPVARIANT } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'

workflow BAM_SNP_VARIANT_CALLING {
    take:
    sorted_bam
    sorted_bai
    fasta
    fasta_fai
    dict
    dbsnp
    dbsnp_tbi
    intervals

    main:
     //deepvariant
    if (params.tools.split(',').contains('deepvariant')) {
        deepvar_input_ch = sorted_bam.join(sorted_bai)
            .map{tuple ->
            def metadata = tuple[0]
            def bam = tuple[1]
            def bai = tuple[2]
            [metadata, bam, bai, intervals]
        }

        DEEPVARIANT_RUNDEEPVARIANT(deepvar_input_ch,
                                    fasta,
                                    fasta_fai,
                                    [[], []]
                                    )
        vcf_ch = DEEPVARIANT_RUNDEEPVARIANT.out.vcf.join(DEEPVARIANT_RUNDEEPVARIANT.out.vcf_tbi)


    }

    //gatk4_haplotypecaller
    if (params.tools.split(',').contains('gatk4')) {
        gatk4_input_ch = sorted_bam.join(sorted_bai)
            .map{tuple ->
            def metadata = tuple[0]
            def bam = tuple[1]
            def bai = tuple[2]
            [metadata, bam, bai, intervals, []]
        }

        GATK4_HAPLOTYPECALLER(gatk4_input_ch,
                            fasta,
                            fasta_fai,
                            dict,
                            dbsnp,
                            dbsnp_tbi
                            )

        vcf_ch = GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi)
    }

    emit:
    vcf_ch
}