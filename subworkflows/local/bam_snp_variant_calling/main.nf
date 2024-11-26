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
    if (params.snv_caller === 'deepvariant') {
        deepvar_input_ch = sorted_bam.join(sorted_bai)
            .map{metadata, bam, bai ->
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
    if (params.snv_caller === 'gatk4') {
        gatk4_input_ch = sorted_bam.join(sorted_bai)
            .map{metadata, bam, bai ->
                [metadata, bam, bai, intervals, []]
        }

        GATK4_HAPLOTYPECALLER(gatk4_input_ch,
                            fasta,
                            fasta_fai,
                            dict,
                            dbsnp.map{it -> [[:], it]},
                            dbsnp_tbi.map{it -> [[:], it]}
                            )

        vcf_ch = GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi)
    }

    emit:
    vcf_ch
}
