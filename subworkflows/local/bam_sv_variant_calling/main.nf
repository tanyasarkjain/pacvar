include { PBSV_DISCOVER } from '../../../modules/nf-core/pbsv/discover/main'
include { PBSV_CALL } from '../../../modules/nf-core/pbsv/call/main'
include { BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index/main'
include { TABIX_BGZIP } from '../../../modules/nf-core/tabix/bgzip/main'

workflow BAM_SV_VARIANT_CALLING {
    take:
    sorted_bam
    sorted_bai
    fasta
    fasta_fai

    main:
    //call the structural variants
    PBSV_DISCOVER(sorted_bam, fasta)
    PBSV_CALL(PBSV_DISCOVER.out.svsig, fasta)

    //zip and index
    TABIX_BGZIP(PBSV_CALL.out.vcf)
    BCFTOOLS_INDEX(TABIX_BGZIP.out.output)

    vcf_ch = TABIX_BGZIP.out.output.join(BCFTOOLS_INDEX.out.csi)

    emit:
    vcf_ch
}
