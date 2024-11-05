include { PBSV_DISCOVER } from '../../../modules/nf-core/pbsv/discover/main'
include { PBSV_CALL } from '../../../modules/nf-core/pbsv/call/main'

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
}
