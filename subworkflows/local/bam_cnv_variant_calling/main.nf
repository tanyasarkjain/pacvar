include { HIFICNV } from '../../../modules/local/hificnv'

workflow BAM_CNV_VARIANT_CALLING {
    take:
    sorted_bam
    sorted_bai
    fasta
    fasta_fai

    main:
    //call the copy number variants
    HIFICNV(sorted_bam,
            sorted_bai,
            fasta)

    //sort and index the resulting bam files
}
