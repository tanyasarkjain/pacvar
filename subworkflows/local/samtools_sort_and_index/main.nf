include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'

workflow SAMTOOLS_SORT_AND_INDEX {
    take:
    aligned_bam
    fasta

    main:

    // Perform sorting
    SAMTOOLS_SORT(aligned_bam, fasta)

    // Perform indexing on the sorted BAM
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    // Emit the outputs
    emit:
    sorted_bam = SAMTOOLS_SORT.out.bam
    sorted_bai = SAMTOOLS_INDEX.out.bai
}
