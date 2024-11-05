include { TRGT_GENOTYPE } from '../../../modules/local/trgt/genotype'
include { TRGT_PLOT } from '../../../modules/local/trgt/plot'
include { BCFTOOLS_SORT } from '../../../modules/nf-core/bcftools/sort/main'
include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow  REPEAT_CHARACTERIZATION{

    take:
    sorted_bam
    sorted_bai
    fasta
    fasta_fai
    bed
    id

    main:
    //genotype the repeat region
    TRGT_GENOTYPE(sorted_bam,
                    sorted_bai,
                    fasta,
                    fasta_fai,
                    bed)

    //sort the resulting spanning bam
    SAMTOOLS_SORT(TRGT_GENOTYPE.out.spanning_bam,
                    fasta)

    //index the resulting bam
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    //sort the resulting vcf
    BCFTOOLS_SORT(TRGT_GENOTYPE.out.vcf)

    //plot the vcf file -- for a specified id
    TRGT_PLOT(SAMTOOLS_SORT.out.bam,
                SAMTOOLS_INDEX.out.bai,
                BCFTOOLS_SORT.out.vcf,
                fasta,
                fasta_fai,
                bed,
                id)
}