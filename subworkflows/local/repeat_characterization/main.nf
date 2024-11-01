// include { TRGT_GENOTYPE } from '../../../modules/local/trgt/genotype'
// include { TRGT_PLOT } from '../../../modules/local/trgt/plot'
// include { BCFTOOLS_SORT } from '../../../modules/local/bcftools/sort'
// include { SAMTOOLS_SORT } from '../../../modules/local/samtools/sort'

workflow  REPEAT_CHARACTERIZATION{

    take:
    sorted_bam
    sorted_bai
    fasta
    fasta_fai
    repeat_region

    //genotype the repeat region
    TRGT_GENOTYPE(sorted_bam,
                    fasta,
                    bed)

    //sort the resulting spanning ba
    SAMTOOLS_SORT(TRGT_GENOTYPE.out.spanning_bam,
                    fasta)

    //sort the resulting spanning ba
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam,
                    fasta)

    //sort the resulting vcf
    SAMTOOLS_INDEX(TRGT_GENOTYPE.out.vcf,
                    fasta)

    //plot the vcf file
    TRGT_GENOTYPE(sorted_bam,
                    fasta,
                    bed)

}