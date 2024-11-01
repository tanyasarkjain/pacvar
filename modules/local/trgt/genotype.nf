
process TRGT_GENOTYPE {
    tag "$meta.id"
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trgt:1.2.0--h9ee0642_0':
        'biocontainers/trgt:1.2.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bed)

    output:
    tuple val(vcf), path("*.vcf.gz"), emit: vcf
    tuple val("*.spanning.bam"), emit: spanning_bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    trgt \\
        genotype \\
        $args \\
        --genome ${fasta} \\
        --repeats ${bed} \\
        --reads ${bam} \\
        --output-prefix ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trgt: \$(trgt --version |& sed '1!d ; s/trgt //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.spanning.bam
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trgt: \$(trgt --version |& sed '1!d ; s/trgt //')
    END_VERSIONS
    """
}
