process HIFICNV {
    tag "$meta.id"
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hificnv:1.0.1--h9ee0642_0':
        'biocontainers/hificnv:1.0.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*depth.bw"), emit: depth_bw
    tuple val(meta), path("*.copynum.bedgraph"), emit: bedgraph
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    hificnv \\
        --bam ${bam} \\
        --ref ${fasta} \\
        ${args} \\
        --threads ${task.cpus} \\
        --output-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(hificnv --version |& sed '1!d ; s/hificnv //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.depth.bw
    touch ${prefix}.copynum.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(hificnv --version |& sed '1!d ; s/hificnv //')
    END_VERSIONS
    """
}