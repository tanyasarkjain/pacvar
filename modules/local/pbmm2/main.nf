process PBMM2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.14.99--h9ee0642_0' :
        'biocontainers/pbmm2:1.14.99--h9ee0642_0' }"


    //Takes in a BAM file
    meta = null
    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("*_aligned.bam") , emit: aligned_bam_ch

    //controls when the task is executed -
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    
    prefix="${bam}_aligned"
    pbmm2 align $fasta $bam \${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$( pbmm2 --version | head -n1 | sed 's/pbmm2 //g' | sed 's/ (.\\+//g' )
    END_VERSIONS
    """
}
