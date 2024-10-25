process FLEXIFORMATTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flexi-formatter%3A1.0.1--pyhdfd78af_0':
        'https://quay.io/biocontainers/flexi-formatter' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_tagged.bam"), emit: bam
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    flexi_formatter \\
        ${bam} \\
        ${args} \\
        | samtools sort -@ $task.cpus -o ${prefix}_tagged.bam -
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexi_formatter: \$(flexi_formatter version |& sed '1!d ; s/flexi_formatter version //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_tagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexi_formatter: \$(flexi_formatter version |& sed '1!d ; s/flexi_formatter version //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
