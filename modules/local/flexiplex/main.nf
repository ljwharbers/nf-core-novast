process FLEXIPLEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flexiplex:1.01--py38h2123bcc_1':
        'biocontainers/flexiplex:1.01--py38h2123bcc_1' }"

    input:
    tuple val(meta), path(reads)
    val(adapter_5prime)
    val(adapter_3prime)
    val(spacer_3prime)
    val(ed_adapter_5prime)
    val(ed_adapter_3prime)
    val(barcode_length)
    val(umi_length)

    output:
    tuple val(meta), path("*flexiplex.fastq")        , emit: reads
    tuple val(meta), path("bc_*_reads_barcodes.txt")    , emit: barcodes
    tuple val(meta), path("umi_*_reads_barcodes.txt")   , emit: umis
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}${meta.part ? "_part_${meta.part}" : ''}"
    def barcode = '?' * barcode_length
    def umi = '?' * umi_length
    def half_threads = task.cpus / 2
    """
    # First run flexiplex in discovery move
    flexiplex \\
        -x "${adapter_3prime}" \\
        -b "${barcode}" \\
        -x  "${spacer_3prime}" \\
        -f ${ed_adapter_3prime} \\
        -p ${task.cpus} \\
        -n discovery_${prefix} \\
        ${reads}

    # Run flexiplex with provided discovered barcodes and pipe to get UMIs
    flexiplex \\
        -x "${adapter_3prime}" \\
        -b "${barcode}" \\
        -x  "${spacer_3prime}" \\
        -k discovery_${prefix}_barcodes_counts.txt \\
        -f ${ed_adapter_3prime} \\
        -p ${half_threads} \\
        -n bc_${prefix} \\
        ${reads} | \\
    flexiplex \\
        -x "${adapter_5prime}" \\
        -u "${umi}" \\
        -x "?" \\
        -b "" \\
        -k "?" \\
        -f ${ed_adapter_5prime} \\
        -p ${half_threads} \\
        -n umi_${prefix} \\
        > ${prefix}_flexiplex.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex --help |& sed '1!d ; s/FLEXIPLEX //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.flexiplex.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex --help |& sed '1!d ; s/FLEXIPLEX //')
    END_VERSIONS
    """
}
