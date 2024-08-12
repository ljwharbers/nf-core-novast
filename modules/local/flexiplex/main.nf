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
    val(ed_adapter_5prime)
    val(ed_adapter_3prime)
    val(barcode_length)
    val(umi_length)

    output:
    tuple val(meta), path("*flexiplex.fastq.gz")    , emit: reads
    tuple val(meta), path("combined_*barcodes.txt") , emit: barcodes
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapter_3prime_1 = adapter_3prime.substring(0, adapter_3prime.length() - 4)
    def adapter_3prime_2 = adapter_3prime.substring(adapter_3prime.length() - 4)
    def barcode = '?' * barcode_length
    def umi = '?' * umi_length
    def threads = task.cpus / 2
    """
    gunzip -c $reads \\
        | flexiplex \\
        -x $adapter_3prime_1 \\
        -b $adapter_3prime_2 \\
        -u $barcode \\
        -f $ed_adapter_3prime \\
        -e 1 \\
        -k $adapter_3prime_2 \\
        -p $threads \\
        $args \\
        -n firstpass_${prefix} \\
        | flexiplex \\
        -x $adapter_5prime \\
        -u $umi \\
        -x "?" \\
        -b "" \\
        -k "?" \\
        -f $ed_adapter_5prime \\
        -n combined_${prefix} \\
        -p $threads \\
        $args2 \\
        | gzip > ${prefix}.flexiplex.fastq.gz

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
