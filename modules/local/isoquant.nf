process ISOQUANT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::isoquant=3.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.6.1--hdfd78af_0':
        'biocontainers/isoquant:3.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_gtf), path(gtf)
    tuple val(meta_fa), path(fasta)
    tuple val(meta_fai), path(fai)
    val group_category

    output:
    tuple val(meta), path("*.gene*_counts.tsv")         , emit: gene_count_mtx
    tuple val(meta), path("*.transcript*_counts.tsv")   , emit: transcript_count_mtx
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // setting custom home via export (see issue #30)
    if ( !group_category?.trim() ){
        """
        export HOME=\$(pwd)

        isoquant.py ${args} \\
                    --threads $task.cpus \\
                    --datatype nanopore \\
                    --reference $fasta \\
                    --genedb $gtf \\
                    --bam $bam \\
                    -o .
        
        # Move output files from OUT/ to base directory and changing OUT. to their prefix
        find OUT -type f -name 'OUT.*' -exec sh -c 'mv "\$0" "\${0%/*}/${prefix}.\${0##*/OUT.}"' {} \\;
        mv OUT/${prefix}* .
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            isoquant: \$(isoquant.py -v | sed 's#IsoQuant ##')
        END_VERSIONS
        """
    } else {
        """
        export HOME=\$(pwd)

        isoquant.py ${args} \\
                    --threads $task.cpus \\
                    --data_type nanopore \\
                    --reference $fasta \\
                    --genedb $gtf \\
                    --bam $bam \\
                    -o . \\
                    --read_group $group_category

        # Move output files from OUT/ to base directory and changing OUT. to their prefix
        find OUT -type f -name 'OUT.*' -exec sh -c 'mv "\$0" "\${0%/*}/${prefix}.\${0##*/OUT.}"' {} \\;
        mv OUT/${prefix}* .
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            isoquant: \$(isoquant.py -v | sed 's#IsoQuant ##')
        END_VERSIONS
        """

    }
}