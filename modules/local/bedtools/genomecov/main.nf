process BEDTOOLS_GENOMECOV {
    tag "${meta.id}"

    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bedgraph"), emit: bedgraph
    path  "versions.yml",                emit: versions

    script:
    def args = task.ext.args ?: '-bga'
    """
    bedtools genomecov \\
        ${args} \\
        -ibam ${bam} \\
        > ${bam.baseName}.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
