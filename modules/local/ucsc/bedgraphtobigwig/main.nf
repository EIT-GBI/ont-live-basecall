process UCSC_BEDGRAPHTOBIGWIG {
    tag "${meta.id}"

    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h2a80c09_0'

    input:
    tuple val(meta), path(bedgraph)
    path  chrom_sizes

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path  "versions.yml",          emit: versions

    script:
    """
    bedGraphToBigWig \\
        ${bedgraph} \\
        ${chrom_sizes} \\
        ${bedgraph.baseName}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc-bedgraphtobigwig: 445
    END_VERSIONS
    """
}
