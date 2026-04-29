process DORADO_BASECALLER {
    tag "${meta.id}"
    label 'process_gpu'

    container 'nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e'

    input:
    tuple val(meta), path(pod5)
    path  reference
    val   model

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml",           emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    dorado basecaller \\
        -x cuda:all \\
        --reference ${reference} \\
        ${args} \\
        ${model} \\
        ${pod5} \\
        > ${pod5.baseName}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -n1)
        model: ${model}
    END_VERSIONS
    """
}
