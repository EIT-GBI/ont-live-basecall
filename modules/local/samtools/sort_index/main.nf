process SAMTOOLS_SORT_INDEX {
    tag "${meta.id}"

    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"),     emit: bam
    tuple val(meta), path("*.sorted.bam.bai"), emit: bai
    path  "versions.yml",                      emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools sort \\
        ${args} \\
        -@ ${task.cpus} \\
        -o ${bam.baseName}.sorted.bam \\
        ${bam}

    samtools index \\
        -@ ${task.cpus} \\
        ${bam.baseName}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
