process GETCHROMSIZES {
    tag "${fasta.simpleName}"

    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    path fasta

    output:
    path "*.sizes",      emit: sizes
    path "*.fai",        emit: fai
    path fasta,          emit: fa
    path "versions.yml", emit: versions

    script:
    """
    samtools faidx ${fasta}
    cut -f1,2 ${fasta}.fai > ${fasta.baseName}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
