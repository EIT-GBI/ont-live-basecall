include { DORADO_BASECALLER } from './modules/local/dorado/basecaller/main.nf'

workflow {
    if (!params.reference) {
        error "Please provide a reference file via --reference <path/to/ref.fa|.mmi>"
    }

    ch_pod = channel.watchPath("${params.folderpath}*.pod5")
        .map { pod -> tuple([id: pod.baseName], pod) }

    DORADO_BASECALLER(
        ch_pod,
        file(params.reference),
        params.dorado_model
    )
}
