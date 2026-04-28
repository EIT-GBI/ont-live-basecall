include { DORADO_BASECALLER     } from './modules/local/dorado/basecaller/main.nf'
include { SAMTOOLS_SORT_INDEX   } from './modules/local/samtools/sort_index/main.nf'
include { GETCHROMSIZES         } from './modules/local/getchromsizes/main.nf'
include { BEDTOOLS_GENOMECOV    } from './modules/local/bedtools/genomecov/main.nf'
include { UCSC_BEDGRAPHTOBIGWIG } from './modules/local/ucsc/bedgraphtobigwig/main.nf'

workflow {
    if (!params.reference) {
        error "Please provide a reference file via --reference <path/to/ref.fa|.mmi>"
    }

    ch_reference = file(params.reference)

    if (params.chrom_sizes) {
        ch_chrom_sizes = file(params.chrom_sizes)
    } else {
        GETCHROMSIZES(ch_reference)
        ch_chrom_sizes = GETCHROMSIZES.out.sizes.first()
    }

// ch_pod = channel.watchPath("${params.folderpath}*.{pod5,done}")
//    .until { it.name.endsWith('.done') }   // or 'final_summary.txt', etc.

    ch_pod = channel.watchPath("${params.folderpath}*.pod5")
        .map { pod -> tuple([id: pod.baseName], pod) }

    DORADO_BASECALLER(
        ch_pod,
        ch_reference,
        params.dorado_model
    )

    SAMTOOLS_SORT_INDEX(DORADO_BASECALLER.out.bam)
    BEDTOOLS_GENOMECOV(SAMTOOLS_SORT_INDEX.out.bam)
    UCSC_BEDGRAPHTOBIGWIG(BEDTOOLS_GENOMECOV.out.bedgraph, ch_chrom_sizes)
}
