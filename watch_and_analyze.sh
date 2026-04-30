#!/bin/bash
# watch_and_analyze.sh
# Live incremental merging + variant calling for a running Nanopore
# sequencing run of a bacterial isolate. Designed to scale to thousands
# of input files by merging only newly-arrived BAMs/BigWigs each cycle.

set -euo pipefail

# --- Config ---
IN_DIR="/mnt/gbi-shared/tmp/labdemo/"
BAM_DIR="${IN_DIR}/bams"
BW_DIR="${IN_DIR}/bigwigs"
OUT_DIR="${IN_DIR}/merged"

MERGED_BAM="${OUT_DIR}/merged.bam"
MERGED_BW="${OUT_DIR}/merged.bw"
VARIANTS="${OUT_DIR}/variants.vcf.gz"

# State files: which inputs have already been folded into the merged output
PROCESSED_BAMS="${OUT_DIR}/.processed_bams.txt"
PROCESSED_BWS="${OUT_DIR}/.processed_bigwigs.txt"

# Sentinel from Nextflow workflow.onComplete
SENTINEL="${IN_DIR}/PIPELINE_DONE"

REFERENCE="${IN_DIR}/reference/genome.fa"
CHROM_SIZES="${IN_DIR}/reference/genome.chrom.sizes"

# Clair3
CLAIR3_MODEL="models/r1041_e82_400bps_sup_v500"
SAMPLE_NAME="isolate"

# Coverage threshold below which we skip variant calling entirely
MIN_COVERAGE_FOR_VARIANTS=5

THREADS=4
MERGE_INTERVAL=60
VARIANT_INTERVAL=120

RUN_CLAIR3="run_clair3.sh"

# --- Setup ---
mkdir -p "$OUT_DIR"
touch "$PROCESSED_BAMS" "$PROCESSED_BWS"

last_variant_run=0
last_nvar="-"
last_variant_status="not yet run"

# --- Helpers ---

# Print a list of files in $1 (glob pattern) that aren't listed in state file $2
new_files() {
    local pattern="$1" state_file="$2"
    shopt -s nullglob
    local all=( $pattern )
    shopt -u nullglob
    [ ${#all[@]} -eq 0 ] && return 0
    # comm requires sorted input
    printf '%s\n' "${all[@]}" | sort > "${state_file}.current"
    sort "$state_file" > "${state_file}.sorted"
    comm -23 "${state_file}.current" "${state_file}.sorted"
    rm -f "${state_file}.current" "${state_file}.sorted"
}

print_status() {
    shopt -s nullglob
    local bams=("${BAM_DIR}"/*.bam)
    local bws=("${BW_DIR}"/*.bw)
    shopt -u nullglob
    local n_bams=${#bams[@]}
    local n_bws=${#bws[@]}
    local n_merged_bams n_merged_bws
    n_merged_bams=$(wc -l < "$PROCESSED_BAMS")
    n_merged_bws=$(wc -l < "$PROCESSED_BWS")

    local n_reads="-" mean_cov="-"
    if [ -f "$MERGED_BAM" ]; then
        n_reads=$(samtools view -c "$MERGED_BAM" 2>/dev/null || echo "-")
        if [ -f "${REFERENCE}.fai" ]; then
            local genome_size total_bases
            genome_size=$(awk '{s+=$2} END {print s}' "${REFERENCE}.fai")
            total_bases=$(samtools depth -a "$MERGED_BAM" 2>/dev/null \
                          | awk '{s+=$3} END {print s+0}')
            if [ "$genome_size" -gt 0 ] 2>/dev/null; then
                mean_cov=$(awk -v t="$total_bases" -v g="$genome_size" \
                           'BEGIN {printf "%.1fx", t/g}')
            fi
        fi
    fi

    local since="never"
    if [ "$last_variant_run" -gt 0 ]; then
        since="$(( ($(date +%s) - last_variant_run) ))s ago"
    fi

    echo "==============================================================="
    echo " [$(date '+%H:%M:%S')] LIVE STATUS"
    echo "  Inputs:    ${n_bams} BAMs   ${n_bws} BigWigs"
    echo "  Merged:    ${n_merged_bams}/${n_bams} BAMs   ${n_merged_bws}/${n_bws} BigWigs"
    echo "  Reads:     ${n_reads}   coverage ${mean_cov}"
    echo "  Variants:  ${last_nvar}   (${last_variant_status}, last attempt: ${since})"
    echo "==============================================================="
}

# Returns mean coverage as a plain number (no 'x' suffix), or 0 if no merged BAM
current_coverage() {
    [ ! -f "$MERGED_BAM" ] && { echo 0; return; }
    [ ! -f "${REFERENCE}.fai" ] && { echo 0; return; }
    local genome_size total_bases
    genome_size=$(awk '{s+=$2} END {print s}' "${REFERENCE}.fai")
    total_bases=$(samtools depth -a "$MERGED_BAM" 2>/dev/null \
                  | awk '{s+=$3} END {print s+0}')
    awk -v t="$total_bases" -v g="$genome_size" \
        'BEGIN {if (g>0) printf "%.1f", t/g; else print 0}'
}

merge_bams() {
    # Find new BAMs not yet merged
    local new_bams
    mapfile -t new_bams < <(new_files "${BAM_DIR}/*.bam" "$PROCESSED_BAMS")
    [ ${#new_bams[@]} -eq 0 ] && return 0

    # Validate each candidate
    local valid_new=()
    for bam in "${new_bams[@]}"; do
        if samtools quickcheck "$bam" 2>/dev/null; then
            valid_new+=("$bam")
        else
            echo "[$(date '+%H:%M:%S')] Skipping invalid BAM (will retry next cycle): $bam"
        fi
    done
    [ ${#valid_new[@]} -eq 0 ] && return 0

    echo "[$(date '+%H:%M:%S')] Adding ${#valid_new[@]} new BAMs to merged output..."

    # Inputs to samtools merge: previous merged file (if any) + new BAMs
    local merge_inputs=()
    [ -f "$MERGED_BAM" ] && merge_inputs+=("$MERGED_BAM")
    merge_inputs+=("${valid_new[@]}")

    if samtools merge -f -@ "$THREADS" "${MERGED_BAM}.tmp" "${merge_inputs[@]}"; then
        samtools index "${MERGED_BAM}.tmp"
        mv "${MERGED_BAM}.tmp"     "${MERGED_BAM}"
        mv "${MERGED_BAM}.tmp.bai" "${MERGED_BAM}.bai"
        # Record which BAMs are now in the merged file
        printf '%s\n' "${valid_new[@]}" >> "$PROCESSED_BAMS"
        echo "[$(date '+%H:%M:%S')] BAM merge done (+${#valid_new[@]} files)."
    else
        echo "[$(date '+%H:%M:%S')] BAM merge failed, keeping previous ${MERGED_BAM}."
        rm -f "${MERGED_BAM}.tmp" "${MERGED_BAM}.tmp.bai"
        return 1
    fi
}

merge_bigwigs() {
    local new_bws
    mapfile -t new_bws < <(new_files "${BW_DIR}/*.bw" "$PROCESSED_BWS")
    [ ${#new_bws[@]} -eq 0 ] && return 0

    echo "[$(date '+%H:%M:%S')] Adding ${#new_bws[@]} new BigWigs to merged output..."

    local merge_inputs=()
    [ -f "$MERGED_BW" ] && merge_inputs+=("$MERGED_BW")
    merge_inputs+=("${new_bws[@]}")

    local tmp_bg="${MERGED_BW}.tmp.bg"
    local tmp_sorted="${MERGED_BW}.tmp.sorted.bg"

    if bigWigMerge "${merge_inputs[@]}" "$tmp_bg" \
       && LC_ALL=C sort -k1,1 -k2,2n "$tmp_bg" > "$tmp_sorted" \
       && bedGraphToBigWig "$tmp_sorted" "$CHROM_SIZES" "${MERGED_BW}.tmp"; then
        mv "${MERGED_BW}.tmp" "$MERGED_BW"
        rm -f "$tmp_bg" "$tmp_sorted"
        printf '%s\n' "${new_bws[@]}" >> "$PROCESSED_BWS"
        echo "[$(date '+%H:%M:%S')] BigWig merge done (+${#new_bws[@]} files)."
    else
        echo "[$(date '+%H:%M:%S')] BigWig merge failed, keeping previous ${MERGED_BW}."
        rm -f "$tmp_bg" "$tmp_sorted" "${MERGED_BW}.tmp"
        return 1
    fi
}

call_variants() {
    if [ ! -f "$MERGED_BAM" ]; then
        last_variant_status="skipped (no merged BAM yet)"
        return 0
    fi

    # Check coverage — skip if too low (early in the run)
    local cov
    cov=$(current_coverage)
    if awk -v c="$cov" -v m="$MIN_COVERAGE_FOR_VARIANTS" \
           'BEGIN {exit !(c < m)}'; then
        echo "[$(date '+%H:%M:%S')] Coverage ${cov}x < ${MIN_COVERAGE_FOR_VARIANTS}x, skipping variant calling."
        last_variant_status="skipped (low coverage: ${cov}x)"
        return 0
    fi

    echo "[$(date '+%H:%M:%S')] Running Clair3 (haploid) at ${cov}x coverage..."

    local tmp_dir="${OUT_DIR}/clair3_tmp"
    rm -rf "$tmp_dir"
    mkdir -p "$tmp_dir"

    if ! "$RUN_CLAIR3" \
            --bam_fn="$MERGED_BAM" \
            --ref_fn="$REFERENCE" \
            --threads="$THREADS" \
            --platform="ont" \
            --model_path="$CLAIR3_MODEL" \
            --output="$tmp_dir" \
            --sample_name="$SAMPLE_NAME" \
            --haploid_sensitive \
            --include_all_ctgs \
            --no_phasing_for_fa \
            --use_gpu \
            > "${tmp_dir}/clair3.log" 2>&1; then
        echo "[$(date '+%H:%M:%S')] Clair3 failed (see ${tmp_dir}/clair3.log) — keeping previous VCF."
        last_variant_status="failed (see ${tmp_dir}/clair3.log)"
        return 1
    fi

    if [ ! -s "${tmp_dir}/merge_output.vcf.gz" ]; then
        echo "[$(date '+%H:%M:%S')] Clair3 produced no output — keeping previous VCF."
        last_variant_status="failed (no output)"
        return 1
    fi

    cp "${tmp_dir}/merge_output.vcf.gz"     "${VARIANTS}.tmp"
    cp "${tmp_dir}/merge_output.vcf.gz.tbi" "${VARIANTS}.tmp.tbi"
    mv "${VARIANTS}.tmp"     "$VARIANTS"
    mv "${VARIANTS}.tmp.tbi" "${VARIANTS}.tbi"

    last_nvar=$(zcat "$VARIANTS" | grep -vc "^#" || true)
    last_variant_status="success"
    echo "[$(date '+%H:%M:%S')] Clair3 done. ${last_nvar} variants in ${VARIANTS}."
}

# --- Main loop ---
echo "[$(date)] Starting watcher."
echo "[$(date)] Watching ${BAM_DIR} and ${BW_DIR}"
echo "[$(date)] Will exit when ${SENTINEL} appears."
echo "[$(date)] Merge every ${MERGE_INTERVAL}s, variants every ${VARIANT_INTERVAL}s (min ${MIN_COVERAGE_FOR_VARIANTS}x)"

while [ ! -f "$SENTINEL" ]; do
    print_status

    merge_bams     || true
    merge_bigwigs  || true

    now=$(date +%s)
    if [ $(( now - last_variant_run )) -ge "$VARIANT_INTERVAL" ] && [ -f "$MERGED_BAM" ]; then
        if call_variants; then
            last_variant_run=$now
        else
            # Don't update last_variant_run on failure — retry next cycle
            :
        fi
    fi

    sleep "$MERGE_INTERVAL"
done

# --- Final pass after sentinel ---
echo "[$(date '+%H:%M:%S')] Sentinel detected. Doing final pass..."
merge_bams     || true
merge_bigwigs  || true
call_variants  || true
print_status
echo "[$(date '+%H:%M:%S')] Done. Exiting."