#!/usr/bin/env bash
set -euo pipefail

# Per-clone sequencing QC
#
# Requirements:
#   - samtools available on PATH
#   - coordinate-sorted, duplicate-marked BAM files
#
# Expected input layout:
#   bams/
#     sample1.final.bam
#     sample2.final.bam
#
# Usage:
#   chmod +x run_samtools_qc.sh
#   ./run_samtools_qc.sh
#
# Outputs:
#   qc/<sample>.flagstat.txt
#   qc/<sample>.coverage.txt
#   qc/mean_depth_summary.tsv
#
# Table S2 mean depth is calculated as the reference-span-weighted mean of
# the chromosome-level `meandepth` values reported by `samtools coverage`,
# restricted to chr1-chr22, chrX, and chrY:
#
#   sum(meandepth_i * (endpos_i - startpos_i + 1))
#   ------------------------------------------------
#          sum(endpos_i - startpos_i + 1)
#
# Edit BAM_DIR and OUT_DIR if your files are stored elsewhere.

BAM_DIR="bams"
OUT_DIR="qc"
SUMMARY="$OUT_DIR/mean_depth_summary.tsv"

mkdir -p "$OUT_DIR"

shopt -s nullglob
bam_files=("$BAM_DIR"/*.final.bam)

if (( ${#bam_files[@]} == 0 )); then
    echo "ERROR: No *.final.bam files found in $BAM_DIR" >&2
    exit 1
fi

printf 'sample\tmean_depth\n' > "$SUMMARY"

for bam in "${bam_files[@]}"; do
    sample=$(basename "$bam" .final.bam)
    echo "Processing $sample"

    samtools flagstat "$bam" \
        > "$OUT_DIR/${sample}.flagstat.txt"

    samtools coverage "$bam" \
        > "$OUT_DIR/${sample}.coverage.txt"

    mean_depth=$(
        awk 'BEGIN { FS="\t" }
             NR > 1 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {
                 span = $3 - $2 + 1
                 weighted_sum += $7 * span
                 total_span += span
             }
             END {
                 if (total_span == 0) exit 1
                 printf "%.2f", weighted_sum / total_span
             }' "$OUT_DIR/${sample}.coverage.txt"
    ) || {
        echo "ERROR: Could not calculate canonical-chromosome mean depth for $sample" >&2
        exit 1
    }

    printf '%s\t%s\n' "$sample" "$mean_depth" >> "$SUMMARY"
done

echo "QC complete. Results written to $OUT_DIR/"
echo "Weighted mean-depth summary: $SUMMARY"
