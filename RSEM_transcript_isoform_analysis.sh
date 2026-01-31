#!/usr/bin/env bash

# =============================================================================
# Title:       RSEM transcript isoform quantification
# Author:      Sarit Weissmann
# Last Update: January 2026
# Description: Quantifies transcript isoform expression using RSEM after STAR alignment to transcriptome.         
# =============================================================================

set -euo pipefail

# Configuration
GENOME_FASTA="species.softmasked.fa"
GENE_ANNOTATION="species_RefGen.gene_exons.gff3"
THREADS=10
RESULT_PREFIX="RSEM_results"

# Paths
WORK_DIR="/Path to STAR alignment *.*_star_out_Aligned.toTranscriptome.out.bam files"
GENOME_ANNOTATION_DIR="/PAth to annotation files (.fa, .gff3)"
OUT_DIR="/Path to directory you want the results to live in"
RSEM_INDEX_DIR="/Path to the RSEM index"

# Build RSEM index
if [[ ! -f "${RSEM_INDEX_DIR}/${RESULT_PREFIX}.seq" ]]; then
    echo "Building RSEM reference..."

    rsem-prepare-reference \
        --gff3 "${GENOME_ANNOTATION_DIR}/${GENE_ANNOTATION}" \
        -p "${THREADS}" \
        "${GENOME_ANNOTATION_DIR}/${GENOME_FASTA}" \
        "${RSEM_INDEX_DIR}/${RESULT_PREFIX}"
else
    echo "RSEM index already exists, moving on"
fi

# Quantify each sample
shopt -s nullglob

found=false

for file in "${WORK_DIR}"/*_star_out_Aligned.toTranscriptome.out.bam; do
    found=true

    sample=$(basename "$file" _star_out_Aligned.toTranscriptome.out.bam)
    echo "Processing sample: ${sample}"

    rsem-calculate-expression \
        --bam \
        --no-bam-output \
        -p "${THREADS}" \
        --forward-prob 0 \
        "${file}" \
        "${RSEM_INDEX_DIR}/${RESULT_PREFIX}" \
        "${OUT_DIR}/${RESULT_PREFIX}_${sample}"
done

if ! $found; then
    echo "No *_star_out_Aligned.toTranscriptome.out.bam files found in ${WORK_DIR}"
    exit 1
fi

echo "All samples were successfully processed."
exit 0

