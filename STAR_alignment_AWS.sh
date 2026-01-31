#!/usr/bin/env bash
# =============================================================================
# Title: Plant RNA-seq STAR Alignment Pipeline (Single-end) - AWS EC2
# Author: Sarit Weissmann
# Latest Update: January 2026
# Description:
#   Plant RNA-seq STAR Alignment Pipeline (Single-end) - AWS EC2:
#    - Downloads sorghum reference + annotation from S3
#    - Builds STAR index with splice junctions from GFF3
#    - Aligns single-end FASTQ files sequentially
#    - Outputs genome + transcriptome BAMs + indexes them +  GeneCounts
#    - Uploads results to S3
#    - Deletes all data from EC2
#    - Designed for r7a.2xlarge (8 vCPU, 64 GiB RAM)
# =============================================================================

set -euo pipefail

# Configuration
S3_BUCKET="<bucket_name>"
S3_BASE_PATH="<path_to_fastq_files_inside_bucket/>"                     

GENOME_FASTA="<species.softmasked.fa>"
GENE_ANNOTATION="<species.refGenxyz.gene_exons.gff3>"

THREADS=8
RAM_LIMIT_GB=64                             

# STAR indexing 
SJDB_OVERHANG=99                            
GTF_PARENT_TAG="Parent"                     

# Paths
WORK_DIR="</home/ubuntu/working_directory_on_EC2>"
GENOME_INDEX_DIR="${WORK_DIR}/star_index"
RESULTS_S3_PREFIX="${S3_BASE_PATH}outputs/"

mkdir -p "${WORK_DIR}" "${GENOME_INDEX_DIR}"
cd "${WORK_DIR}" || { echo "Cannot cd to ${WORK_DIR}"; exit 1; }

# Download reference files from S3
echo "Downloading reference files"

for file in "${GENOME_FASTA}" "${GENE_ANNOTATION}"; do
    if [[ ! -f "${file}" ]]; then
        echo "Downloading ${file} ..."
        aws s3 cp "s3://${S3_BUCKET}/${S3_BASE_PATH}${file}" . || {
            echo "ERROR: Failed to download ${file}"; exit 1;
        }
    else
        echo "Already have: ${file}"
    fi
done

# Build the STAR index for the species
echo "Checking STAR genome index ..."

if [[ ! -f "${GENOME_INDEX_DIR}/SA" ]]; then
    echo "Building STAR genome index ..."

    LIMIT_RAM_BYTES=$((RAM_LIMIT_GB * 1000 * 1000 * 1000))

    STAR \
        --runThreadN "${THREADS}" \
        --runMode genomeGenerate \
        --genomeDir "${GENOME_INDEX_DIR}" \
        --quantMode GeneCounts \
        --genomeFastaFiles "${GENOME_FASTA}" \
        --sjdbGTFfile "${GENE_ANNOTATION}" \
        --sjdbOverhang "${SJDB_OVERHANG}" \
        --sjdbGTFtagExonParentTranscript "${GTF_PARENT_TAG}" \
        --sjdbGTFfeatureExon exon \
        --limitGenomeGenerateRAM "${LIMIT_RAM_BYTES}" \
        --genomeSAindexNbases 14 \
        --genomeChrBinNbits 18 \
        --outTmpDir "${GENOME_INDEX_DIR}/_tmp" 2>&1 | tee star_index.log

    if [[ ! -f "${GENOME_INDEX_DIR}/SA" ]]; then
        echo "ERROR: STAR index build failed — please check your star_index.log"
        exit 1
    fi
else
    echo "STAR index already exists — skipping build"
fi

# Find all fastq files (.gz|.fastq|.fq) in my S3 bucket
echo "Listing FASTQ files in s3://${S3_BUCKET}/${S3_BASE_PATH}"

aws s3 ls --recursive "s3://${S3_BUCKET}/${S3_BASE_PATH}" \
    | awk '/\.(fastq|fq)(\.gz)?$/ {print $4}' \
    | sed "s|^${S3_BASE_PATH}||" \
    | grep -v '^$' > fastq_list.txt

if [[ ! -s fastq_list.txt ]]; then
    echo "ERROR: No .fastq / .fq / .fastq.gz / .fq.gz files found in S3 path"
    exit 1
fi

echo "Found $(wc -l < fastq_list.txt) FASTQ files"

# Process fastq files one by one (Cost restrictions. Can be done in parallel)
while IFS= read -r s3_key; do
    # Clean sample name (remove any known FASTQ extensions)
    filename=$(basename "${s3_key}")
    sample="${filename}"
    for ext in ".fastq.gz" ".fq.gz" ".fastq" ".fq"; do
        sample="${sample%"$ext"}"
    done

    if [[ -z "${sample}" || "${sample}" = "${filename}" ]]; then
        echo "Skipping invalid filename: ${s3_key}"
        continue
    fi

    echo "" 
    echo " Processing: ${sample} (${s3_key})"

    # Keep original extension for local file
    local_fastq="${filename}"

    # Download fastq from S3
    aws s3 cp "s3://${S3_BUCKET}/${S3_BASE_PATH}${s3_key}" "${local_fastq}" || {
        echo "Download failed for ${s3_key} ... skipping"
        continue
    }

    outdir="${sample}_star_out"
    mkdir -p "${outdir}"

    # STAR alignment 
    LIMIT_RAM_BYTES=$((RAM_LIMIT_GB * 1000 * 1000 * 1000))

    STAR \
        --runThreadN "${THREADS}" \
        --genomeDir "${GENOME_INDEX_DIR}" \
        --readFilesIn "${local_fastq}" \
        --outFileNamePrefix "${outdir}/" \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM NM MD \
        --limitBAMsortRAM "${LIMIT_RAM_BYTES}" \
        --sjdbGTFfile "${GENE_ANNOTATION}" \
        2>&1 | tee "${outdir}/star.log"

    # Index BAM files
    for bam in "${outdir}/Aligned.sortedByCoord.out.bam" \
               "${outdir}/Aligned.toTranscriptome.out.bam"; do
        if [[ -f "${bam}" ]]; then
            samtools index "${bam}" || echo "Warning: indexing failed for ${bam}"
        fi
    done

    # Keep log files
    cp "${outdir}/Log.final.out" "${sample}_Log.final.out" 2>/dev/null || true

    # Upload results to my S3 bucket
    aws s3 sync "${outdir}/" "s3://${S3_BUCKET}/${RESULTS_S3_PREFIX}${outdir}/" \
        --no-progress || echo "Warning: upload had issues"

    # Clean fastq and results from EC2 (Again, saves costs)
    rm -f "${local_fastq}"
    rm -rf "${outdir}"

    echo "→ Finished ${sample} — results in s3://${S3_BUCKET}/${RESULTS_S3_PREFIX}${outdir}/"

done < fastq_list.txt

# Delete reference + index to free space
rm -rf "${GENOME_FASTA}" "${GENE_ANNOTATION}" "${GENOME_INDEX_DIR}" fastq_list.txt

# When the analysis is done:
echo "" 
echo "All alignments completed!"
echo "Results are in:"
echo "s3://${S3_BUCKET}/${RESULTS_S3_PREFIX}"

