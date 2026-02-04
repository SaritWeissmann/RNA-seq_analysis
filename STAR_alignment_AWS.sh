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

# CONFIGURATION
S3_BUCKET="<bucket_name>"
S3_BASE_PATH="<path_to_fastq_files_inside_bucket/>"                     

GENOME_FASTA="<species.softmasked.fa>"
GENE_ANNOTATION="<species.refGenxyz.gene_exons.gff3>"

THREADS=8
RAM_LIMIT_GB=64                             

# STAR indexing 
SJDB_OVERHANG=99                            
GTF_PARENT_TAG="Parent"                     


# PATHS
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


# Build the STAR index
echo "Checking STAR genome index ..."

if [[ ! -f "${GENOME_INDEX_DIR}/SA" ]]; then
    echo "Building STAR genome index with GTF ..."
    # Use 1024 for binary GiB calculation
    LIMIT_RAM_BYTES=$((RAM_LIMIT_GB * 1024 * 1024 * 1024))

    STAR \
        --runThreadN "${THREADS}" \
        --runMode genomeGenerate \
        --genomeDir "${GENOME_INDEX_DIR}" \
        --genomeFastaFiles "${GENOME_FASTA}" \
        --sjdbGTFfile "${GENE_ANNOTATION}" \
        --sjdbOverhang "${SJDB_OVERHANG}" \
        --limitGenomeGenerateRAM "${LIMIT_RAM_BYTES}" \
        --genomeSAindexNbases 14 \
        --outTmpDir "${GENOME_INDEX_DIR}/_tmp"

    if [[ ! -f "${GENOME_INDEX_DIR}/SA" ]]; then
        echo "ERROR: STAR index build failed"
        exit 1
    fi
else
    echo "STAR index already exists — skipping build"
fi


# Find FASTQ files
echo "Listing FASTQ files in s3://${S3_BUCKET}/${S3_BASE_PATH}"

aws s3 ls --recursive "s3://${S3_BUCKET}/${S3_BASE_PATH}" \
    | awk '/\.(fastq|fq)(\.gz)?$/ {print $4}' \
    | sed "s|^${S3_BASE_PATH}||" \
    | grep -v '^$' > fastq_list.txt

if [[ ! -s fastq_list.txt ]]; then
    echo "ERROR: No FASTQ files found"
    exit 1
fi

echo "Found $(wc -l < fastq_list.txt) FASTQ files"


# Process files
while IFS= read -r s3_key; do
    filename=$(basename "${s3_key}")
    sample="${filename}"
    for ext in ".fastq.gz" ".fq.gz" ".fastq" ".fq"; do
        sample="${sample%"$ext"}"
    done

    if [[ -z "${sample}" || "${sample}" = "${filename}" ]]; then
        continue
    fi

    echo "" 
    echo " Processing: ${sample} (${s3_key})"

    local_fastq="${filename}"
    aws s3 cp "s3://${S3_BUCKET}/${S3_BASE_PATH}${s3_key}" "${local_fastq}" || continue

    outdir="${sample}_star_out"
    mkdir -p "${outdir}"

    # STAR alignment 
    LIMIT_RAM_BYTES=$((RAM_LIMIT_GB * 1024 * 1024 * 1024))

    STAR \
        --runThreadN "${THREADS}" \
        --genomeDir "${GENOME_INDEX_DIR}" \
        --readFilesIn "${local_fastq}" \
        --outFileNamePrefix "${outdir}/" \
        --quantMode TranscriptomeSAM GeneCounts \
        --quantTranscriptomeSAMoutput BanSingleEnd_ExtendSoftclip \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMstrandField intronMotif \
        --sjdbGTFfile "${GENE_ANNOTATION}" \
        --limitBAMsortRAM "${LIMIT_RAM_BYTES}"

    genome_bam="${outdir}/Aligned.sortedByCoord.out.bam"
    if [[ -f "${genome_bam}" ]]; then
        echo "Indexing genome BAM..."
        samtools index "${genome_bam}"
    fi

    # Keep log files
    cp "${outdir}/Log.final.out" "${sample}_Log.final.out" 2>/dev/null || true

    # Upload results
    aws s3 sync "${outdir}/" "s3://${S3_BUCKET}/${RESULTS_S3_PREFIX}${outdir}/" --no-progress

    # Clean up
    rm -f "${local_fastq}"
    rm -rf "${outdir}"

    echo "→ Finished ${sample}"

done < fastq_list.txt

# Final Cleanup
rm -rf "${GENOME_FASTA}" "${GENE_ANNOTATION}" "${GENOME_INDEX_DIR}" fastq_list.txt

echo "All alignments completed!"
