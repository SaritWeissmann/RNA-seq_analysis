#!/bin/bash
set -euo pipefail
set -x   # remove later once it works

# Input reads
read1=$1
read2=$2

# Extract sample base name (robust)
base=$(basename "$read1" | sed 's/_R[12].*//')

echo "Sample base: $base"

# Path to software
BOWTIE2=/home/usr/bowtie2-2.4.4/bowtie2
INDEX=usr/Genomes/Indices/BT_index/Sbicolor_PI_156178/Sbicolor_PI_156178
SAMTOOLS=/home/usr/samtools-1.7/samtools

# Alignment + BAM generation pipeline
$BOWTIE2 \
  --threads 15 \
  -x "$INDEX" \
  -1 "$read1" \
  -2 "$read2" \
| $SAMTOOLS view -@ 15 -b - \
| $SAMTOOLS sort -@ 15 -m 4G -T "${base}.tmp" -o "${base}.srt.bam"

# Index BAM
$SAMTOOLS index "${base}.srt.bam"

echo "Done: ${base}.srt.bam"
