#!/usr/bin/env bash
# Minimal end-to-end demo of the mapping -> primer trimming -> consensus
# calling -> Freyja variants/demix pipeline used in ../src/run_samples.sh
# and ../src/run_consensus.sh, applied to a handful of raw reads pulled
# directly from NCBI SRA.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

REF="../assets/wg_reference.fasta"
BED="../assets/measles_primers.bed"

mkdir -p data bams trimmed sequences variants depths outputs

while read -r acc; do
    # skip blank lines and comments
    [ -z "$acc" ] && continue
    [[ "$acc" == \#* ]] && continue

    r1="data/${acc}_1.fastq.gz"
    r2="data/${acc}_2.fastq.gz"

    if [ -f "$r1" ] && [ -f "$r2" ]; then
        echo "[${acc}] reads already downloaded"
    else
        echo "[${acc}] downloading from SRA"
        prefetch "$acc" -O data
        fasterq-dump --split-3 --outdir data "data/${acc}/${acc}.sra"
        gzip -f "data/${acc}_1.fastq" "data/${acc}_2.fastq"
        rm -f "data/${acc}.fastq"
        rm -rf "data/${acc}"
    fi

    if [ -f "outputs/${acc}.demixed.tsv" ]; then
        echo "[${acc}] already processed"
        continue
    fi

    echo "[${acc}] mapping"
    minimap2 -ax sr -O16,36 "$REF" "$r1" "$r2" | samtools view -bS - >"bams/${acc}.bam"
    samtools sort -o "bams/${acc}.sorted.bam" "bams/${acc}.bam"
    samtools index "bams/${acc}.sorted.bam"

    echo "[${acc}] trimming primers"
    ivar trim -x 3 -e -m 80 -i "bams/${acc}.sorted.bam" -b "$BED" -p "trimmed/${acc}.bam"
    samtools sort -o "trimmed/${acc}.sorted.bam" "trimmed/${acc}.bam"
    samtools index "trimmed/${acc}.sorted.bam"

    echo "[${acc}] consensus calling"
    samtools mpileup -d 1000 -A -Q 0 "trimmed/${acc}.sorted.bam" | ivar consensus -p "sequences/${acc}" -q 20 -t 0.9

    echo "[${acc}] Freyja variants"
    freyja variants "trimmed/${acc}.sorted.bam" --variants "variants/${acc}.tsv" --depths "depths/${acc}.tsv" --ref "$REF"

    echo "[${acc}] Freyja demix"
    freyja demix "variants/${acc}.tsv" "depths/${acc}.tsv" --output "outputs/${acc}.demixed.tsv" --depthcutoff 5 --autoadapt --pathogen MEASLES
done < accessions.txt

echo "Aggregating demixed outputs"
freyja aggregate outputs/ --output agg_demixed.tsv --ext tsv

echo "Done. See demo/agg_demixed.tsv for combined lineage abundance estimates."
