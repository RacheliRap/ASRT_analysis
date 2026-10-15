#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH -N 1
#SBATCH -o logs/slurm.prepare.%j.out
#SBATCH -t 2-00:00
#SBATCH --mail-type=FAIL,END

# ==============================================================================
# 1_prepare_and_submit.sh
#
# Prepares the SNP files, the N-masked reference genome and the STAR index for
# allele-specific ATAC-seq and RNA-seq in CAST/EiJ x BALB/cJ hybrid cells, trims
# the ATAC reads, then submits 2_map_and_demultiplex.sh for every sample.
#   ATAC-seq: paired-end, trimmed
#   RNA-seq : single-end, not trimmed
#
# Allele convention used throughout:
#   genome1 = BALB/cJ   (reference-like allele, C57BL/6 base at the SNPs used)
#   genome2 = CAST/EiJ
#
# Usage: edit the paths below, uncomment the steps to run, then from the
# repository root (the logs/ folder must exist before submitting):
#   mkdir -p logs
#   sbatch scripts/1_prepare_and_submit.sh
# ==============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# Paths (edit for your system)
# ------------------------------------------------------------------------------
PROJECT=/path/to/async_repli_project
DATA=$PROJECT/data                                   # SNP files, masked genome, STAR index
ATAC_FASTQ_DIR=$DATA/ATAC-seq/sample_fastq           # paired-end *_R1_*.fastq.gz / *_R2_*.fastq.gz
RNA_FASTQ_DIR=$DATA/RNA-seq/sample_fastq             # single-end *.fastq.gz
SNPSPLIT=$PROJECT/tools/SNPsplit_v0.3.2
REF_FASTA_DIR=/path/to/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta  # per-chromosome FASTA (GRCm38)
GTF=/path/to/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf

# Mouse Genomes Project v5 per-strain SNP files (GRCm38, dbSNP142)
CAST_VCF=$DATA/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf
BALB_VCF=$DATA/BALB_cJ.mgp.v5.snps.dbSNP142.vcf

MASK_VCF=$DATA/BALB_cJ.CAST_EiJ.SNPs_for_mask.vcf        # all positions to N-mask
SPLIT_VCF=$DATA/CAST_EiJ.no_BALB_cJ.mgp.v5.snps.dbSNP142.vcf  # SNPs used by SNPsplit
SPLIT_TXT=$DATA/CAST_EiJ.no_BALB_cJ.snpsplit.txt              # same, in SNPsplit format
MASKED_DIR=$DATA/CAST_EiJ_N-masked
STAR_INDEX=$DATA/BALB_cJ.masked_CAST.STAR_index

module load STAR
module load bedtools2/2.31.0
mkdir -p logs

# ------------------------------------------------------------------------------
# STEP 1: SNP files
#   - Masking file: every position where CAST or BALB differs from the reference,
#     so that neither allele is penalised during mapping.
#   - SNPsplit file: CAST SNPs at positions where BALB carries the reference base.
#     At these positions reference base = BALB (genome1), alternative = CAST
#     (genome2). BALB-specific and shared positions are masked but are not used
#     for allele assignment, since they cannot distinguish the two alleles as
#     reference vs alternative.
# ------------------------------------------------------------------------------
function create_snps() {
    echo "[STEP 1] Building masking and SNPsplit SNP files"
    grep '^#' "$CAST_VCF" > "$DATA/CAST_EiJ.mgp.v5.snps.dbSNP142.header"

    # CAST SNPs where BALB = reference  -> masked, and used by SNPsplit
    bedtools intersect -v -header -a "$CAST_VCF" -b "$BALB_VCF" > "$SPLIT_VCF"
    # BALB SNPs where CAST = reference  -> masked only
    bedtools intersect -v -a "$BALB_VCF" -b "$CAST_VCF" > "$DATA/BALB_cJ.no_CAST.mgp.v5.snps.dbSNP142.vcf"
    # BALB SNPs at positions also variable in CAST -> masked only
    bedtools intersect -u -a "$BALB_VCF" -b "$CAST_VCF" > "$DATA/BALB_cJ.overlap_CAST.mgp.v5.snps.dbSNP142.vcf"

    # masking file = union of the three sets, sorted, with the VCF header
    { grep -v '^#' "$SPLIT_VCF"
      cat "$DATA/BALB_cJ.no_CAST.mgp.v5.snps.dbSNP142.vcf" \
          "$DATA/BALB_cJ.overlap_CAST.mgp.v5.snps.dbSNP142.vcf"; } | \
        bedtools sort | \
        cat "$DATA/CAST_EiJ.mgp.v5.snps.dbSNP142.header" - > "$MASK_VCF"
    # SNPsplit format: ID chr pos strand ref/alt  (ref -> genome1 BALB, alt -> genome2 CAST)
    grep -v '^#' "$SPLIT_VCF" | awk 'BEGIN{OFS="\t"} {print NR, $1, $2, 1, $4"/"$5}' > "$SPLIT_TXT"
    echo "[STEP 1] Done"
}

# ------------------------------------------------------------------------------
# STEP 2: N-masked reference genome (all positions in $MASK_VCF set to N)
# ------------------------------------------------------------------------------
function genome_preparation() {
    echo "[STEP 2] Building N-masked genome"
    rm -rf "$MASKED_DIR"
    cd "$DATA"
    "$SNPSPLIT"/SNPsplit_genome_preparation \
        --nmasking \
        --strain CAST_EiJ \
        --reference_genome "$REF_FASTA_DIR" \
        --vcf_file "$MASK_VCF"
    cd - > /dev/null
    ls "$MASKED_DIR"/chr*.fa > /dev/null 2>&1 || { echo "[ERROR] No N-masked FASTA files in $MASKED_DIR" >&2; exit 1; }
    cat $(ls "$MASKED_DIR"/chr*.fa | sort -V) > "$MASKED_DIR"/genome.fa
    echo "[STEP 2] Done"
}

# ------------------------------------------------------------------------------
# STEP 3: STAR index of the N-masked genome
# (the same index is used for the RNA-seq data, hence the annotation GTF)
# ------------------------------------------------------------------------------
function genome_index() {
    echo "[STEP 3] Building STAR index"
    mkdir -p "$STAR_INDEX"
    STAR --runThreadN 8 \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$MASKED_DIR"/genome.fa \
        --sjdbGTFfile "$GTF"
    echo "[STEP 3] Done"
}

# ------------------------------------------------------------------------------
# STEP 4: ATAC read trimming (both mates). RNA-seq reads are not trimmed.
# ------------------------------------------------------------------------------
function trim_fastq() {
    echo "[STEP 4] Trimming ATAC reads"
    command -v fastx_trimmer > /dev/null || { echo "[ERROR] fastx_trimmer (FASTX-Toolkit) not on PATH" >&2; exit 1; }
    local files=( "$ATAC_FASTQ_DIR"/*_R[12]_*.fastq.gz )
    [ -e "${files[0]}" ] || { echo "[ERROR] No ATAC FASTQ files in $ATAC_FASTQ_DIR" >&2; exit 1; }
    for f in "${files[@]}"; do
        [[ "$f" == *trim* ]] && continue
        out="${f%.fastq.gz}.trim.fastq.gz"
        zcat "$f" | fastx_trimmer -Q33 -f 16 -l 103 -z -o "$out"
    done
    echo "[STEP 4] Done"
}

# ------------------------------------------------------------------------------
# STEP 5: submit mapping + allele assignment for every sample
#   the assay (atac / rna) is passed explicitly to script 2
# ------------------------------------------------------------------------------
function submit_atac() {
    echo "[STEP 5] Submitting ATAC-seq samples (paired-end)"
    local files=( "$ATAC_FASTQ_DIR"/*_R1_*.trim.fastq.gz )
    [ -e "${files[0]}" ] || { echo "[ERROR] No trimmed ATAC R1 files in $ATAC_FASTQ_DIR (run trim_fastq first)" >&2; exit 1; }
    for r1 in "${files[@]}"; do
        echo "[SUBMIT] atac $r1"
        sbatch scripts/2_map_and_demultiplex.sh atac "$ATAC_FASTQ_DIR" "$r1" "$STAR_INDEX" "$SPLIT_TXT"
    done
}

function submit_rna() {
    echo "[STEP 5] Submitting RNA-seq samples (single-end)"
    local files=( "$RNA_FASTQ_DIR"/*.fastq.gz )
    [ -e "${files[0]}" ] || { echo "[ERROR] No RNA FASTQ files in $RNA_FASTQ_DIR" >&2; exit 1; }
    for r1 in "${files[@]}"; do
        echo "[SUBMIT] rna $r1"
        sbatch scripts/2_map_and_demultiplex.sh rna "$RNA_FASTQ_DIR" "$r1" "$STAR_INDEX" "$SPLIT_TXT"
    done
}

# ==============================================================================
# MAIN — uncomment the steps to run (steps 1-3 are run once)
# ==============================================================================
# create_snps
# genome_preparation
# genome_index
# trim_fastq
submit_atac
submit_rna
