#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 0-06:00
#SBATCH -o logs/slurm.bin_counts.%j.out
#SBATCH --mail-type=FAIL,END

# ==============================================================================
# 3_allelic_bin_counts.sh
#
# Counts allele-specific ATAC reads (SNPsplit output of script 2) in 100 kb
# genomic bins annotated by asynchronous-replication (ASRT) status, producing
# the count table read by 4_plot_allelic_ratio.R.
#
# Usage:
#   sbatch scripts/3_allelic_bin_counts.sh <FASTQ_DIR> <SAMPLE_REP1> <SAMPLE_REP2>
#
#   FASTQ_DIR    directory passed to script 2 (contains mapped/)
#   SAMPLE_REPn  sample name of each cl3 replicate as printed by script 2 ([SAMPLE] line)
#
# Output:
#   data/coverage_genomes.cl3_atac.bed  (tab-separated, no header)
#     col 1-4 : chr, start, end, annotation (All / C57BL / CAST)
#     col 5-8 : read counts rep1_genome1, rep1_genome2, rep2_genome1, rep2_genome2
#   SNPsplit genome1 = BALB/cJ, genome2 = CAST/EiJ
# ==============================================================================
set -euo pipefail

FASTQ_DIR="$1"
REP1="$2"
REP2="$3"

DATA="data"
BINS="$DATA/coverage_genomes.Britnys_fig.mm10.bed"              # 100 kb bins
ASRT_LIST="$DATA/async_master_list.GRCm38.bed"                  # col 4 = early allele: CAST, or C57BL (= BALB/c, reference-like allele)
ANNOT="$DATA/coverage_genomes.Britnys_fig.mm10.ASRT_annotation.bed"
OUT="$DATA/coverage_genomes.cl3_atac.bed"

module load samtools
module load bedtools2/2.31.0

# ------------------------------------------------------------------------------
# STEP 1: annotate each bin by ASRT status
#   bins not overlapping an ASRT region -> "All"
#   bins overlapping an ASRT region     -> early-replicating allele from col 4
#   (-wao reports bins without overlap with start = -1 in the ASRT columns)
# The output keeps the bin order of $BINS; the count table and the R script
# rely on this order.
# ------------------------------------------------------------------------------
echo "[STEP 1] Annotating bins with ASRT status"
for f in "$BINS" "$ASRT_LIST"; do
    [ -s "$f" ] || { echo "[ERROR] Missing input: $f" >&2; exit 1; }
done
bedtools intersect -a "$BINS" -b "$ASRT_LIST" -wao | \
    awk 'BEGIN{OFS="\t"} $5=="-1" {print $1,$2,$3,"All"} $5!="-1" {print $1,$2,$3,$7}' \
    > "$ANNOT"
# a bin overlapping >1 ASRT region would be duplicated; abort if so
dups=$(cut -f1-3 "$ANNOT" | sort | uniq -d | wc -l)
if [ "$dups" -ne 0 ]; then
    echo "[ERROR] $dups bins overlap more than one ASRT region" >&2; exit 1
fi
echo "[STEP 1] $(wc -l < "$ANNOT") bins: $(cut -f4 "$ANNOT" | sort | uniq -c | tr '\n' ' ')"

# ------------------------------------------------------------------------------
# STEP 2: coordinate-sort and index the SNPsplit per-allele BAMs
# ------------------------------------------------------------------------------
BAMS=()
for s in "$REP1" "$REP2"; do
    for g in genome1 genome2; do
        in=$(ls "$FASTQ_DIR"/mapped/"$s"*mpq20.nodups."$g".bam 2>/dev/null | head -n1 || true)
        if [ -z "$in" ]; then
            echo "[ERROR] No SNPsplit $g BAM found for sample $s" >&2; exit 1
        fi
        sorted="${in%.bam}.sorted.bam"
        echo "[STEP 2] Sorting $in"
        samtools sort -@ 4 -o "$sorted" "$in"
        samtools index "$sorted"
        BAMS+=("$sorted")
    done
done

# ------------------------------------------------------------------------------
# STEP 3: count reads per bin (bin order of $ANNOT is preserved; R relies on it)
# BAM order: rep1_genome1 rep1_genome2 rep2_genome1 rep2_genome2
# ------------------------------------------------------------------------------
echo "[STEP 3] Counting reads per bin"
bedtools multicov -bams "${BAMS[@]}" -bed "$ANNOT" > "$OUT"

echo "[SUCCESS] Wrote $OUT ($(wc -l < "$OUT") bins)"
