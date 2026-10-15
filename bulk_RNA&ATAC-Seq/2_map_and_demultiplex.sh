#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 00:10:00
#SBATCH -o logs/slurm.submit_sample.%j.out
#SBATCH --mail-type=FAIL

# ==============================================================================
# 2_map_and_demultiplex.sh
#
# Per-sample allele-specific ATAC-seq processing (paired-end). Submitted by
# 1_prepare_and_submit.sh; each step runs as a dependent SLURM job.
#
#   1. STAR mapping of R1 + R2 to the N-masked genome (no soft-clipping)
#   2. MAPQ >= 20, properly paired reads, ENCODE blacklist removal
#   3. duplicate removal (Picard MarkDuplicates)
#   4. allele assignment (SNPsplit): genome1 = BALB/cJ, genome2 = CAST/EiJ
#
# Usage: 2_map_and_demultiplex.sh <FASTQ_DIR> <R1.trim.fastq.gz> <STAR_INDEX> <SNPSPLIT_SNP_TXT>
# Output: <FASTQ_DIR>/mapped/<sample>.Aligned.sortedByCoord.out.mpq20.nodups.genome{1,2}.bam
# ==============================================================================
set -euo pipefail

FASTQ_DIR="$1"
FIRST_MATE="$2"                       # <sample>_R1_<...>.trim.fastq.gz
STAR_INDEX="$3"
SNP_TXT="$4"                          # SNPsplit-format CAST SNPs where BALB = reference (from script 1)
SECOND_MATE="${FIRST_MATE/_R1_/_R2_}"

PROJECT=/path/to/async_repli_project
SNPSPLIT=$PROJECT/tools/SNPsplit_v0.3.2
BLACKLIST=$PROJECT/data/mm10-blacklist.v2.bed   # chromosome names must match the BAM (1 vs chr1)
PICARD_JAR=/path/to/picard-tools-2.1.0/picard.jar

[ -s "$SNP_TXT" ] || { echo "[ERROR] SNP file not found: $SNP_TXT" >&2; exit 1; }
[ -f "$SECOND_MATE" ] || { echo "[ERROR] Mate 2 not found: $SECOND_MATE" >&2; exit 1; }

SAMPLE=$(basename "$FIRST_MATE" .trim.fastq.gz | sed 's/_R1_/_/')
MAPPED="$FASTQ_DIR/mapped"
BASE="$MAPPED/$SAMPLE.Aligned.sortedByCoord.out"
mkdir -p "$MAPPED" logs
echo "[SAMPLE] $SAMPLE : $FIRST_MATE + $SECOND_MATE"

# ------------------------------------------------------------------------------
# STEP 1: mapping
# ------------------------------------------------------------------------------
MAP_JOB=$(sbatch --parsable --mem=40G -c 8 -t 20:00:00 -o logs/slurm.STAR.$SAMPLE.%j.out \
  --wrap="STAR --genomeDir $STAR_INDEX \
      --runThreadN 8 \
      --readFilesIn $FIRST_MATE $SECOND_MATE \
      --readFilesCommand zcat \
      --outFileNamePrefix $MAPPED/$SAMPLE. \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMunmapped None \
      --alignEndsType EndToEnd \
      --outFilterScoreMinOverLread 0 \
      --outFilterMatchNminOverLread 0 \
      --outFilterMatchNmin 0 \
      --outFilterMismatchNmax 2 \
      --outSAMattributes NH HI NM MD nM AS")
echo "[JOB] STAR: $MAP_JOB"

# ------------------------------------------------------------------------------
# STEP 2: MAPQ >= 20 (uniquely mapped), properly paired, blacklist removed
# ------------------------------------------------------------------------------
MPQ_JOB=$(sbatch --parsable --dependency=afterok:$MAP_JOB --mem=10G -c 2 -t 04:00:00 \
  -o logs/slurm.mpq.$SAMPLE.%j.out \
  --wrap="samtools view -b -q 20 -f 2 $BASE.bam | samtools sort -@ 2 - | \
          bedtools intersect -a - -b $BLACKLIST -v > $BASE.mpq20.bam")
echo "[JOB] MAPQ/blacklist: $MPQ_JOB"

# ------------------------------------------------------------------------------
# STEP 3: duplicate removal
# ------------------------------------------------------------------------------
DUP_JOB=$(sbatch --parsable --dependency=afterok:$MPQ_JOB --mem=20G -c 1 -t 04:00:00 \
  -o logs/slurm.rmDups.$SAMPLE.%j.out \
  --wrap="java -Xmx4g -jar $PICARD_JAR MarkDuplicates \
      INPUT=$BASE.mpq20.bam \
      OUTPUT=$BASE.mpq20.nodups.bam \
      METRICS_FILE=$BASE.dup_metrics.txt \
      REMOVE_DUPLICATES=true")
echo "[JOB] MarkDuplicates: $DUP_JOB"

# ------------------------------------------------------------------------------
# STEP 4: allele assignment
# SNPsplit SNP format: ID  chr  pos  strand  ref/alt
#   ref = C57BL/6 = BALB/cJ base  -> genome1 (BALB/cJ)
#   alt = CAST/EiJ base           -> genome2 (CAST/EiJ)
# --paired requires name-sorted input; SNPsplit sorts the coordinate-sorted
# BAM itself (hence no --no_sort).
# ------------------------------------------------------------------------------
SPLIT_JOB=$(sbatch --parsable --dependency=afterok:$DUP_JOB --mem=25G -c 2 -t 04:00:00 \
  -o logs/slurm.snpSplit.$SAMPLE.%j.out \
  --wrap="cd $MAPPED && $SNPSPLIT/SNPsplit \
      --snp_file $SNP_TXT \
      --paired \
      --conflicting \
      --weird \
      $BASE.mpq20.nodups.bam")
echo "[JOB] SNPsplit: $SPLIT_JOB"
