#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 00:10:00
#SBATCH -o logs/slurm.submit_sample.%j.out
#SBATCH --mail-type=FAIL

# ==============================================================================
# 2_map_and_demultiplex.sh
#
# Per-sample allele-specific processing of ATAC-seq (paired-end) or RNA-seq
# (single-end) reads. Submitted by 1_prepare_and_submit.sh; each step runs as a
# dependent SLURM job.
#
#                        ATAC-seq                       RNA-seq
#   1. STAR mapping      R1 + R2                        R1
#   2. filtering         MAPQ >= 20, properly paired,   MAPQ >= 20,
#                        ENCODE blacklist removed       ENCODE blacklist removed
#   3. duplicates        removed (Picard)               kept
#   4. SNPsplit          --paired                       single-end
#   genome1 = BALB/cJ, genome2 = CAST/EiJ
#
# Usage:
#   2_map_and_demultiplex.sh <atac|rna> <FASTQ_DIR> <R1.fastq.gz> <STAR_INDEX> <SNPSPLIT_SNP_TXT>
# Output (in <FASTQ_DIR>/mapped/):
#   ATAC: <sample>.Aligned.sortedByCoord.out.mpq20.nodups.genome{1,2}.bam
#   RNA : <sample>.Aligned.sortedByCoord.out.mpq20.genome{1,2}.bam
# ==============================================================================
set -euo pipefail

ASSAY="$1"
FASTQ_DIR="$2"
FIRST_MATE="$3"                       # ATAC: <sample>_R1_<...>.trim.fastq.gz ; RNA: <sample>...fastq.gz
STAR_INDEX="$4"
SNP_TXT="$5"                          # SNPsplit-format CAST SNPs where BALB = reference (from script 1)

# absolute paths, since the SNPsplit job changes directory
FASTQ_DIR=$(realpath -m "$FASTQ_DIR"); FIRST_MATE=$(realpath -m "$FIRST_MATE")
STAR_INDEX=$(realpath -m "$STAR_INDEX"); SNP_TXT=$(realpath -m "$SNP_TXT")

# Tools; the dependent jobs below inherit this environment
module load STAR samtools bedtools2/2.31.0 java
PROJECT=/path/to/async_repli_project
SNPSPLIT=$PROJECT/tools/SNPsplit_v0.3.2
BLACKLIST=$PROJECT/data/mm10-blacklist.v2.bed   # chromosome names must match the BAM (1 vs chr1)
PICARD_JAR=/path/to/picard-tools-2.1.0/picard.jar

case "$ASSAY" in
    atac|rna) ;;
    *) echo "[ERROR] First argument must be 'atac' or 'rna', got '$ASSAY'" >&2; exit 1 ;;
esac
[ -f "$FIRST_MATE" ] || { echo "[ERROR] Read file not found: $FIRST_MATE" >&2; exit 1; }
[ -s "$SNP_TXT" ]    || { echo "[ERROR] SNP file not found: $SNP_TXT" >&2; exit 1; }

if [ "$ASSAY" = atac ]; then
    SECOND_MATE="${FIRST_MATE/_R1_/_R2_}"
    [ -f "$SECOND_MATE" ] || { echo "[ERROR] Mate 2 not found: $SECOND_MATE" >&2; exit 1; }
    READS="$FIRST_MATE $SECOND_MATE"
    PAIR_FILTER="-f 2"                # properly paired reads only
    SNPSPLIT_MODE="--paired"
else
    READS="$FIRST_MATE"
    PAIR_FILTER=""
    SNPSPLIT_MODE=""
fi

SAMPLE=$(basename "$FIRST_MATE")
SAMPLE=${SAMPLE%.fastq.gz}
SAMPLE=${SAMPLE%.trim}
SAMPLE=$(echo "$SAMPLE" | sed 's/_R1_/_/')
MAPPED="$FASTQ_DIR/mapped"
BASE="$MAPPED/$SAMPLE.Aligned.sortedByCoord.out"
mkdir -p "$MAPPED" logs
echo "[SAMPLE] $SAMPLE ($ASSAY): $READS"

# ------------------------------------------------------------------------------
# STEP 1: mapping
# ------------------------------------------------------------------------------
MAP_JOB=$(sbatch --parsable --mem=40G -c 8 -t 20:00:00 -o logs/slurm.STAR.$SAMPLE.%j.out \
  --wrap="STAR --genomeDir $STAR_INDEX \
      --runThreadN 8 \
      --readFilesIn $READS \
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
# STEP 2: MAPQ >= 20 (uniquely mapped) [+ properly paired for ATAC], blacklist removed
# ------------------------------------------------------------------------------
MPQ_JOB=$(sbatch --parsable --dependency=afterok:$MAP_JOB --mem=10G -c 2 -t 04:00:00 \
  -o logs/slurm.mpq.$SAMPLE.%j.out \
  --wrap="samtools view -b -q 20 $PAIR_FILTER $BASE.bam | samtools sort -@ 2 - | \
          bedtools intersect -a - -b $BLACKLIST -v > $BASE.mpq20.bam")
echo "[JOB] MAPQ/blacklist: $MPQ_JOB"

# ------------------------------------------------------------------------------
# STEP 3: duplicate removal (ATAC only)
# ------------------------------------------------------------------------------
if [ "$ASSAY" = atac ]; then
    DUP_JOB=$(sbatch --parsable --dependency=afterok:$MPQ_JOB --mem=20G -c 1 -t 04:00:00 \
      -o logs/slurm.rmDups.$SAMPLE.%j.out \
      --wrap="java -Xmx4g -jar $PICARD_JAR MarkDuplicates \
          INPUT=$BASE.mpq20.bam \
          OUTPUT=$BASE.mpq20.nodups.bam \
          METRICS_FILE=$BASE.dup_metrics.txt \
          REMOVE_DUPLICATES=true")
    echo "[JOB] MarkDuplicates: $DUP_JOB"
    SPLIT_INPUT="$BASE.mpq20.nodups.bam"
    SPLIT_AFTER=$DUP_JOB
else
    SPLIT_INPUT="$BASE.mpq20.bam"
    SPLIT_AFTER=$MPQ_JOB
fi

# ------------------------------------------------------------------------------
# STEP 4: allele assignment
# SNPsplit SNP format: ID  chr  pos  strand  ref/alt
#   ref = C57BL/6 = BALB/cJ base  -> genome1 (BALB/cJ)
#   alt = CAST/EiJ base           -> genome2 (CAST/EiJ)
# For ATAC, --paired requires name-sorted input; SNPsplit sorts the
# coordinate-sorted BAM itself (hence no --no_sort).
# ------------------------------------------------------------------------------
SPLIT_JOB=$(sbatch --parsable --dependency=afterok:$SPLIT_AFTER --mem=25G -c 2 -t 04:00:00 \
  -o logs/slurm.snpSplit.$SAMPLE.%j.out \
  --wrap="cd $MAPPED && $SNPSPLIT/SNPsplit \
      --snp_file $SNP_TXT \
      $SNPSPLIT_MODE \
      --conflicting \
      --weird \
      $SPLIT_INPUT")
echo "[JOB] SNPsplit: $SPLIT_JOB"
