#!/bin/bash

#SBATCH -o logs/slurm.sample_node.%j.out
#SBATCH --mail-type=FAIL

# Capture operational variables passed from Orchestrator script
CLONE_DIR="$1"
FIRST_MATE="$2"
STAR_INDEX="$3"
VCF_FILE="$4"

# Configuration Utilities
export PROJECT_DIR="./async_repli_project"
export SNPSPLIT_DIR="$PROJECT_DIR/tools/SNPsplit_v0.3.2"
export BLACKLIST_BED="$PROJECT_DIR/data/mm10-blacklist.v2.bed"
export PICARD_JAR="/sci/labs/itamarsi/britnyb/lab_files/bioinformatics/picard-tools-2.1.0/picard.jar"

STAR_POSTFIX="trim.Aligned.sortedByCoord.out"

echo "[RUNNING] Processing data track file: $FIRST_MATE"
SAMPLE_NAME=$(basename "$FIRST_MATE" _1.fastq.gz)

# Make sure directory outputs structure exists cleanly
mkdir -p "$CLONE_DIR/mapped"

# ==============================================================================
# STEP 1: MAPPING VIA STAR
# ==============================================================================
echo "[STAGE 1] Launching core alignment script via STAR..."

MAPPING_JOB_ID=$(sbatch --wrap="STAR --genomeDir $STAR_INDEX \
    --runThreadN 8 \
    --readFilesIn $FIRST_MATE \
    --outFileNamePrefix $CLONE_DIR/mapped/$SAMPLE_NAME. \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped None \
    --alignEndsType EndToEnd \
    --readFilesCommand zcat \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMatchNmin 0 \
    --outFilterMismatchNmax 2 \
    --outSAMattributes NH HI NM MD nM AS" \
    --mem=40GB -c 8 -t 20:00:00 -N 1 -o logs/slurm.STAR.$SAMPLE_NAME.%j.out | cut -f 4 -d' ')

echo "[JOB ASSIGNED] STAR Mapping Job ID: $MAPPING_JOB_ID"

# ==============================================================================
# STEP 2: MAPQ FILTERING & BLACKLIST ELIMINATION
# ==============================================================================
echo "[STAGE 2] Queueing MAPQ20 alignment filtration and blacklist filtering..."

MPQ_JOB_ID=$(sbatch --dependency=afterok:"$MAPPING_JOB_ID" --wrap="samtools view -bq 20 $CLONE_DIR/mapped/$SAMPLE_NAME.$STAR_POSTFIX.bam | \
    samtools sort -@ 2 - |\
    bedtools intersect -a - -b $BLACKLIST_BED -v \
    > $CLONE_DIR/mapped/$SAMPLE_NAME.$STAR_POSTFIX.mpq20.bam" \
    --mem=10GB -c 2 -t 04:00:00 -o logs/slurm.mpq.$SAMPLE_NAME.%j.out | cut -f 4 -d' ')

echo "[JOB ASSIGNED] MAPQ Filter Job ID: $MPQ_JOB_ID"

# ==============================================================================
# STEP 3: DUPLICATE REMOVAL VIA PICARD
# ==============================================================================
echo "[STAGE 3] Running duplicate fragment elimination via Picard MarkDuplicates..."

DUPS_JOB_ID=$(sbatch --dependency=afterok:"$MPQ_JOB_ID" --wrap="java -Xmx4g -jar $PICARD_JAR MarkDuplicates \
    INPUT=$CLONE_DIR/mapped/$SAMPLE_NAME.$STAR_POSTFIX.mpq20.bam \
    OUTPUT=$CLONE_DIR/mapped/$SAMPLE_NAME.$STAR_POSTFIX.mpq20.nodups.bam \
    METRICS_FILE=$CLONE_DIR/mapped/$SAMPLE_NAME.$STAR_POSTFIX.metrics.txt \
    REMOVE_DUPLICATES=true" \
    --mem=20GB -c 1 -t 04:00:00 -o logs/slurm.rmDups.$SAMPLE_NAME.%j.out | cut -f 4 -d' ')

echo "[JOB ASSIGNED] MarkDuplicates Job ID: $DUPS_JOB_ID"

# ==============================================================================
# STEP 4: ALLELE-SPECIFIC SORTING VIA SNPSPLIT
# ==============================================================================
echo "[STAGE 4] Converting variant configurations and parsing SNPSplit arrays..."

# Format inputs to match SNPsplit naming standards
awk '{OFS="\t"}; {print NR, $1, $2, 1, $4"/"$5}' "$VCF_FILE" > "$VCF_FILE.phased.txt"

SNPSPLIT_JOB_ID=$(sbatch --dependency=afterok:"$DUPS_JOB_ID" --wrap="$SNPSPLIT_DIR/SNPsplit \
    --snp_file $VCF_FILE.phased.txt \
    $CLONE_DIR/mapped/$SAMPLE_NAME*mpq20.nodups.bam \
    --conflicting \
    --weird \
    --paired \
    --no_sort" \
    --mem=25GB -c 2 -t 04:00:00 -o logs/slurm.snpSplit.$SAMPLE_NAME.%j.out | cut -f 4 -d' ')

echo "[JOB ASSIGNED] SNPsplit Processing Job ID: $SNPSPLIT_JOB_ID"
echo "[SUCCESS] All processing jobs for sample $SAMPLE_NAME successfully submitted."