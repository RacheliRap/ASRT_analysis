#!/bin/bash

#SBATCH -n 8
#SBATCH --mem=50G
#SBATCH -N 8
#SBATCH -o logs/slurm.orchestrator.%N.%j.out
#SBATCH -t 0-25:00
#SBATCH --mail-type=FAIL,END

# ==============================================================================
# ENVIRONMENT ENVIRONMENT & PATH CONFIGURATION
# ==============================================================================
export PROJECT_DIR="./async_repli_project"
export DATA_DIR="$PROJECT_DIR/data/ATAC-Seq"
export SNPSPLIT_DIR="$PROJECT_DIR/tools/SNPsplit_v0.3.2"
export REF_GENOME="/sci/data/reference_data/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index"

# Direct paths matching your experimental parameters
export BRITNY_SNPS="/sci/labs/itamarsi/rachelirap/britnyb/icore-data/lab_files/tor_seq_data"
export CLONE="clone2_merged"

# Export shared libraries path
export LD_LIBRARY_PATH="/usr/local/bioinfo/gsl/:$LD_LIBRARY_PATH"

# Load Core Cluster Modules
module load STAR
module load bedtools2/2.31.0

# ==============================================================================
# PIPELINE FUNCTIONS
# ==============================================================================

function trim_fastq() {
    echo "[INFO] Status: Running read adapter trimming via fastx_trimmer..."
    rm -f $DATA_DIR/$CLONE/fastq/*trim*
    
    # Safely decompress active fastq paths if necessary
    if ls $DATA_DIR/$CLONE/fastq/*.gz &>/dev/null; then
        gunzip $DATA_DIR/$CLONE/fastq/*.gz
    fi

    for f in $DATA_DIR/$CLONE/fastq/*Control2*001.fastq; do
        if [ -f "$f" ]; then
            local name=$(basename -s .fastq "$f")
            fastx_trimmer -f 16 -l 103 -i "$f" -o "$DATA_DIR/$CLONE/fastq/$name.trim.fastq"
        fi
    done
}

function refined_SNPs() {
    echo "[INFO] Status: Compiling and refining high-confidence variant tables..."
    mkdir -p $DATA_DIR/$CLONE
    
    local snps_file="$BRITNY_SNPS/$CLONE/pileups/combined_pileup_snps.nopar.vcf"
    local refined_snps="$BRITNY_SNPS/$CLONE/pileups/combined_pileup_snps.nopar.bed.refined"
    local header_file="$DATA_DIR/snpSplit.vcf_heder.vcf"
    
    cp "$header_file" "$DATA_DIR/$CLONE/combined_pileup_snps.snpSplit_header.vcf"
    awk 'NR > 54' "$snps_file" >> "$DATA_DIR/$CLONE/combined_pileup_snps.snpSplit_header.vcf"
    
    awk '{OFS="\t"}; NR > 54 {print $1, $2-1, $2, 1, $4"/"$5, NR}' "$snps_file" > "$DATA_DIR/$CLONE/combined_pileup_snps.snpSplit_format.bed"
    bedtools intersect -a "$DATA_DIR/$CLONE/combined_pileup_snps.snpSplit_format.bed" -b "$refined_snps" -wa > "$DATA_DIR/$CLONE/combined_pileup_snps.snpSplit_format.refined.bed"
    awk '{OFS="\t"}; {print $6, $1, $3, 1, $5}' "$DATA_DIR/$CLONE/combined_pileup_snps.snpSplit_format.refined.bed" > "$DATA_DIR/$CLONE/combined_pileup_snps.snpSplit_format.refined.txt"
    
    echo "[INFO] Status: VCF modifications complete."
}

function genome_preparation() {
    echo "[INFO] Status: Assembling N-Masked Reference Framework..."
    rm -f $DATA_DIR/CAST_EiJ_N-masked/*
    cd $DATA_DIR
	
    $SNPSPLIT_DIR/SNPsplit_genome_preparation \
        --nmasking \
        --strain CAST_EiJ \
        --reference_genome $REF_GENOME \
        --vcf_file $DATA_DIR/$CLONE/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf

    cat $(ls $DATA_DIR/CAST_EiJ_N-masked/chr*.fa | sort -V) > $DATA_DIR/CAST_EiJ_N-masked/genome.fa
}

function genome_index() {
    echo "[INFO] Status: Generating STAR Reference Indices..."
    mkdir -p $DATA_DIR/$CLONE/STAR_index
    
    STAR --runThreadN 8 \
        --runMode genomeGenerate \
        --genomeDir $DATA_DIR/$CLONE/STAR_index \
        --genomeFastaFiles $DATA_DIR/CAST_EiJ_N-masked/genome.fa 
}

function fastqToSNPs() {
    echo "[INFO] Status: Re-compressing fastq targets and scheduling cluster jobs..."
    
    awk 'NR > 69 {print $0}' $DATA_DIR/$CLONE/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf > $DATA_DIR/$CLONE/CAST_EiJ.mgp.v5.snps.dbSNP142.noHeader.vcf
    
    local vcf_no_header="$DATA_DIR/$CLONE/CAST_EiJ.mgp.v5.snps.dbSNP142.noHeader.vcf"
    local star_index="$DATA_DIR/$CLONE/STAR_index"
    
    # Package active trim runs before starting sequence matching
    if ls $DATA_DIR/$CLONE/fastq/*trim.fastq &>/dev/null; then
        gzip $DATA_DIR/$CLONE/fastq/*trim.fastq
    fi

    mkdir -p logs/

    for file in $DATA_DIR/$CLONE/fastq/*R1*trim*gz; do
        if [ -f "$file" ]; then
            echo "[SUBMIT] Launching dynamic tracking nodes for sample: $file"
            sbatch scripts/2_map_and_demultiplex.sh "$DATA_DIR/$CLONE" "$file" "$star_index" "$vcf_no_header"
        fi
    done
}

# ==============================================================================
# MAIN EXECUTION ROUTINE
# ==============================================================================
# Uncomment parameters relative to initialization dependencies
# trim_fastq
# refined_SNPs
# genome_preparation
# genome_index
fastqToSNPs