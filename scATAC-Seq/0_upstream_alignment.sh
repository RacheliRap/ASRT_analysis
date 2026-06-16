#!/bin/bash

#SBATCH --ntasks=1          # Single task execution
#SBATCH --cpus-per-task=16  # Symmetric multiprocessing cores
#SBATCH --mem=128G          # Memory allocation
#SBATCH -o logs/slurm.%N.%j.out
#SBATCH -t 50:00:00         # Walltime threshold (D-HH:MM:MM)
#SBATCH --mail-type=FAIL,END

# ==============================================================================
# ENVIRONMENT ENVIRONMENT & CONFIGURATION
# ==============================================================================
# Define project structures (adjust these paths to reflect your cluster mount points)
export PROJECT_DIR="./rt_splicing_project"
export DATA_DIR="$PROJECT_DIR/data/scATAC_ESC"
export SNPSPLIT_DIR="$PROJECT_DIR/tools/SNPsplit_v0.3.2"
export CELLRANGER_DIR="$PROJECT_DIR/tools/cellranger-atac-2.1.0"

# Files & References
export REF_GENOME="/sci/data/reference_data/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index"
export TARGET_GTF="/sci/data/reference_data/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/masked.gtf"
export VCF_FILE="$PROJECT_DIR/data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"
export PICARD_JAR="$PROJECT_DIR/tools/picard.jar"

export ID="scATAC_ESC_all"

# Load Required Dynamic Environment Modules
module load cutadapt/4.4
module load samtools/1.17
module load bedtools2/2.31.0
module load star/2.7.11b-x86_64-gcc-12.2.0-u5gm272

export PATH="$CELLRANGER_DIR:$PATH"

# ==============================================================================
# WORKFLOW FUNCTIONS
# ==============================================================================

function create_snps() {
    echo "[INFO] Running: Variant masking file creation..."
    awk 'NR < 70' $VCF_FILE > $DATA_DIR/CAST_EiJ.mgp.v5.snps.dbSNP142.header
    
    bedtools intersect -v -a $VCF_FILE -b $DATA_DIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf > $DATA_DIR/CAST_EiJ.no_129S1_SvImJ.vcf
    bedtools intersect -v -b $VCF_FILE -a $DATA_DIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf > $DATA_DIR/129S1_SvImJ.no_CAST.vcf
    bedtools intersect -b $VCF_FILE -a $DATA_DIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf -wo > $DATA_DIR/129S1_SvImJ.overlap_CAST.vcf

    cat $DATA_DIR/CAST_EiJ.no_129S1_SvImJ.vcf \
        $DATA_DIR/129S1_SvImJ.no_CAST.vcf \
        $DATA_DIR/129S1_SvImJ.overlap_CAST.vcf | bedtools sort | \
        cat $DATA_DIR/CAST_EiJ.mgp.v5.snps.dbSNP142.header - > $DATA_DIR/129S1_SvImJ.CAST_EiJ.SNPs_for_mask.vcf
}

function genome_preparation() {
    echo "[INFO] Running: N-Masked Reference Assembly Construction..."
    mkdir -p $DATA_DIR/CAST_EiJ_N-masked
    rm -f $DATA_DIR/CAST_EiJ_N-masked/*
    
    cd $DATA_DIR
    $SNPSPLIT_DIR/SNPsplit_genome_preparation \
        --nmasking \
        --strain CAST_EiJ \
        --reference_genome $REF_GENOME \
        --vcf_file $DATA_DIR/129S1_SvImJ.CAST_EiJ.SNPs_for_mask.vcf

    cat $(ls $DATA_DIR/CAST_EiJ_N-masked/chr*.fa | sort -V) > $DATA_DIR/CAST_EiJ_N-masked/genome.fa
}

function genome_index() {
    echo "[INFO] Running: Compiling CellRanger Target Indices..."
    local gtf_name=$(basename -s .gtf $TARGET_GTF)
    grep -vE '^J|^G' $TARGET_GTF > $DATA_DIR/$gtf_name.filter.gtf
    
    cellranger-atac mkref --config=config_cellRanger.ESC.txt
}

function cell_ranger() {
    echo "[INFO] Running: Initial Biallelic Alignment Count Pipeline..."
    cd $DATA_DIR
    cellranger-atac count --id=$ID \
        --sample=JB_1107,JB_1107_2 \
        --reference=$PROJECT_DIR/analysis/GRCm38_ESC \
        --fastqs=$DATA_DIR/fastqShallow,$DATA_DIR/fastqDeep \
        --localcores=16 \
        --localmem=128
}

function addTags() {
    echo "[INFO] Running: Picard MD/UQ/NM sequence tag refinement..."
    local masked_genome=$DATA_DIR/CAST_EiJ_N-masked/genome.fa
   
    for bam in $DATA_DIR/$ID/outs/*.bam; do
        local out_bam=$bam.setTags.bam

        java -Xmx128G -jar $PICARD_JAR SetNmMdAndUqTags \
            --REFERENCE_SEQUENCE $masked_genome \
            --INPUT $bam \
            --OUTPUT $out_bam \
            --USE_JDK_DEFLATER true \
            --USE_JDK_INFLATER true

        # Clean hard-clipped artifacts via dynamic CIGAR filtering
        samtools view -@ 16 -h $out_bam \
            | awk '{if($0 ~ /^@/ || $6 ~ /^[0-9]*[MIDN]+[0-9MIDN]*$/) {print $0}}' \
            | samtools view -Sb - > tmp_clean.bam && mv tmp_clean.bam $out_bam
    done
}

function snpSplit() {
    echo "[INFO] Running: SNPsplit Haplotypic Alignment Sorting..."
    grep -v -E "^#" $VCF_FILE | awk '{OFS="\t"}; {print NR, $1, $2, 1, $4"/"$5}' > $VCF_FILE.phased.txt
		
    $SNPSPLIT_DIR/SNPsplit --snp_file $VCF_FILE.phased.txt \
        $DATA_DIR/$ID/outs/*setTags.bam \
        --paired \
        --conflicting \
        --weird 
}

function alleleToFastq() {
    echo "[INFO] Running: De-multiplexing sorted BAM variants back into Fastq..."
    samtools view -h $DATA_DIR/$ID/outs/*genome1*bam | awk 'BEGIN{FS="\t";OFS="\t"} {if(substr($1,1,1)!="@"){print $1}}' | LC_ALL=C sort -u > $DATA_DIR/read_ids.genome1.txt
    samtools view -h $DATA_DIR/$ID/outs/*genome2*bam | awk 'BEGIN{FS="\t";OFS="\t"} {if(substr($1,1,1)!="@"){print $1}}' | LC_ALL=C sort -u > $DATA_DIR/read_ids.genome2.txt

    datasets=("fastqShallow" "fastqDeep")
    for dataset in "${datasets[@]}"; do
        local allelic_dir="${DATA_DIR}/${dataset}Allelic"
        mkdir -p "${allelic_dir}/genome1" "${allelic_dir}/genome2"
        rm -f "${allelic_dir}/genome1/*" "${allelic_dir}/genome2/*"

        for f in "${DATA_DIR}/${dataset}/JB_"*; do
            local name=$(basename -s .fastq.gz "$f")
            for g in "genome1" "genome2"; do
                zcat "$f" | seqtk subseq - "${DATA_DIR}/read_ids.$g.txt" | gzip > "${allelic_dir}/$g/$name.fastq.gz"
            done
        done
    done
}

function re_cellranger() {
    echo "[INFO] Running: Secondary Allele-Specific Quantitation Profiles..."
    cd $DATA_DIR
    cellranger-atac count --id="${ID}_genome1" \
        --description="Maternal_Genome1" \
        --sample=JB_1107,JB_1107_2 \
        --fastqs=$DATA_DIR/fastqShallowAllelic/genome1,$DATA_DIR/fastqDeepAllelic/genome1 \
        --reference=$PROJECT_DIR/analysis/GRCm38_preB \
        --localcores=16 \
        --localmem=128

    cellranger-atac count --id="${ID}_genome2" \
        --description="Paternal_Genome2" \
        --sample=JB_1107,JB_1107_2 \
        --fastqs=$DATA_DIR/fastqShallowAllelic/genome2,$DATA_DIR/fastqDeepAllelic/genome2 \
        --reference=$PROJECT_DIR/analysis/GRCm38_preB \
        --localcores=16 \
        --localmem=128
}

function reanalyze() {
    echo "[INFO] Running: Re-analysis mapping to a standardized consensus shared peak-set..."
    cd $DATA_DIR
    cellranger-atac reanalyze --id="${ID}_genome1_biPeaks" \
        --peaks=$DATA_DIR/$ID/outs/peaks.bed \
        --fragments=$DATA_DIR/"${ID}_genome1"/outs/fragments.tsv.gz \
        --reference=$PROJECT_DIR/analysis/GRCm38_preB \
        --localcores=16 \
        --localmem=128

    cellranger-atac reanalyze --id="${ID}_genome2_biPeaks" \
        --peaks=$DATA_DIR/$ID/outs/peaks.bed \
        --fragments=$DATA_DIR/"${ID}_genome2"/outs/fragments.tsv.gz \
        --reference=$PROJECT_DIR/analysis/GRCm38_preB \
        --localcores=16 \
        --localmem=128
}

# ==============================================================================
# MAIN EXECUTION ROUTINE
# ==============================================================================
# Uncomment parameters as necessary depending on your checkpoint restarts
create_snps
genome_preparation
genome_index
cell_ranger
addTags
snpSplit
alleleToFastq 
re_cellranger
reanalyze

echo "[SUCCESS] Upstream Pipeline Run Completed."
