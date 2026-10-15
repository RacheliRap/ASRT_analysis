# Allele-specific ATAC-seq and RNA-seq signal in asynchronously replicating regions

Code to reproduce the allelic ATAC-seq and RNA-seq density figures for clone 3
(cl3), a CAST/EiJ × BALB/cJ hybrid mouse cell line: the distribution of
log2(BALB/c / CAST) signal in 100 kb bins, split by which allele replicates early.

| | ATAC-seq | RNA-seq |
|---|---|---|
| Reads | paired-end | single-end |
| Trimming | `fastx_trimmer -f 16 -l 103` | none |
| Filtering | MAPQ ≥ 20, properly paired, blacklist | MAPQ ≥ 20, blacklist |
| Duplicate removal | yes | no |
| Replicates | 2 | 2 |

The assay is given explicitly (`atac` or `rna`) to scripts 2, 3 and 4.

**Allele convention:** SNPsplit genome1 = BALB/cJ, genome2 = CAST/EiJ.
In `async_master_list.GRCm38.bed` and the bin annotation derived from it, the
BALB/c allele is labelled `C57BL`, because BALB/c carries the C57BL/6 reference
base at the SNPs used for allele assignment. Figure and statistics output use BALB/c.

Generated files (not in the repository): `data/coverage_genomes.Britnys_fig.mm10.ASRT_annotation.bed`
(script 3, step 1), `data/coverage_genomes.cl3_atac.bed` and `data/coverage_genomes.cl3_rna.bed`
(script 3), and `results/` (script 4).

## Repository layout

```
scripts/
  1_prepare_and_submit.sh     SNP files, N-masked genome, STAR index, trimming; submits script 2
  2_map_and_demultiplex.sh    per-sample mapping, filtering, deduplication (ATAC), allele assignment
  3_allelic_bin_counts.sh     bin annotation and allelic read counts per 100 kb bin
  4_plot_allelic_ratio.R      normalisation, filtering, figure and statistics
data/                         input files (see below)
results/                      created by script 4
```

## How to run

```bash
# edit the paths at the top of scripts 1 and 2 first; run everything from the repository root
mkdir -p logs                                      # SLURM log folder; must exist before submitting
sbatch scripts/1_prepare_and_submit.sh            # uncomment steps as needed; submits script 2

sbatch scripts/3_allelic_bin_counts.sh atac <ATAC_FASTQ_DIR> <cl3_atac_rep1> <cl3_atac_rep2>
sbatch scripts/3_allelic_bin_counts.sh rna  <RNA_FASTQ_DIR>  <cl3_rna_rep1>  <cl3_rna_rep2>

Rscript scripts/4_plot_allelic_ratio.R atac
Rscript scripts/4_plot_allelic_ratio.R rna
```

All scripts are run from the repository root. Script 2 is submitted by script 1;
scripts 3 and 4 are run once all samples have finished. The two `sbatch` calls of script 3
(ATAC and RNA) can run at the same time.

## Scripts

### 1. `1_prepare_and_submit.sh` — references and sample submission

Run once per data set; each step is a function, enabled at the bottom of the script.

1. **SNP files** (`create_snps`), from Mouse Genomes Project v5 SNPs (GRCm38, dbSNP142):
   - *Masking file*: union of CAST-only SNPs, BALB-only SNPs and positions variable in
     both strains. Every position where either allele differs from the reference is
     masked, so that neither allele is favoured during mapping.
     (`SNPsplit_genome_preparation` masks the high-confidence, homozygous calls, FI = 1.)
   - *Allele-assignment file*: CAST/EiJ SNPs at positions where BALB/cJ carries the
     reference base, converted to SNPsplit format (`ID chr pos strand ref/alt`).
     SNPsplit scores a read only at masked (N) positions listed in this file: the
     reference base counts for genome1 (BALB/c), the CAST base for genome2. Masked
     BALB-only and shared positions are not listed and therefore do not assign reads.
   One set of SNP files, genome and index is used for both ATAC-seq and RNA-seq.
2. **N-masked genome** (`genome_preparation`): `SNPsplit_genome_preparation --nmasking`
   on GRCm38 with the masking file.
3. **STAR index** (`genome_index`) of the N-masked genome. The annotation GTF is
   included because the same index is used for the RNA-seq data.
4. **Trimming** (`trim_fastq`): `fastx_trimmer -f 16 -l 103` on both ATAC mates.
   RNA-seq reads are not trimmed.
5. **Submission** (`submit_atac`, `submit_rna`): one `2_map_and_demultiplex.sh` job
   per ATAC read pair (`atac`) and per RNA-seq read file (`rna`).

### 2. `2_map_and_demultiplex.sh` — per-sample processing

`2_map_and_demultiplex.sh <atac|rna> <FASTQ_DIR> <R1.fastq.gz> <STAR_INDEX> <SNP_TXT>`.
Each step is a dependent SLURM job:

| Step | Tool | ATAC-seq | RNA-seq |
|---|---|---|---|
| Mapping | STAR | R1 + R2 | R1 |
| | | end-to-end alignment (no soft-clipping), ≤ 2 mismatches | same |
| Filtering | samtools, bedtools | MAPQ ≥ 20 (unique), properly paired, ENCODE mm10 blacklist removed | MAPQ ≥ 20, blacklist removed |
| Deduplication | Picard MarkDuplicates | `REMOVE_DUPLICATES=true` | not done |
| Allele assignment | SNPsplit 0.3.2 | `--paired --conflicting --weird`; SNPsplit sorts the BAM by read name itself | `--conflicting --weird` |

Output (genome1 = BALB/c, genome2 = CAST):
- ATAC: `<sample>.Aligned.sortedByCoord.out.mpq20.nodups.genome{1,2}.bam`
- RNA: `<sample>.Aligned.sortedByCoord.out.mpq20.genome{1,2}.bam`

### 3. `3_allelic_bin_counts.sh` — counts per bin

`3_allelic_bin_counts.sh <atac|rna> <FASTQ_DIR> <SAMPLE_REP1> <SAMPLE_REP2>`; sample
names are those printed by script 2.

1. Annotates each 100 kb bin with the early-replicating allele from
   `async_master_list.GRCm38.bed` (`C57BL` = BALB/c, `CAST`), or `All` for bins
   outside asynchronously replicating regions:
   ```bash
   bedtools intersect -a coverage_genomes.Britnys_fig.mm10.bed -b async_master_list.GRCm38.bed -wao | \
     awk 'BEGIN{OFS="\t"} $5=="-1" {print $1,$2,$3,"All"} $5!="-1" {print $1,$2,$3,$7}' \
     > coverage_genomes.Britnys_fig.mm10.ASRT_annotation.bed
   ```
   The script stops if any bin overlaps more than one asynchronous region.
2. Sorts and indexes the genome1 / genome2 BAMs of the two cl3 replicates.
3. Counts reads per bin with `bedtools multicov`.

Output: `data/coverage_genomes.cl3_<atac|rna>.bed` — columns `chr start end annotation
rep1_BALB rep1_CAST rep2_BALB rep2_CAST`, in the bin order of the annotation file.

### 4. `4_plot_allelic_ratio.R` — figure and statistics

`Rscript scripts/4_plot_allelic_ratio.R <atac|rna>`; the same analysis is applied to both assays.

- Each allele/replicate track is normalised to its library size; normalised signal is
  summed over replicates per allele, and log2(BALB/c / CAST) computed per bin.
- Bins are removed if total coverage is > 3 SD from the mean or ≤ 10 reads, or if they
  overlap an olfactory receptor gene cluster, an imprinted gene, or lie on chrX.
  Imprinted-gene coordinates are taken from Ensembl release 102 (the last GRCm38
  release) so that they match the mm10 bins.
- Imprinted-gene coordinates are downloaded with biomaRt, so script 4 needs internet access.
- No pseudocount is added; bins with zero reads in one allele give ±Inf ratios (their
  number is reported in the statistics file).
- BALB/c-early and CAST-early bins are compared by one-sided Kolmogorov–Smirnov and
  Welch t-tests (H1: higher BALB/c / CAST ratio in BALB/c-early bins).

Output: `results/cl3_<ATAC|RNA>_allelic_ratio_density.pdf` and
`results/cl3_<ATAC|RNA>_allelic_ratio_stats.txt` (group sizes, means, tests, `sessionInfo()`).

## Input files (`data/`)

| File | Description | Source |
|---|---|---|
| `CAST_EiJ.mgp.v5.snps.dbSNP142.vcf`, `BALB_cJ.mgp.v5.snps.dbSNP142.vcf` | Per-strain SNPs vs GRCm38 (script 1) | [Mouse Genomes Project](https://www.mousegenomes.org/) v5 |
| `mm10-blacklist.v2.bed` | ENCODE blacklist (script 2) | [Boyle lab](https://github.com/Boyle-Lab/Blacklist) |
| `coverage_genomes.Britnys_fig.mm10.bed` | 100 kb genomic bins (mm10) | <!-- TODO: describe origin --> |
| `async_master_list.GRCm38.bed` | Asynchronously replicating regions; column 4 = early-replicating allele (`CAST` / `C57BL` = BALB/c) | This study |
| `olf_clusters.500000.mm10.bed` | Olfactory receptor gene clusters, excluded | <!-- TODO: describe origin --> |
| `imprinted_genes.csv` | Mouse imprinted genes; `converted_alias` = Ensembl gene ID | [Geneimprint](https://www.geneimprint.com/site/genes-by-species) |

Raw sequencing data: <!-- TODO: GEO / SRA accession -->.

## Software

STAR, samtools, bedtools 2.31.0, Picard 2.1.0, SNPsplit 0.3.2, FASTX-Toolkit;
R with ggplot2, GenomicRanges, rtracklayer and biomaRt (R package versions are
written to the statistics file).
