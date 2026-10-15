# Allele-specific ATAC-seq signal in asynchronously replicating regions

Code to reproduce the allelic ATAC-seq density figure for clone 3 (cl3), a
CAST/EiJ × BALB/cJ hybrid mouse cell line: the distribution of
log2(BALB/c / CAST) ATAC signal in 100 kb bins, split by which allele replicates
early.

**Allele convention:** SNPsplit genome1 = BALB/cJ, genome2 = CAST/EiJ.
In `async_master_list.GRCm38.bed` and the bin annotation derived from it, the
BALB/c allele is labelled `C57BL`, because BALB/c carries the C57BL/6 reference
base at the SNPs used for allele assignment. Figure and statistics output use BALB/c.

Generated files (not in the repository): `data/coverage_genomes.Britnys_fig.mm10.ASRT_annotation.bed`
(script 3, step 1), `data/coverage_genomes.cl3_atac.bed` (script 3) and `results/` (script 4).

## Repository layout

```
scripts/
  1_prepare_and_submit.sh     SNP files, N-masked genome, STAR index, trimming; submits script 2
  2_map_and_demultiplex.sh    per-sample mapping, filtering, deduplication, allele assignment
  3_allelic_bin_counts.sh     bin annotation and allelic read counts per 100 kb bin
  4_plot_allelic_ratio.R      normalisation, filtering, figure and statistics
data/                         input files (see below)
results/                      created by script 4
```

## How to run

```bash
# edit the paths at the top of scripts 1 and 2 first
sbatch scripts/1_prepare_and_submit.sh                       # uncomment steps as needed
sbatch scripts/3_allelic_bin_counts.sh <FASTQ_DIR> <cl3_rep1> <cl3_rep2>
Rscript scripts/4_plot_allelic_ratio.R                       # from the repository root
```

All scripts are run from the repository root.

## Scripts

### 1. `1_prepare_and_submit.sh` — references and sample submission

Run once per data set; each step is a function, enabled at the bottom of the script.

1. **SNP files** (`create_snps`), from Mouse Genomes Project v5 SNPs (GRCm38, dbSNP142):
   - *Masking file*: all CAST/EiJ SNPs plus BALB/cJ-specific SNPs. Every position where
     either allele differs from the reference is masked, so that neither allele is
     favoured during mapping.
   - *Allele-assignment file*: CAST/EiJ SNPs at positions where BALB/cJ carries the
     reference base. At these positions reference = BALB (genome1) and alternative =
     CAST (genome2). Converted to SNPsplit format (`ID chr pos strand ref/alt`).
2. **N-masked genome** (`genome_preparation`): `SNPsplit_genome_preparation --nmasking`
   on GRCm38 with the masking file.
3. **STAR index** (`genome_index`) of the N-masked genome. The annotation GTF is
   included because the same index is used for the RNA-seq data.
4. **Trimming** (`trim_fastq`): `fastx_trimmer -f 16 -l 103` on both mates.
5. **Submission** (`submit_samples`): one `2_map_and_demultiplex.sh` job per read pair.

### 2. `2_map_and_demultiplex.sh` — per-sample processing (paired-end)

Each step is a dependent SLURM job:

| Step | Tool | Settings |
|---|---|---|
| Mapping | STAR | R1 + R2, end-to-end alignment (no soft-clipping), ≤ 2 mismatches |
| Filtering | samtools, bedtools | MAPQ ≥ 20 (unique), properly paired, ENCODE mm10 blacklist removed |
| Deduplication | Picard MarkDuplicates | `REMOVE_DUPLICATES=true` |
| Allele assignment | SNPsplit 0.3.2 | `--paired --conflicting --weird`, allele-assignment SNP file from step 1 |

Output: `<sample>.Aligned.sortedByCoord.out.mpq20.nodups.genome1.bam` (BALB/c) and
`.genome2.bam` (CAST).

### 3. `3_allelic_bin_counts.sh` — counts per bin

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

Output: `data/coverage_genomes.cl3_atac.bed` — columns `chr start end annotation
rep1_BALB rep1_CAST rep2_BALB rep2_CAST`, in the bin order of the annotation file.

### 4. `4_plot_allelic_ratio.R` — figure and statistics

- Each allele/replicate track is normalised to its library size; normalised signal is
  summed over replicates per allele, and log2(BALB/c / CAST) computed per bin.
- Bins are removed if total coverage is > 3 SD from the mean or ≤ 10 reads, or if they
  overlap an olfactory receptor gene cluster, an imprinted gene, or lie on chrX.
  Imprinted-gene coordinates are taken from Ensembl release 102 (the last GRCm38
  release) so that they match the mm10 bins.
- No pseudocount is added; bins with zero reads in one allele give ±Inf ratios (their
  number is reported in the statistics file).
- BALB/c-early and CAST-early bins are compared by one-sided Kolmogorov–Smirnov and
  Welch t-tests (H1: higher BALB/c / CAST ratio in BALB/c-early bins).

Output: `results/cl3_ATAC_allelic_ratio_density.pdf`,
`results/cl3_ATAC_allelic_ratio_stats.txt` (group sizes, means, tests, `sessionInfo()`).

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
