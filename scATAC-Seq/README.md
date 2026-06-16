# Allele-Specific Chromatin Accessibility Pipeline for Hybrid Mouse

This repository contains a complete, reproducible single-cell ATAC-seq pipeline designed to analyze parent-of-origin variation in hybrid mouse models ($CAST \times B6$). The workflow uses a multi-phase design, starting with data preprocessing in R (**Signac/Seurat**), moving to latent cell-state classification in Python via an **Expectation-Maximization (EM) Algorithm**, and returning to R for downstream clustering, statistical testing, and genomic visualization.

---

## Data Flow & Execution Steps

The pipeline operates as an integrated "round-trip" ecosystem across three scripts.

---

### Phase 0: Upstream Read Alignment & Sorting
**Script:** `scripts/0_upstream_alignment.sh`

Execute this script within an HPC SLURM queue environment. It resolves parental allocation mapping bias by creating N-masked synthetic genomes and isolating cell fragments by haplotypic origin.

- **Variant Processing:** Re-engineers reference configurations by resolving shared variations between $129S1$ and $CAST$ strains.
- **Genome Masking:** Calls `SNPsplit_genome_preparation` to flag known single-nucleotide variant zones as `N`-masked positions.
- **Haplotype Demultiplexing:** Evaluates Picard alignments, drops hard-clipped anomalies, and parses read IDs using variant anchors back into distinct maternal (`genome1`) and paternal (`genome2`) FASTQ arrays via `seqtk`.
- **Consensus Re-analysis:** Re-runs `cellranger-atac reanalyze` to evaluate bi-allelic peaks across both parental backgrounds.

---

### Phase 1: Preprocessing & Coordination
**Script:** `scripts/1_preprocessing.Rmd`

Execute this script first in RStudio. It ingests the raw 10x `.h5` matrices, single-cell metadata, and fragment tracks. It then:

- Discards background noise and low-quality cells using dynamic, metric-based distribution bounds.
- Intersects individual parental tracks against coordinates from the master `async_master_list` BED file.
- Combines and structures cross-allelic unique hits into consolidated count matrices.

**Outputs generated** — saves the integrated workspace state to `data/scATAC_preB_all_bi.RData` and exports two clean CSV matrices to feed the EM step:

| File | Description |
|------|-------------|
| `data/c_merged_bi_preB.csv` | Cell-by-region count matrix |
| `data/regions_numeric_merged_bi_preB.csv` | Numeric region identifiers |

---

### Phase 2: Latent EM Classification
**Script:** `scripts/2_em_classification.ipynb`

Open this notebook (`scATAC_seq_EM_final_publication.ipynb`) in Google Colab or your local Jupyter instance. This notebook:

- Loads the `c_merged_bi_preB.csv` and `regions_numeric_merged_bi_preB.csv` files generated in Phase 1.
- Models read counts using region- and state-specific Poisson distributions to account for specific allelic biases (e.g., CAST-early vs. B6-early regions).
- Runs an Expectation-Maximization (EM) algorithm to identify three cell states:

| State | Label | Description |
|-------|-------|-------------|
| State 1 | CAST > B6 | Cells with preferential chromatin accessibility on CAST alleles |
| State 2 | B6 > CAST | Cells with preferential chromatin accessibility on B6 alleles |
| State 3 | Noise | Low-quality cells or cells without a clear allelic preference |

**Outputs generated** — saves the final assignments matrix as `data/cell_states_bi_preB.csv`.

---

### Phase 3: Validation, Embeddings, and Browser Tracks
**Script:** `scripts/3_downstream_analysis.Rmd`

Return to R and execute the final script to map the EM classifications back to the single-cell objects for visualization and profiling:

- **Statistical Testing:** Performs independent Kolmogorov-Smirnov (KS) and Student's t-tests across allelic logs to track parent-of-origin skewing.
- **Genomic Browser Tracks:** Formats and exports spatial coverage data into a `data/coverage_keep_bi.bedgraph` file, ready to be visualized on the UCSC Genome Browser or IGV.
- **Semi-Supervised Manifold Alignment:** Runs custom LSI/SVD structural reductions on the filtered ratio matrix and passes the categorical EM classification vectors directly into `uwot::umap` to generate guided cluster charts.
