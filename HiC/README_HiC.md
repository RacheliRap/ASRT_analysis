# Hi-C Contact Enrichment Analysis

Permutation-based statistical testing of chromatin contact enrichment between genomic groups.

## Overview

This notebook tests whether genomic regions from different groups show preferential same-group contacts using Hi-C data. The analysis converts Hi-C matrices to standard formats, calculates contact frequencies, performs permutation testing, and visualizes results as volcano plots.

## Requirements

```bash
# Install via pip (in notebook or terminal)
pip install numpy scikit-learn cooler cooltools hicmatrix hicexplorer hictk
```

**Software needed:**
- Python 3.7+
- hictk (file conversion)
- cooler (matrix manipulation)
- HiCExplorer (Hi-C operations)
- Standard Python: numpy, scipy, pandas, matplotlib, seaborn, joblib

**System requirements:**
- 12-16 GB RAM
- Multi-core CPU recommended

## Input Files

**1. Hi-C matrix:** `.hic` format (from 4DN, ENCODE, or your own data)

**2. Genomic intervals:** BED-like format with group labels
```
chr1    12459919    12859919    CAST
chr1    24143161    24443161    CAST
chr1    32693155    33293155    C57BL
```

## Usage

**1. Configure parameters** in the CONFIGURATION cell:
```python
DATA_DIR = "/path/to/data"
HIC_FILES = ["your_file.hic"]
RESOLUTION = 50000  # 50kb bins
N_PERMUTATIONS = 1000
```

**2. Run all cells** - The notebook will:
- Convert .hic → .cool → .h5 format
- Calculate within-group and between-group contacts
- Perform permutation testing (1000 iterations)
- Generate volcano plots

## Output

- **Converted matrices:** `*.cool` and `*.h5` files
- **Results CSV:** Enrichment scores, p-values, confidence intervals
- **Volcano plot:** PNG showing significant regions

## Key Parameters

- `RESOLUTION`: Bin size (default 50000 = 50kb)
- `N_PERMUTATIONS`: Number of random permutations (default 1000)
- `N_JOBS`: CPU cores for parallel processing (default 2)

## Typical Runtime

~30 minutes for one Hi-C file at 50kb resolution with 1000 permutations.

## Troubleshooting

**Memory error?** Reduce `N_PERMUTATIONS` or increase `RESOLUTION`

**File not found?** Check `DATA_DIR` path is correct

**Slow performance?** Increase `N_JOBS` to use more CPU cores

## Citation

Please cite the tools used:
- **hictk**: Marçais & Kingsford (2011) Bioinformatics
- **cooler**: Abdennur & Mirny (2020) Bioinformatics  
- **HiCExplorer**: Wolff et al. (2018) Nucleic Acids Research

## Contact

Rachel Rapoport | rachel.rapoport@mail.huji.ac.il
