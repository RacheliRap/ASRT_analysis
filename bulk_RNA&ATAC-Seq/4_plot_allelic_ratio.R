# ==============================================================================
# 4_plot_allelic_ratio.R
#
# Allelic ATAC-seq or RNA-seq signal in clone 3 (CAST/EiJ x BALB/cJ hybrid) in 100 kb bins,
# split by asynchronous replication status: bins in regions where the BALB/c
# allele replicates early, where the CAST allele replicates early, and all other
# bins.
#
# Allele naming: SNPsplit genome1 = BALB/cJ, genome2 = CAST/EiJ.
# In the ASRT annotation file the BALB/c allele is labelled "C57BL", since BALB/c
# carries the C57BL/6 reference base at the SNPs used for allele assignment.
# The figure uses BALB/c.
#
# Input  (data/):
#   coverage_genomes.cl3_<atac|rna>.bed                    - from script 3
#   coverage_genomes.Britnys_fig.mm10.ASRT_annotation.bed  - annotated 100 kb bins (script 3)
#   olf_clusters.500000.mm10.bed                           - olfactory receptor clusters (excluded)
#   imprinted_genes.csv                                    - imprinted genes (excluded),
#                                                            column `converted_alias` = Ensembl gene ID
# Output (results/):
#   cl3_<ATAC|RNA>_allelic_ratio_density.pdf
#   cl3_<ATAC|RNA>_allelic_ratio_stats.txt
#
# All coordinates are GRCm38/mm10 with Ensembl chromosome names (1, 2, ..., X).
# Run from the repository root:
#   Rscript scripts/4_plot_allelic_ratio.R atac
#   Rscript scripts/4_plot_allelic_ratio.R rna
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(GenomicRanges)
  library(rtracklayer)
  library(biomaRt)
})

assay <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(assay) || !assay %in% c("atac", "rna"))
  stop("Usage: Rscript scripts/4_plot_allelic_ratio.R <atac|rna>")
assay_label <- c(atac = "ATAC", rna = "RNA")[[assay]]
prefix <- paste0("cl3_", assay_label, "_allelic_ratio")

data_dir <- "data"
out_dir  <- "results"
dir.create(out_dir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Regions excluded from the analysis
# ------------------------------------------------------------------------------
# Olfactory receptor gene clusters (monoallelic expression)
olf <- import.bed(file.path(data_dir, "olf_clusters.500000.mm10.bed"))

# Imprinted genes (parent-of-origin allelic bias).
# Coordinates are taken from Ensembl release 102, the last release on GRCm38,
# so that they match the mm10 bins.
imprinted_genes <- read.csv(file.path(data_dir, "imprinted_genes.csv"))
ensembl <- useEnsembl(biomart = "genes",
                      dataset = "mmusculus_gene_ensembl",
                      version = 102)
imp_list <- getBM(attributes = c("chromosome_name", "start_position",
                                 "end_position", "ensembl_gene_id"),
                  filters    = "ensembl_gene_id",
                  values     = unique(as.character(imprinted_genes$converted_alias)),
                  mart       = ensembl)
imp_list_gr <- makeGRangesFromDataFrame(imp_list,
                                        seqnames.field = "chromosome_name",
                                        start.field    = "start_position",
                                        end.field      = "end_position")

# ------------------------------------------------------------------------------
# 2. Allelic read counts per bin
#    genome1 = BALB/cJ, genome2 = CAST/EiJ (SNPsplit)
# ------------------------------------------------------------------------------
counts <- read.delim(file.path(data_dir, paste0("coverage_genomes.cl3_", assay, ".bed")), header = FALSE)
coords <- read.delim(file.path(data_dir, "coverage_genomes.Britnys_fig.mm10.ASRT_annotation.bed"),
                     header = FALSE, col.names = c("chr", "start", "end", "annt"))

# counts and annotation must describe the same bins in the same order
stopifnot(nrow(counts) == nrow(coords),
          all(counts$V1 == coords$chr),
          all(counts$V2 == coords$start),
          all(counts$V3 == coords$end))

data <- counts[, 5:8]
names(data) <- c("rep1_g1", "rep1_g2", "rep2_g1", "rep2_g2")

# Library-size normalisation of each allele/replicate track.
# Note: no pseudocount is added; bins with zero reads in one allele give
# +/-Inf log2 ratios (these fall outside the plotted range).
coverage <- data.frame(total_cov = rowSums(data))
for (col in names(data)) {
  coverage[col] <- data[col]
  coverage[paste(col, "norm", sep = "_")] <- coverage[col] / sum(coverage[col])
}
coverage <- cbind.data.frame(coverage, coords)

# ------------------------------------------------------------------------------
# 3. Bin filtering
# ------------------------------------------------------------------------------
# remove coverage outliers (> 3 SD from the mean) and low-coverage bins
coverage <- coverage[abs(coverage$total_cov - mean(coverage$total_cov)) <
                       3 * sd(coverage$total_cov), ]
coverage <- coverage[coverage$total_cov > 10, ]

# summed normalised signal per allele across replicates, and allelic ratio
coverage$genome1_norm <- rowSums(coverage[grep("g1_norm", names(coverage))])
coverage$genome2_norm <- rowSums(coverage[grep("g2_norm", names(coverage))])
coverage$log2_ratio   <- log2(coverage$genome1_norm / coverage$genome2_norm)

# exclude olfactory clusters, imprinted genes and chromosome X
coords.gr <- makeGRangesFromDataFrame(coverage[, c("chr", "start", "end")])
keep.coords <- !(coords.gr %over% olf) &
               !(coords.gr %over% imp_list_gr) &
               as.character(seqnames(coords.gr)) != "X"
coverage <- coverage[keep.coords, ]

# ------------------------------------------------------------------------------
# 4. Figure
# ------------------------------------------------------------------------------
coverage$annt <- factor(coverage$annt, levels = c("All", "C57BL", "CAST"))

group_labels <- c(All = "Other bins", C57BL = "BALB/c early", CAST = "CAST early")

p <- ggplot(coverage, aes(x = log2_ratio, fill = annt)) +
  geom_density(alpha = 0.65) +
  xlim(c(-2, 2)) +
  xlab(bquote(log[2] ~ "(BALB/c / CAST)" ~ .(assay_label) ~ "signal")) +
  ylab("Density") +
  scale_fill_manual(values = c(All = "dimgrey", C57BL = "red", CAST = "darkblue"),
                    breaks = names(group_labels), labels = unname(group_labels)) +
  theme_classic() +
  labs(fill = NULL) +
  ggtitle(paste("cl3", assay_label))

ggsave(file.path(out_dir, paste0(prefix, "_density.pdf")), p, width = 5, height = 4)

# ------------------------------------------------------------------------------
# 5. Statistics: BALB/c-early vs CAST-early bins
#    H1: BALB/c / CAST ratio is higher in BALB/c-early than in CAST-early bins.
#    (ks.test "less": the BALB/c-early CDF lies below the CAST-early CDF, i.e.
#     BALB/c-early values are stochastically larger.)
# ------------------------------------------------------------------------------
ratio.ASRT <- droplevels(coverage[coverage$annt != "All", ])
ratio.ASRT$annt <- factor(ratio.ASRT$annt, levels = c("C57BL", "CAST"))

sink(file.path(out_dir, paste0(prefix, "_stats.txt")))
cat("Assay:", assay_label, "\n")
cat("Groups: C57BL = BALB/c early, CAST = CAST early, All = other bins (not asynchronous)\n")
cat("Ratio: log2(BALB/c / CAST)\n\n")
cat("Bins retained after filtering:", nrow(coverage), "\n")
print(table(coverage$annt))
cat("\nNon-finite log2 ratios (zero counts in one allele) per group:\n")
print(tapply(!is.finite(coverage$log2_ratio), coverage$annt, sum))
cat("\nMean log2(BALB/c / CAST) per group:\n")
print(tapply(coverage$log2_ratio, coverage$annt, mean))
print(ks.test(log2_ratio ~ annt, data = ratio.ASRT, alternative = "less"))
print(t.test(log2_ratio ~ annt, data = ratio.ASRT, alternative = "greater"))
cat("\n")
print(sessionInfo())
sink()
