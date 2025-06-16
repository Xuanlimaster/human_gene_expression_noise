suppressPackageStartupMessages({
  library(tidyverse)
  library(scater)
})

dir_liver <- "human_adult_liver"
dir_kidney <- "human_adult_kidney"
dir_lung <- "human_adult_lung"

sce_liver <- readRDS(file.path(dir_liver, "sce_liver.rds"))
sce_kidney <- readRDS(file.path(dir_kidney, "sce_kidney.rds"))
sce_lung <- readRDS(file.path(dir_lung, "sce_lung.rds"))

# Process tissue QC metrics:
#   - Creates output directory structure for QC results
#   - Calculates per-cell QC metrics including library size, gene counts, and mitochondrial percentages
#   - Automatically detects outlier cells using MAD-based thresholds
#   - Generates summary statistics of filtered cells
#   - Saves raw metrics and filtering thresholds to CSV files
process_tissue_qc <- function(sce, dir) {
  # Create QC output directory if not exists
  #   - dir_tissue: Parent directory path
  #   - qc_path: Subdirectory for storing QC plots ("QC plots")
  qc_path <- file.path(dir, "QC plots")
  if (!dir.exists(qc_path)) {
    dir.create(qc_path, recursive = TRUE)
  }
  
  # Detect outlier cells using Median Absolute Deviation (MAD):
  #   - lib: Filters cells with unusually low library size (total UMIs)
  #   - gene: Filters cells with too few detected genes
  #   - mito: Filters cells with excessive mitochondrial percentage
  qc_stats <- sce %>% perCellQCMetrics() %>% as.data.frame() %>%
    cbind(colData(sce)[, c("subsets_Mito_sum",
                           "subsets_Mito_detected",
                           "subsets_Mito_percent")])
  
  # Detect outlier cells using Median Absolute Deviation (MAD):
  #   - lib: Filters cells with unusually low library size (total UMIs)
  #   - gene: Filters cells with too few detected genes
  #   - mito: Filters cells with excessive mitochondrial percentage
  thresholds <- list(
    lib = isOutlier(qc_stats$sum, log = TRUE, type = "lower"),
    gene = isOutlier(qc_stats$detected, log = TRUE, type = "lower"),
    mito = isOutlier(qc_stats$subsets_Mito_percent, type = "higher")
  )
  
  # Create aggregated summary of filtering results:
  #   - Metric: Type of QC metric being evaluated
  #   - Threshold: Calculated cutoff value for each metric
  #   - Cells_Filtered: Count of cells removed by each filter
  summary_df <- data.frame(
    Metric = c("Total UMIs", "Detected Genes", "Mitochondrial Percentage"),
    Threshold = c(
      attr(thresholds$lib, "thresholds")["lower"],
      attr(thresholds$gene, "thresholds")["lower"],
      attr(thresholds$mito, "thresholds")["higher"]
    ),
    Cells_Filtered = c(
      sum(thresholds$lib),
      sum(thresholds$gene),
      sum(thresholds$mito)
    )
  )
  
  # Save results to CSV files
  write.csv(qc_stats, file.path(qc_path, "qc_metrics.csv"), row.names = TRUE)
  write.csv(summary_df, file.path(qc_path, "qc_thresholds.csv"), row.names = FALSE)
}

process_tissue_qc(sce_liver, dir_liver)
process_tissue_qc(sce_kidney, dir_kidney)
process_tissue_qc(sce_lung, dir_lung)

# Quality Control Filtering (Based on plots and calculated thresholds)
qc.lib.liver <- sce_liver$sum > 2000          # Keep cells with >2000 total UMIs
qc.nexprs.liver <- sce_liver$detected > 500   # Keep cells detecting >500 genes
qc.mito.liver <- sce_liver$subsets_Mito_percent < 15    # MT% < 15%
keep_genes_liver <- rowSums(counts(sce_liver) > 0) > 20 # At least 20 cells

qc.lib.kidney <- sce_kidney$sum > 1000        # Keep cells with >2000 total UMIs
qc.nexprs.kidney <- sce_kidney$detected > 500 # Keep cells detecting >500 genes
qc.mito.kidney <- sce_kidney$subsets_Mito_percent < 10    # MT% < 10%
keep_genes_kidney <- rowSums(counts(sce_kidney) > 0) > 20 # At least 20 cells

qc.lib.lung <- sce_lung$sum > 1000          # Keep cells with >2000 total UMIs
qc.nexprs.lung <- sce_lung$detected > 500   # Keep cells detecting >500 genes
qc.mito.lung <- sce_lung$subsets_Mito_percent < 10    # MT% < 10%
keep_genes_lung <- rowSums(counts(sce_lung) > 0) > 20 # At least 20 cells

# Filter SingleCellExperiment object
sce_liver_filtered <- sce_liver[
  keep_genes_liver,
  qc.lib.liver & qc.nexprs.liver & qc.mito.liver
]

sce_kidney_filtered <- sce_kidney[
  keep_genes_kidney,
  qc.lib.kidney & qc.nexprs.kidney & qc.mito.kidney
]

sce_lung_filtered <- sce_lung[
  keep_genes_lung,
  qc.lib.lung & qc.nexprs.lung & qc.mito.lung
]

# Save filtered SCE objects as RDS files to their respective directories
saveRDS(
  sce_liver_filtered,
  file = file.path(dir_liver, "sce_liver_filtered.rds")
)
saveRDS(
  sce_kidney_filtered,
  file = file.path(dir_kidney, "sce_kidney_filtered.rds")
)
saveRDS(
  sce_lung_filtered,
  file = file.path(dir_lung, "sce_lung_filtered.rds")
)