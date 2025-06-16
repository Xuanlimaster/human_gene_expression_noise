suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(scDblFinder)
})

dir_liver <- "human_adult_liver"
dir_kidney <- "human_adult_kidney"
dir_lung <- "human_adult_lung"

sce_liver_filtered <- readRDS(file.path(dir_liver, "sce_liver_filtered.rds"))
sce_kidney_filtered <- readRDS(file.path(dir_kidney, "sce_kidney_filtered.rds"))
sce_lung_filtered <- readRDS(file.path(dir_lung, "sce_lung_filtered.rds"))

# Normalize and process filtered single-cell RNA-seq data:
#   - Perform log-normalization of counts for downstream analysis
#   - Model gene variance to decompose technical/biological variation
#   - Select highly variable genes (HVGs) for dimensionality reduction
#   - Perform PCA on HVGs using IRLBA algorithm for efficiency
#   - Detect and filter doublets using scDblFinder
#   - Calculate UMAP embedding for visualization with tuned parameters
sce_normalize <- function(sce_tissue_filtered) {
  # Data normalization using scater::logNormCounts
  #   - Transforms raw counts to log-normalized counts
  #   - Stores results in new 'logcounts' assay slot
  sce_tissue_filtered_normalized <- logNormCounts(sce_tissue_filtered)
  
  # Decompose technical and biological variation:
  #   - Uses scran::modelGeneVar to estimate per-gene variance
  #   - Separates technical noise from biological variation
  dec <- modelGeneVar(sce_tissue_filtered_normalized)
  
  # Select top 2000 highly variable genes (HVGs):
  #   - Identifies genes with highest biological variation
  #   - Focuses downstream analysis on biologically relevant signals
  hvg <- getTopHVGs(dec, n = 2000)
  
  # Dimensional reduction pipeline (PCA -> UMAP):
  #   - PCA: Reduces dimensions while preserving maximum variance
  #     - Uses IRLBA algorithm for efficient computation on large datasets
  #     - Only uses selected HVGs to reduce technical noise
  #     - Stores results in reducedDim(sce, "PCA")
  sce_tissue_filtered_normalized <- runPCA(
    sce_tissue_filtered_normalized,
    subset_row = hvg,
    exprs_values = "logcounts",
    BSPARAM = BiocSingular::IrlbaParam()
  )
  
  # Detect and filter doublets using scDblFinder:
  #   - Identifies potential doublets based on artificial nearest-neighbor profiles
  #   - Uses batch information to account for batch-specific effects
  #   - Adds doublet classification metadata to colData
  sce_tissue_filtered_normalized <- scDblFinder(
    sce_tissue_filtered_normalized,
    samples = "BatchInfo"
  )
  
  # Filter to retain only singlets:
  #   - Removes cells classified as doublets (scDblFinder.class != "singlet")
  #   - Reduces potential artifacts from multiple cells captured in one droplet
  singlet <- sce_tissue_filtered_normalized$scDblFinder.class == "singlet"
  sce_tissue_filtered_normalized <- sce_tissue_filtered_normalized[, singlet] 
  
  # Calculate UMAP embedding:
  #   - Non-linear dimensionality reduction for visualization
  #   - Uses PCA results as input (50 dimensions by default)
  #   - Parameters tuned for balance between local/global structure:
  #     - n_neighbors = 15: Balances local/global structure preservation
  #     - min_dist = 0.1: Controls cluster compactness
  #     - metric = "cosine": Suitable for scRNA-seq data
  reducedDim(sce_tissue_filtered_normalized, "UMAP") <- calculateUMAP(
    sce_tissue_filtered_normalized, 
    dimred = "PCA",
    n_neighbors = 15,
    min_dist = 0.1,
    metric = "cosine"
  )
  
  # Return fully processed SingleCellExperiment object
  sce_tissue_filtered_normalized
}

sce_liver_filtered_normalized <- sce_normalize(sce_liver_filtered)
sce_kidney_filtered_normalized <- sce_normalize(sce_kidney_filtered)
sce_lung_filtered_normalized <- sce_normalize(sce_lung_filtered)

# Save normalized filtered SCE objects as RDS files
saveRDS(
  sce_liver_filtered_normalized,
  file = file.path(dir_liver, "sce_liver_filtered_normalized.rds")
)
saveRDS(
  sce_kidney_filtered_normalized,
  file = file.path(dir_kidney, "sce_kidney_filtered_normalized.rds")
)
saveRDS(
  sce_lung_filtered_normalized,
  file = file.path(dir_lung, "sce_lung_filtered_normalized.rds")
)