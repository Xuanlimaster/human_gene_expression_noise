suppressPackageStartupMessages({
  library(tidyverse)
  library(scater)
  library(batchelor)
})

dir_liver <- "human_adult_liver"
dir_kidney <- "human_adult_kidney"
dir_lung <- "human_adult_lung"

sce_liver_filtered_normalized <- 
  readRDS(file.path(dir_liver, "sce_liver_filtered_normalized.rds"))
sce_kidney_filtered_normalized <- 
  readRDS(file.path(dir_kidney, "sce_kidney_filtered_normalized.rds"))
sce_lung_filtered_normalized <-
  readRDS(file.path(dir_lung, "sce_lung_filtered_normalized.rds"))

# Process and quantify batch effects in normalized single-cell data:
#   - Creates output directory for QC results if not exists
#   - Visualizes batch effects using UMAP clustering
#   - Quantifies batch effects through centroid analysis
#   - Saves visualization plots and quantitative metrics
process_batch_effects <- function(sce_normalized,dir) {
  # Create QC output directory if not exists
  #   - dir_tissue: Parent directory path
  #   - qc_path: Subdirectory for storing QC plots ("QC plots")
  qc_path <- file.path(dir, "QC plots")
  if (!dir.exists(qc_path)) {
    dir.create(qc_path, recursive = TRUE)
  }
  
  # Plot 4: Batch effect visualization - UMAP colored by batch
  #   - Generates UMAP plot colored by BatchInfo metadata
  #   - Uses reduced dimensions stored in "UMAP" slot
  #   - Applies classic theme with centered bold title
  #   - Saves high-resolution PNG (600dpi) with specified dimensions
  p4 <- plotReducedDim(
    sce_normalized,
    dimred = "UMAP",
    colour_by = "BatchInfo",
    point_size = 1.5
  ) +
    ggtitle("UMAP by Batch") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Save UMAP charts
  #   - Stores plot in QC subdirectory
  #   - Filename: "UMAP Batch Clustering.png"
  #   - Dimensions: 8x6 inches
  umap_path <- file.path(qc_path, "UMAP Batch Clustering.png")
  ggsave(umap_path, p4, width = 8, height = 6, dpi = 600)
  
  # Batch Effect Quantification
  #   - Extracts UMAP coordinates from reducedDim slot
  #   - Calculates per-batch metrics:
  #     - Centroid positions (mean UMAP1/UMAP2)
  #     - Average distance from centroid (dispersion measure)
  #   - Adds tissue identifier from directory name
  #   - Rescales batch effect score to 0-1 range for comparability
  batch_score <- reducedDim(sce_normalized, "UMAP") %>% 
    as.data.frame() %>%
    mutate(Batch = sce_normalized$BatchInfo) %>% 
    group_by(Batch) %>% 
    summarise(
      Centroid_X = mean(UMAP1),
      Centroid_Y = mean(UMAP2),
      Avg_Distance = mean(sqrt((UMAP1 - mean(UMAP1))^2 + (UMAP2 - mean(UMAP2))^2)),
      .groups = "drop"
    ) %>%
    mutate(
      Tissue = basename(dir),
      Batch_Effect_Score = scales::rescale(Avg_Distance, to = c(0, 1))
    )
  
  # Save batch scoring forms
  #   - Stores quantitative metrics in CSV format
  #   - Filename: "batch_scores.csv"
  #   - Contains columns: Batch, Centroid_X/Y, Avg_Distance, Tissue, Score
  score_path <- file.path(qc_path, "batch_scores.csv")
  write_csv(batch_score, score_path)
}

process_batch_effects(sce_liver_filtered_normalized, dir_liver)
process_batch_effects(sce_kidney_filtered_normalized, dir_kidney)
process_batch_effects(sce_lung_filtered_normalized, dir_lung)