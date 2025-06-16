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

# Create and save QC plots for a given SCE object
generate_qc_plots <- function(sce, dir) {
  # Create QC output directory if not exists
  #   - dir_tissue: Parent directory path
  #   - qc_path: Subdirectory for storing QC plots ("QC plots")
  qc_path <- file.path(dir, "QC plots")
  if (!dir.exists(qc_path)) {
    dir.create(qc_path, recursive = TRUE)
  }
  
  # Plot 1: Library Size vs Detected Genes
  #   - x: Total UMI counts (log10-scaled)
  #   - y: Number of detected genes per cell
  #   - color: Batch information (technical replicates)
  #   - Visualizes potential low-quality cells (bottom-left quadrant)
  p1 <- plotColData(
    sce,
    x = "sum",                    # Total UMI counts per cell
    y = "detected",               # Number of detected genes
    colour_by = "BatchInfo"
  ) +
    labs(
      title = "Library Size vs Genes Detected",
      x = "Total UMI Counts (×10³)",
      y = "Number of Detected Genes",
      color = "ENA Sample"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(face = "bold")
    ) +
    scale_x_continuous(labels = function(x) x/1000) +
    scale_color_viridis_d() # Color-blind friendly palette
  # Save Plot 1
  p1_path <- file.path(qc_path, "BatchEffect_UMI_vs_Genes.png")
  ggsave(p1_path, p1, width = 10, height = 6, dpi = 600)
  
  # Plot 2: Mitochondrial Gene Percentage
  #   - x: Total UMI counts (log10-scaled)
  #   - y: Percentage of mitochondrial reads (MT%)
  #   - Red dashed line: Typical QC threshold (10%)
  #   - High MT% indicates cell stress/apoptosis
  p2 <- plotColData(
    sce, 
    x = "sum",                   # Total UMI counts
    y = "subsets_Mito_percent",  # MT gene percentage
    colour_by = "BatchInfo"
  ) +
    labs(
      title = "Mitochondrial Gene Percentage",
      x = "Total UMI Counts (×10³)",
      y = "MT Gene % of Total Counts",
      color = "ENA Sample"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(face = "bold")
    ) +
    scale_x_continuous(labels = function(x) x/1000) +
    scale_color_viridis_d() +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") # QC threshold
  # Save Plot 2
  p2_path <- file.path(qc_path, "BatchEffect_UMI_vs_MTpercent.png")
  ggsave(p2_path, p2, width = 10, height = 6, dpi = 600)
  
  # Plot 3: Combined QC Metrics
  #   - x: Total UMI counts (log10-scaled)
  #   - y: Number of detected genes per cell
  #   - color: Mitochondrial gene percentage (MT%)
  #     - Blue: Low MT% (<5%) - healthy cells
  #     - Grey: Intermediate MT% (5-10%) - borderline cells
  #     - Red: High MT% (>10%) - potentially stressed/dying cells
  #   - Visualizes the relationship between library size, gene detection, and mitochondrial content
  #   - Helps identify low-quality cells in bottom-left quadrant (low UMI + low genes + high MT%)
  p3 <- plotColData(
    sce,
    x = "sum",                         # Total UMI counts per cell
    y = "detected",                    # Number of detected genes
    colour_by = "subsets_Mito_percent" # Mitochondrial percentage
  ) + 
    labs(
      title = "Library Size vs Genes Detected",
      subtitle = "Colored by Mitochondrial Gene Percentage (MT%)",
      x = "Total UMI Counts (×10³)",
      y = "Number of Detected Genes",
      color = "MT Gene %"
    ) +
    scale_color_gradient2(
      low = "blue", 
      mid = "grey", 
      high = "red", 
      midpoint = 10
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(face = "bold")
    )
  # Save Plot 3
  p3_path <- file.path(qc_path, "QC_Filtering_Matrix.png")
  ggsave(p3_path, p3, width = 10, height = 6, dpi = 600)
}

# Generate plots for each tissue
generate_qc_plots(sce_liver, dir_liver)
generate_qc_plots(sce_kidney, dir_kidney)
generate_qc_plots(sce_lung, dir_lung)