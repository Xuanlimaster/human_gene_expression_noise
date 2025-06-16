suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(biomaRt)
  library(scater)
  library(SingleCellExperiment)
  })

dir_liver <- "human_adult_liver"
dir_kidney <- "human_adult_kidney"
dir_lung <- "human_adult_lung"

# Load human scRNA-seq data from MatrixMarket format (.mtx)
#   - 'cells' file contains cell barcodes (colnames)
#   - 'features' file contains gene identifiers (rownames)
counts_liver <- ReadMtx(
  mtx = file.path(dir_liver, "E-MTAB-10553.aggregated_counts.mtx.gz"),
  cells = file.path(dir_liver, "E-MTAB-10553.aggregated_counts.mtx_cols.gz"),
  features = file.path(dir_liver, "E-MTAB-10553.aggregated_counts.mtx_rows.gz")
  )

counts_kidney <- ReadMtx(
  mtx = file.path(dir_kidney, "E-CURD-119.aggregated_counts.mtx.gz"),
  cells = file.path(dir_kidney, "E-CURD-119.aggregated_counts.mtx_cols.gz"),
  features = file.path(dir_kidney, "E-CURD-119.aggregated_counts.mtx_rows.gz")
  )

counts_lung <- ReadMtx(
  mtx = file.path(dir_lung, "E-CURD-126.aggregated_counts.mtx.gz"),
  cells = file.path(dir_lung, "E-CURD-126.aggregated_counts.mtx_cols.gz"),
  features = file.path(dir_lung, "E-CURD-126.aggregated_counts.mtx_rows.gz")
  )

meta_liver <- fread(file.path(dir_liver, "E-MTAB-10553.cell_metadata.tsv")) %>%
  .[, 1:9] %>%
  mutate(BioSD_SAMPLE = str_extract(id, "^[^-]+")) %>%
  left_join(
    read.delim(file.path(dir_liver, "E-MTAB-10553.sdrf.txt")) %>%
      dplyr::select(
        BioSD_SAMPLE = "Comment.BioSD_SAMPLE.",
        BatchInfo = "Comment.ENA_SAMPLE."  # Comment.ENA_SAMPLE. as BatchInfo
        ) %>%
      distinct(BioSD_SAMPLE, BatchInfo, .keep_all = FALSE), 
    by = "BioSD_SAMPLE"
    ) %>%
  column_to_rownames("id")                 # Set cell IDs as rownames

meta_kidney <- fread(file.path(dir_kidney, "E-CURD-119.cell_metadata.tsv")) %>%
  .[, 1:12] %>%
  mutate(BioSD_SAMPLE = str_extract(id, "^[^-]+")) %>%
  left_join(
    read.delim(file.path(dir_kidney, "E-CURD-119.sdrf.txt")) %>%
      dplyr::select(
        BioSD_SAMPLE = "Comment.BioSD_SAMPLE.",
        BatchInfo = "Comment.ENA_SAMPLE."  # Comment.ENA_SAMPLE. as BatchInfo
        ) %>%
      distinct(BioSD_SAMPLE, BatchInfo, .keep_all = FALSE), 
    by = "BioSD_SAMPLE"
    ) %>%
  column_to_rownames("id")                 # Set cell IDs as rownames

meta_lung <- fread(file.path(dir_lung, "E-CURD-126.cell_metadata.tsv")) %>% 
  .[, 1:7] %>%
  mutate(BioSD_SAMPLE = str_extract(id, "^[^-]+")) %>%
  left_join(
    read.delim(file.path(dir_lung, "E-CURD-126.sdrf.txt")) %>%
      dplyr::select(
        BioSD_SAMPLE = "Comment.BioSD_SAMPLE.",
        BatchInfo = "Comment..ENA_SAMPLE." # Comment..ENA_SAMPLE. as BatchInfo
        ) %>%
      distinct(BioSD_SAMPLE, BatchInfo, .keep_all = FALSE), 
    by = "BioSD_SAMPLE"
    ) %>%
  column_to_rownames("id")                 # Set cell IDs as rownames

# Create SingleCellExperiment object:
#   - Ensure counts and metadata are matched by cell IDs
#   - Store raw counts in 'assays' slot
#   - Store cell metadata in 'colData' slot
#   - Retain only normal individuals
#   - Calculate and add mitpchondrial gene proportion
create_sce_obj <- function(counts, meta){
  counts <- counts[, match(rownames(meta), colnames(counts))]
  sce <- SingleCellExperiment(assays = list(counts = counts),
                              colData = meta)
  sce <- sce[, sce$disease == "normal"]
  # Connect to Ensembl and fetch mitochondrial genes
  ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  mt_genes <- getBM(
    attributes = c("ensembl_gene_id", "chromosome_name"),
    filters = "chromosome_name",
    values = "MT",
    mart = ensembl
  )
  # Validate gene expression across tissues
  valid_mt_genes <- intersect(mt_genes$ensembl_gene_id, rownames(sce))
  sce <- addPerCellQC(sce, subsets = list(Mito = valid_mt_genes))
  sce
}

sce_liver <- create_sce_obj(counts_liver, meta_liver)
sce_kidney <- create_sce_obj(counts_kidney, meta_kidney)
sce_lung <- create_sce_obj(counts_lung, meta_lung)

# Save SCE objects as RDS files to their respective directories
saveRDS(sce_liver, file = file.path(dir_liver, "sce_liver.rds"))
saveRDS(sce_kidney, file = file.path(dir_kidney, "sce_kidney.rds"))
saveRDS(sce_lung, file = file.path(dir_lung, "sce_lung.rds"))