suppressPackageStartupMessages({
  library(BASiCS)
})

dir_kidney <- "human_adult_kidney"

sce_kidney_filtered_normalized <- 
  readRDS(file.path(dir_kidney, "sce_kidney_filtered_normalized.rds"))

# Run BASiCS MCMC
chain <- BASiCS_MCMC(
  sce_kidney_filtered_normalized,
  N = 30000,                     # Total iterations
  Thin = 15,                     # Keep every 15th sample
  Burn = 15000,                  # Initial samples discarded
  Regression = TRUE,             # Account for technical noise via regression
  WithSpikes = FALSE,            # No spike-in controls available
  PriorParam = BASiCS_PriorParam(
    sce_kidney_filtered_normalized,
    PriorMu = "EmpiricalBayes"   # Empirical Bayes priors for mean expression
  ),
  RunName = "human_adult_kidney", # Unique identifier for output
  Threads = 22,                  # Parallel processing
  StoreChains = TRUE,            # Save complete MCMC chains
  StoreDir = dir_kidney
)