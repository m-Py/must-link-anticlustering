
# Author: Martin Papenberg
# Year: 2025

## This file implements the simulation reported in the manuscript
# "Anticlustering for Sample Allocation To Minimize Batch Effects"
# that compares anticlustering, OSAT and PSBA (propensity score batch assignment)

# It uses the script "functions.R", which implements some helper methods (e.g. wrappers for OSAT and PSBA methods)

library(anticlust) # I used version 0.8.9-1
library(OSAT) # I got version OSAT_1.52.0 from Bioconductor
library(here) # for managing relative paths in the project
source(here("Simulation", "functions.R")) # functions used for the simulation 

sessionInfo()

set.seed(100) # for reproducibility; "uncomment" for running different data sets
# Use 10000 data sets for final simulation. If done in one setting, should run in about 
# 3 days on a personal computer (mostly because PSBA is really slow).
nsim <- 10000

# Unfortunately, the PS implementation needs a seed. It also sets seeds 
# repeatedly, which can easily mess with the logic of the remaining 
# simulation. For this reason I generate initial seeds 
# at the beginning for all simulation iterations (see below).
seeds_ps <- sample(2147483647, nsim)

start_all <- Sys.time()
for (i in 1:nsim) {
  cat("Working iteration", i, "\n")
  cat("Time passed", as.character(round(as.numeric(difftime(Sys.time(), start_all, units = "s")), 2)), "s\n")
  K <- sample(c(2, 4, 10), size = 1)
  sample_sizes <- 50:500 
  N <- sample(sample_sizes[sample_sizes %% K == 0], 1) # sample size
  M  <- sample(2:5, size = 1) # number of variables
  P <- sample(2:5, size = 1) # number of classes for each category
  
  data <- generate_categorical_data(N, M, P)
  
  start <- Sys.time()
  OSAT <- my_osat(
    data,
    K = K, 
    nSim = 5000
  )
  OSAT_t <- as.numeric(difftime(Sys.time(), start, units = "s"))
  
  # Now do anticlustering
  binary_categories <- categories_to_binary(data)
  dists <- as.matrix(dist(binary_categories))
  # Perform anticlustering using default exchange method with default objective ("Euclidean diversity")
  start <- Sys.time()
  ANTICLUST <- anticlustering(
    dists,
    K = K
  )
  ANTICLUST_t <- as.numeric(difftime(Sys.time(), start, units = "s"))
  
  must_link <- sample(N, replace = TRUE)
  while (any(table(must_link) > (N/K))) { # ensure that no must-link group is larger than the group size
    must_link <- sample(N, replace = TRUE)
  }
  
  # Perform must-link anticlustering using 2PML
  start <- Sys.time()
  ANTICLUST2 <- anticlustering(
    dists,
    K = K,
    must_link = must_link,
    method = "2PML"
  )
  ANTICLUST2_t <- as.numeric(difftime(Sys.time(), start, units = "s"))
  
  percent_diversity <- diversity_objective(dists, ANTICLUST2) / diversity_objective(dists, ANTICLUST)
  
  ## Perform Propensity Score Batch Effect
  # ONLY USE FOR K <= 4 (so here, not for K = 10)
  USE_PS <- K <= 4
  if (USE_PS) {
    start <- Sys.time()
    PS <- my_ps_batch_effect(binary_categories, K = K, SEED = seeds_ps[i])
    PS_time <- as.numeric(difftime(Sys.time(), start, units = "s"))
  }

  
  # OSAT paper uses p values to quantify discrepancy in each category between batches, which makes sense
  # We do that as well.
  pvalues_osat <- rep(NA, 5)
  pvalues_anticlust <- rep(NA, 5)
  pvalues_anticlust2 <- rep(NA, 5)
  pvalues_ps <- rep(NA, 5)
  pvalues_osat[1:M] <- sapply(1:M, function(x) chisq.test(table(data[, x], OSAT))$p.value)
  pvalues_anticlust[1:M] <- sapply(1:M, function(x) chisq.test(table(data[, x], ANTICLUST))$p.value)
  pvalues_anticlust2[1:M] <- sapply(1:M, function(x) chisq.test(table(data[, x], ANTICLUST2))$p.value)
  if (USE_PS) {
    pvalues_ps[1:M] <- sapply(1:M, function(x) chisq.test(table(data[, x], PS))$p.value)
  }
  
  results_file <- here("Simulation", "results.csv")
  results_file_exists <- file.exists(results_file)
  
  df <- data.frame(
    N = N, 
    M = M, 
    K = K,
    P = P,
    named_1row_matrix(pvalues_osat, "p_osat"),
    named_1row_matrix(pvalues_anticlust, "p_anticlust"),
    named_1row_matrix(pvalues_anticlust2, "p_anticlust_c"),
    named_1row_matrix(pvalues_ps, "p_ps"),
    OSAT_t = OSAT_t,
    ANTICLUST_t = ANTICLUST_t,
    ANTICLUST_t_c = ANTICLUST2_t,
    PS_t = ifelse(USE_PS, PS_time, NA),
    diversity_unconstrained = diversity_objective(dists, ANTICLUST),
    diversity_constrained = diversity_objective(dists, ANTICLUST2)
  )
  
  if (!results_file_exists) {
    file.create(results_file)
  }
  
  write.table(
    write.table(
      df, 
      file = results_file, 
      append = results_file_exists, 
      col.names = !results_file_exists,
      row.names = FALSE, 
      quote = FALSE,
      sep = ";"
    )
  )
}
