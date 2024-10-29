# oast_vs_anticlust.R

library(anticlust) # I used version 0.8.7, which has must_link argument for `anticlustering()`
library(OSAT) # I got version OSAT_1.52.0 from Bioconductor
library(here) # for managing relative paths in the project

sessionInfo()

# Wrapper function for OSAT method that can be used in the context 
# of the simulation. It returns a vector of group assignments (as does anticlustering())
my_osat <- function(data, K, nSim = 5000) {
  data <- as.data.frame(data)
  gs <- setup.sample(data, optimal = colnames(data))
  data$ID <- 1:nrow(data)
  N <- nrow(data)
  
  perfect_split_possible <- N %% K == 0
  
  myChip <- new("BeadChip", nRows=N/K, nColumns=1, byrow=FALSE,
                comment="Dummy chip.")
  chip_plate <- myPlate <- new("BeadPlate", chip=myChip, nRows = 1L, nColumns=1L)
  # Dummy plate with just 1 chip that already has N rows.
  gc <- setup.container(chip_plate, K, batch='plates')
  
  # DO OSAT
  gSetup <- create.optimized.setup(sample = gs, container = gc, nSim = nSim, fun = "optimal.shuffle")
  # Merge OSAT groupings to initial data frame
  df <- merge(data, gSetup@data$link[, c("OrigRow", "plates")], by.x = "ID", by.y = "OrigRow")
  df <- df[order(df$ID), ]
  
  if (perfect_split_possible) {
    tab <- table(df$plates)
    stopifnot(all(tab[1] == tab)) # ensure that OSAT returns equal sized groups
    stopifnot(length(tab) == K) # ensure that the correct number of clusters is returned
  }
  df$plates
}

# Also create wrapper function for propensity score method
# This uses the R files of the authors, which are available from Github, but are not redistributed here
# because they lack a license. On my local computer, I just downloaded the files and
# put them in a subdirectory "/PS-Batch-Effect/". 
# Here, I check if the files if such a directory exists, and the files are downloaded from Github if not.
files_locally_available <- dir.exists(here::here("Simulation_vs_OSAT", "PS-Batch-Effect"))
if (files_locally_available) {
  source(here::here("Simulation_vs_OSAT", "./PS-Batch-Effect/Almost.Optimal.R"))
  source(here::here("Simulation_vs_OSAT", "./PS-Batch-Effect/R.Allocation.R"))
  source(here::here("Simulation_vs_OSAT", "./PS-Batch-Effect/SR.Allocation.R"))
} else {
  source("https://raw.githubusercontent.com/carryp/PS-Batch-Effect/refs/heads/main/R/Almost.Optimal.R")
  source("https://raw.githubusercontent.com/carryp/PS-Batch-Effect/refs/heads/main/R/R.Allocation.R")
  source("https://raw.githubusercontent.com/carryp/PS-Batch-Effect/refs/heads/main/R/SR.Allocation.R")
}

# In this approach I do not use a stratification variable, which their implementation 
# allows, but only if the stratification variable has 2 levels. In the simulation
# there is (a) no obvious variable which serves for stratification purposes 
# and (b) we use more than 2 levels for our categorical variables. 

my_ps_batch_effect <- function(data, K, nSim = 5000, SEED) {
  data <- data.frame(data) # binary_categories returns a matrix...
  N  <- nrow(data)
  M <- ncol(data)
  group_sizes <- rep(0, 4)
  group_sizes[1:K] <- N/K
  output1 <- R.Allocation(
    Pheno.Data = data,
    N.Iter = nSim,
    Initial.Seed = SEED,
    B1 = group_sizes[1], B2 = group_sizes[2], B3 = group_sizes[3], B4 = group_sizes[4]
  )
  output2 <- Almost.Optimal(
    Covariates = colnames(data),
    Pheno.Dataset = data,
    Randomization.Dataset = output1
  )
  groups <- as.numeric(as.factor(output2$Batch.Assignment.Av))
  # some checking of the solution:
  tab <- table(groups)
  stopifnot(all(tab[1] == tab)) # ensure that PS returns equal sized groups
  stopifnot(length(tab) == K) # ensure that the correct number of clusters is returned
  groups
}

# for converting the p values to columns in data frame
named_1row_matrix <- function(x, prefix) {
  t(matrix(x, ncol = 1, dimnames = list(paste0(prefix, 1:length(x)))))
}


generate_categorical_data <- function(N, M, P) {
  data <- data.frame(matrix(sample(P, replace = TRUE, size = N*M), ncol = M))
  as.data.frame(lapply(data, as.factor))
}

# set.seed(123) # set seed for final simulation 
nsim <- 500

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
  
  # Perform anticlustering using default exchange method + constraints
  start <- Sys.time()
  ANTICLUST2 <- anticlustering(
    dists,
    K = K,
    must_link = must_link
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

  
  # OSAT vignette uses p values to quantify discrepancy in each category between batches
  # (I would also use p values in the simulation in our paper I guess)
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
  
  results_file <- here("Simulation_vs_OSAT", "results.csv")
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


