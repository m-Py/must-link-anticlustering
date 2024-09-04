# oast_vs_anticlust.R

library(anticlust) # I used 0.8.6.9999, development version, which has must_link argument for `anticlustering()`
library(OSAT) # I got version OSAT_1.52.0 from Bioconductor

sessionInfo()

# data = data frame with categorical variables that should be balanced
# 
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
  }
  df$plates
}

# for converting the p values to columns in data frame
named_1row_matrix <- function(x, prefix) {
  t(matrix(x, ncol = 1, dimnames = list(paste0(prefix, 1:length(x)))))
}

# set.seed(123) # set seed for final simulation 
nsim <- 10000

generate_categorical_data <- function(N, M, P) {
  data <- data.frame(matrix(sample(P, replace = TRUE, size = N*M), ncol = M))
  as.data.frame(lapply(data, as.factor))
}


start_all <- Sys.time()
for (i in 1:nsim) {
  cat("Working iteration", i, "\n")
  cat("Time passed", as.character(round(as.numeric(difftime(Sys.time(), start_all, units = "s")), 2)), "s\n")
  K <- sample(c(2, 5, 10), size = 1)
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
  
  # OSAT vignette uses p values to quantify discrepancy in each category between batches
  # (I would also use p values in the simulation in our paper I guess)
  pvalues_osat <- rep(NA, 5)
  pvalues_anticlust <- rep(NA, 5)
  pvalues_anticlust2 <- rep(NA, 5)
  pvalues_osat[1:M] <- sapply(1:M, function(x) chisq.test(table(data[, x], OSAT))$p.value)
  pvalues_anticlust[1:M] <- sapply(1:M, function(x) chisq.test(table(data[, x], ANTICLUST))$p.value)
  pvalues_anticlust2[1:M] <- sapply(1:M, function(x) chisq.test(table(data[, x], ANTICLUST2))$p.value)
  
  results_file <- "results.csv"
  results_file_exists <- file.exists(results_file)
  
  df <- data.frame(
    N = N, 
    M = M, 
    K = K,
    P = P,
    named_1row_matrix(pvalues_osat, "p_osat"),
    named_1row_matrix(pvalues_anticlust, "p_anticlust"),
    named_1row_matrix(pvalues_anticlust2, "p_anticlust_c"),
    OSAT_t = OSAT_t,
    ANTICLUST_t = ANTICLUST_t,
    ANTICLUST_t_c = ANTICLUST2_t,
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


