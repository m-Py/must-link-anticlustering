
# Author: Martin Papenberg
# Year: 2025

## Some functions that are used in the simulation


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
files_locally_available <- dir.exists(here::here("Simulation", "PS-Batch-Effect"))
if (files_locally_available) {
  source(here::here("Simulation", "./PS-Batch-Effect/Almost.Optimal.R"))
  source(here::here("Simulation", "./PS-Batch-Effect/R.Allocation.R"))
  source(here::here("Simulation", "./PS-Batch-Effect/SR.Allocation.R"))
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
