
# Author: Martin Papenberg
# Year: 2025

# This script reproduces the Batch Assignments for the synthetic data set presented in the paper
# (1) Unrestricted Anticlustering
# (2) OSAT
# (3) Must-link anticlustering

# The results/assignments are written to a file "assignments_synthetic_dataset.csv"
# and are further processed in the file 01_Analysis.R to check balance among batches

# The script is reprocucible via a seed  (100)
set.seed(100)

library(OSAT)
library(tableone)
library(anticlust) # at least requires anticlust version 0.8.9 
library(here)
library(dplyr)
source(here("Simulation", "functions.R")) # I have a wrapper for the OSAT method in the Simulation directory


dataset <- read.csv(here("Application", "synthetic_dataset_20250114.csv"))
data_ind <- dataset[!duplicated(dataset$patient_id), ] # only unique patients / individuals

table(table(dataset$patient_id))
table(table(data_ind$patient_id)) # 139 individuals

CreateTableOne(data = dataset, vars = c("endo", "site", "phase", "stage"))
CreateTableOne(data = data_ind, vars = c("endo", "site", "phase", "stage"))

input <- dist(categories_to_binary(dataset[, 3:6]))^2 # squared Euclidean distance on binary coded features

reps <- 500 # 500 restarts for anticlustering methods
K <- 20

start <- Sys.time()
dataset$BatchAnticlust <- anticlustering(
  input, 
  K = K, 
  method = "local-maximum", 
  repetitions = reps
)
Sys.time() - start 

start <- Sys.time()
dataset$BatchAnticlustML <- anticlustering(
  input, 
  K = K, 
  method = "2PML", 
  repetitions = reps, 
  must_link = dataset$patient_id
)
Sys.time() - start 

start <- Sys.time()
dataset$BatchOSAT <- my_osat(dataset[, 3:6], K = 20, nSim = 5000)
Sys.time() - start 

# check balance / diversity (diversity is being optimized via the algorithms)
diversity_objective(input, dataset$BatchAnticlust) # 8244
diversity_objective(input, dataset$BatchOSAT) # 8204
diversity_objective(input, dataset$BatchAnticlustML) # 8182 with 2PML method (250x Phase 1, 250x Phase 2)

## NOW: Encode if a sample is a duplicated within a batch (this information is used to assess balance on individuals' rather than samples' level)

# These `by()` operations are only used (below) to verify the output of the dplyr operation (which is much more useful for further processing)
# sum(unlist(by(data = dataset$patient_id, INDICES = dataset$BatchAnticlust, FUN = duplicated)))
# sum(unlist(by(data = dataset$patient_id, INDICES = dataset$BatchAnticlustML, FUN = duplicated)))
# sum(unlist(by(data = dataset$patient_id, INDICES = dataset$BatchOSAT, FUN = duplicated)))

dataset <- dataset |> 
  group_by(BatchAnticlust) |> 
  reframe(DUPLICATED_ANTICLUST = duplicated(patient_id), patient_id = patient_id, sample_num = sample_num) |> 
  full_join(dataset, by = c("patient_id", "BatchAnticlust", "sample_num"))
stopifnot(
  sum(dataset$DUPLICATED_ANTICLUST) == sum(unlist(by(data = dataset$patient_id, INDICES = dataset$BatchAnticlust, FUN = duplicated)))
)

dataset <- dataset |> 
  group_by(BatchAnticlustML) |> 
  reframe(DUPLICATED_ANTICLUST_ML = duplicated(patient_id), patient_id = patient_id, sample_num = sample_num) |> 
  full_join(dataset, by = c("patient_id", "BatchAnticlustML", "sample_num"))
stopifnot(
  sum(dataset$DUPLICATED_ANTICLUST_ML) == sum(unlist(by(data = dataset$patient_id, INDICES = dataset$BatchAnticlustML, FUN = duplicated)))
)

dataset <- dataset |> 
  group_by(BatchOSAT) |> 
  reframe(DUPLICATED_OSAT = duplicated(patient_id), patient_id = patient_id, sample_num = sample_num) |> 
  full_join(dataset, by = c("patient_id", "BatchOSAT", "sample_num"))
stopifnot(
  sum(dataset$DUPLICATED_OSAT) == sum(unlist(by(data = dataset$patient_id, INDICES = dataset$BatchOSAT, FUN = duplicated)))
)

# Write data to file
write.table(
  dataset, 
  file = here("Application", "assignments_synthetic_dataset.csv"),
  quote = TRUE, 
  row.names = FALSE, 
  sep = ";"
)
