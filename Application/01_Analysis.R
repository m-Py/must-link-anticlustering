
# Author: Martin Papenberg
# Year: 2025

# This script analyses the batch assignments generated via script 00_Assignments.R
# Uses tableone package to assess balance on the levels of (a) samples, (b) individuals
library(tableone)
library(here)

# Read data set that includes the assignments, as well as information on duplication of an individual within samples
dataset <- read.csv(here("Application", "assignments_synthetic_dataset.csv"), sep = ";")


# Balance for the entire data set including duplicate samples
save_t1 <- function(t1, file) {
  write.csv(
    print(
      t1, 
      quote = TRUE,
      noSpaces = TRUE,
      printToggle = FALSE
    ), file = file
  )
}

## Encode per batch if individuals in it are unique to the Batch
# (a) OSAT
sapply(1:20, function(i) !any(dataset$patient_id[dataset$BatchOSAT == i] %in% dataset$patient_id[dataset$BatchOSAT != i]))
# (b) unrestricted anticlustering
sapply(1:20, function(i) !any(dataset$patient_id[dataset$BatchAnticlust == i] %in% dataset$patient_id[dataset$BatchAnticlust != i]))
# (c) must-link constrained anticlustering
sapply(1:20, function(i) !any(dataset$patient_id[dataset$BatchAnticlustML == i] %in% dataset$patient_id[dataset$BatchAnticlustML != i]))


## <- We insert this info manually into the table because aggregating it programatically with the tableone output is quite hard

## Now Generate Tableones to view balance among batches
vars <- c("endo", "site", "phase", "stage")

# Level of samples 
save_t1(
  tableone::CreateCatTable(data = dataset, vars = vars, strata = "BatchAnticlust"),
  file = here("Application", paste0(format(Sys.time(), "%Y %m %d"), " Tableone Unrestricted Anticlustering Assignment - All Samples.csv"))
)
save_t1(
  tableone::CreateCatTable(data = dataset, vars = vars, strata = "BatchOSAT"),
  file = here("Application", paste0(format(Sys.time(), "%Y %m %d"), " Tableone OSAT Assignment - All Samples.csv"))
)
save_t1(
  tableone::CreateCatTable(data = dataset, vars = vars, strata = "BatchAnticlustML"),
  file = here("Application/", paste0(format(Sys.time(), "%Y %m %d"), " Tableone Must-Link Assignment - All Samples.csv"))
)

# Level of individuals
save_t1(
  tableone::CreateCatTable(data = dataset[!dataset$DUPLICATED_ANTICLUST, ], vars = vars, strata = "BatchAnticlust"),
  file = here("Application", paste0(format(Sys.time(), "%Y %m %d"), " Tableone Unrestricted Anticlustering Assignment - Individuals.csv"))
)

save_t1(
  tableone::CreateCatTable(data = dataset[!dataset$DUPLICATED_OSAT, ], vars = vars, strata = "BatchOSAT"),
  file = here("Application", paste0(format(Sys.time(), "%Y %m %d"), " Tableone OSAT Assignment - Individuals.csv"))
)

save_t1(
  tableone::CreateCatTable(data = dataset[!dataset$DUPLICATED_ANTICLUST_ML, ], vars = vars, strata = "BatchAnticlustML"),
  file = here("Application", paste0(format(Sys.time(), "%Y %m %d"), " Tableone Must-Link Assignment - Individuals.csv"))
)

sum(!dataset$DUPLICATED_ANTICLUST_ML) # this is actually the number of individuals; each individual is assigned to only one batch using ML anticlustering
# not true for the other methods:
sum(!dataset$DUPLICATED_ANTICLUST)
sum(!dataset$DUPLICATED_OSAT)
