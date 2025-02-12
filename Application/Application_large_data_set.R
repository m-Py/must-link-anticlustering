
# Author: Martin Papenberg
# Year: 2025


library(tableone)
library(anticlust) # at least requires anticlust version 0.8.9 
library(here)

dataset <- read.csv(here("Application", "synthetic_dataset_10x_20250211.csv"))
data_ind <- dataset[!duplicated(dataset$patient_id), ] # only unique patients / individuals

table(table(dataset$patient_id))

vars_categorical <- c("endo", "site", "phase", "stage")
vars_numeric <- "age"

CreateTableOne(data = dataset, vars = c(vars_categorical, vars_numeric))

categorical_input <- categories_to_binary(dataset[, vars_categorical])
input <- cbind(categorical_input, dataset$age)
# squared Euclidean distance on binary coded features

K <- 20

start <- Sys.time()
dataset$BatchAnticlust <- anticlustering(
  input, 
  K = K
)
Sys.time() - start 

# Use specialized function for large data sets: 
# (optimizes k-means criterion and not diversity, does not incorporate must-link constraints)
start <- Sys.time()
dataset$BatchAnticlust_Fast <- fast_anticlustering(
  input, 
  K = K
)
Sys.time() - start 

start <- Sys.time()
dataset$BatchAnticlustML <- anticlustering(
  input, 
  K = K, 
  method = "2PML",
  must_link = dataset$patient_id, 
  repetitions = 1 # 1 usually leads to significant differences on the features endo and stage; 100 had good results for me but took 13min
)
Sys.time() - start 


print(tableone::CreateTableOne(data = dataset, vars = c(vars_categorical, vars_numeric), strata = "BatchAnticlust"), smd = TRUE)
print(tableone::CreateTableOne(data = dataset, vars = c(vars_categorical, vars_numeric), strata = "BatchAnticlustML"), smd = TRUE)
print(tableone::CreateTableOne(data = dataset, vars = c(vars_categorical, vars_numeric), strata = "BatchAnticlust_Fast"), smd = TRUE)
