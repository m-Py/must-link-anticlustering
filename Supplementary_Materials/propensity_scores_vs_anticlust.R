
library(anticlust)
library(table1)

# Unfortunately, their code / data is not available with a license, so I download their 
# data file in this script instead of re-distributing them here. 
# When their repository is offline, this may no longer work.

# Phenotype file in their repository:
Pheno <- readRDS(file = url("https://github.com/carryp/PS-Batch-Effect/raw/refs/heads/main/Data/Pheno.GSE122288.Clean.Rds"))
# Data originally from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122288 (they did some preprocessing)
Pheno$neonate.gender <- factor(Pheno$neonate.gender)
# The implementation of the algorithm is currently available in https://github.com/carryp/PS-Batch-Effect
# (specifically https://raw.githubusercontent.com/carryp/PS-Batch-Effect/refs/heads/main/Vignette.Rmd)

# I do not reproduce their code here because it lacks a license. But I was exactly able to 
# reproduce their results using the code in the vignette.

# In the paper I reproduce the numbers they have in their rendered vignette (retrieved on 2024 09 25)

## This does the same with anticlustering

# compare to anticlustering
library(anticlust)
input_anticlustering <- cbind(
  scale(Pheno[, c("birth.weight", "gestational.age")]),
  categories_to_binary(Pheno$neonate.gender)
)

set.seed(12)
Pheno$ANTICLUSTERS <- anticlustering(
  input_anticlustering,
  K = 3
)

# Their method uses a stratification split on gender, which we could also do when 
# we pass Pheno$neonate.gender to an additional argument `categories` rather
# than using it as part of the objective (in the first argument). both works here.

# do some re-labeling because their method has group 3 as the largest group, anticlust has group 1
Pheno$ANTICLUSTERS[Pheno$ANTICLUSTERS == 1] <- 4
Pheno$ANTICLUSTERS[Pheno$ANTICLUSTERS == 3] <- 1
Pheno$ANTICLUSTERS[Pheno$ANTICLUSTERS == 4] <- 3
Pheno$ANTICLUSTERS <- factor(paste0("B", Pheno$ANTICLUSTERS))

tabgeo <- table1(~ neonate.gender + birth.weight + gestational.age | ANTICLUSTERS, data = Pheno)

write.table(
  as.data.frame(tabgeo), 
  file = "Results_GSE122288_Anticlust.csv",
  sep = ";"
)

