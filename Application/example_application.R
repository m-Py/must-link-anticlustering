
library(anticlust) 
# currently requires a development version, available via
# remotes::install_github("m-Py/anticlust", ref = "must-link-new")

tt <- read.csv("./Application/sample_data.csv")
N <- nrow(tt)
tt <- tt[order(sample(N)), ]
table(table(tt$patientID)) # good!
summary(tt)

# generate input data for k-means anticlustering
numeric_features <- tt[, c("age", "BMI")]
categorical_features <- tt[, c("race", "ethnicity", "clinical_site", 
       "anatomical_site", "disease_stage", "cycle_phase", "sample_type", "lesion_type")]

kmeans_input <- cbind(numeric_features, categories_to_binary(categorical_features)) 
# standardization might help, i.e. use scale(kmeans_input)

distances <- dist(kmeans_input)^2 # squared distance matrix for k-means anticlustering

cluster_sizes_max <- table(rep_len(1:16, N))
start <- Sys.time()
groups_must_link <- anticlustering(
  distances, 
  K = cluster_sizes_max, 
  must_link = tt$patientID, 
  method= "local-maximum", # method = "local-maximum" may be better but takes longer (about 7min on my computer)
  repetitions = 10
)
Sys.time() - start 

# Compare objectives:
diversity_objective(distances, groups_must_link) # 597304 was the best one I found using 1000 repetitions

# verify that must-link constraints are met:
ID <- tt$patientID
same <- as.numeric(names(table(ID)[table(ID) > 1]))
for (i in same) {
  stopifnot(all(groups_must_link[ID == i] == groups_must_link[ID == i][1]))
}

# check out in data frame
tt$Anticlusters <- groups_must_link

