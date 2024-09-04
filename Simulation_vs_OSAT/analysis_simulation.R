
tt <- read.csv("./Simulation_vs_OSAT/results_N=4543_save.csv", sep = ";")
nrow(tt)
p_values <- tt[, grepl("p_", colnames(tt))]

# global results OSAT vs. anticlust. How often is p value of anticlust better (i.e. larger)?
anticlust_better <- apply(p_values, 1, function(x) x[1:5] < x[6:10])
osat_better <- apply(p_values, 1, function(x) x[1:5] > x[6:10])
same <- apply(p_values, 1, function(x) x[1:5] == x[6:10])

mean(anticlust_better, na.rm = TRUE) |> round(2)
mean(osat_better, na.rm = TRUE) |> round(2)
mean(same, na.rm = TRUE) |> round(2)

# osat vs. constrained anticlust 
anticlust_better <- apply(p_values, 1, function(x) x[1:5] < x[11:15])
osat_better <- apply(p_values, 1, function(x) x[1:5] > x[11:15])
same <- apply(p_values, 1, function(x) x[1:5] == x[11:15])

mean(anticlust_better, na.rm = TRUE) |> round(2)
mean(osat_better, na.rm = TRUE) |> round(2)
mean(same, na.rm = TRUE) |> round(2)

# Average relative objective of constrained assignment:
mean(tt$diversity_constrained / tt$diversity_unconstrained) |> round(3)


# How often is balance not reduced by constraints:
mean(apply(p_values, 1, function(x) x[6:10] <= x[11:15]), na.rm = TRUE) |> round(2)

# Run times. We do not report those in the paper but the difference is actually
# quite substantial (anticlust is implemented in C and OSAT in R. I think this 
# explains the difference. I even think that anticlustering should actually be 
# the computationally more expensive algorithm - but I'm not entirely sure)
mean(tt$OSAT_t)
mean(tt$ANTICLUST_t)
mean(tt$ANTICLUST_t_c)



## Expected number of must-link partners:

Ns <- tt$N

tables <- lapply(lapply(Ns, function(x) sample(x, replace = TRUE)), function(x) table(table(x)))
# make each member the same length:
largest_group <- max(lengths(tables))

new_matrix <- matrix(0, nrow = length(tables), ncol = largest_group )

for (i in seq_along(tables)) {
  new_matrix[i, 1:length(tables[[i]])] <- tables[[i]]
}

ff <- t(apply(new_matrix, 1, function(x) x / sum(x)))
results <- colMeans(ff)
results[4] <- sum(results[4:6])
results <- results[1:4]

names(results) <- c("1", "2", "3", "4 or more")
results |> round(2)
#         1         2         3 4 or more 
#      0.58      0.29      0.10      0.03 
# These results are quite stable

# example application has:
# 0.523560209 0.246073298 0.068062827 0.125654450 0.015706806 0.015706806 0.005235602
