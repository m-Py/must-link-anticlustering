
# Author: Martin Papenberg
# Year: 2025

# This file implements the analysis of the simulation study, making the results
# in the paper reproducible.

library(here)
library(tidyverse)
library(santoku)

tt <- read.csv(here("Simulation", "results.csv"), sep = ";")
tt$ID <- 1:nrow(tt) # ID that identifies a single simulation run / data set in the data. Useful when data is transformed into long format.
nrow(tt)
p_values <- tt[, grepl("p_", colnames(tt))]


# global results OSAT vs. anticlust. How often is p value of anticlust better (i.e. larger)?
anticlust_better <- apply(p_values, 1, function(x) x[paste0("p_anticlust", 1:5)] > x[paste0("p_osat", 1:5)])
osat_better <- apply(p_values, 1, function(x) x[paste0("p_anticlust", 1:5)] < x[paste0("p_osat", 1:5)])
same <- apply(p_values, 1, function(x) x[paste0("p_anticlust", 1:5)] == x[paste0("p_osat", 1:5)])

# How often was anticlust better than OSAT (percentage):
mean(anticlust_better, na.rm = TRUE) |> round(4) * 100
# How often was OSAT better than anticlust (percentage):
mean(osat_better, na.rm = TRUE) |> round(4) * 100
# How often was OSAT better than anticlust (frequency):
sum(osat_better, na.rm = TRUE)
# Total number of comparisons between OSAT and anticlust:
length(osat_better[!is.na(osat_better)]) 
# How often was OSAT equal to anticlust (percentage):
mean(same, na.rm = TRUE) |> round(4) * 100

# global results PS vs. anticlust. How often is p value of anticlust better (i.e. larger)?
anticlust_better <- apply(p_values, 1, function(x) x[paste0("p_anticlust", 1:5)] > x[paste0("p_ps", 1:5)])
ps_better <- apply(p_values, 1, function(x) x[paste0("p_anticlust", 1:5)] < x[paste0("p_ps", 1:5)])
same <- apply(p_values, 1, function(x) x[paste0("p_anticlust", 1:5)] == x[paste0("p_ps", 1:5)])

# How often was anticlust better than PSBA (percentage):
mean(anticlust_better, na.rm = TRUE) |> round(4) * 100
# How often was PSBA better than anticlust (percentage):
mean(ps_better, na.rm = TRUE) |> round(4) * 100
# How often was PSBA better than anticlust (frequency):
sum(ps_better, na.rm = TRUE) 
# Total number of comparisons between PSBA and anticlust (it's lower than for OSAT because PSBA was not used for K = 10)
length(ps_better[!is.na(ps_better)])
# How often was PSBA equal to anticlust (percentage):
mean(same, na.rm = TRUE) |> round(2)


### Average relative objective of constrained assignment:
mean(tt$diversity_constrained / tt$diversity_unconstrained) |> round(3)

# How often is balance not reduced by constraints:
mean(apply(p_values, 1, function(x) x[paste0("p_anticlust", 1:5)] <= x[paste0("p_anticlust_c", 1:5)]), na.rm = TRUE) |> round(2)

# Average run times by method:
runtimes <- colMeans(tt[, grepl("_t", colnames(tt))], na.rm = TRUE) |> round(2)
# run times relative to anticlustering run time:
(runtimes / runtimes["ANTICLUST_t"]) |> round()

## Some more sophisticated analyes using grouping via dplyr:

# Create long format data:
dfl <- tt |> 
  select(ID, N, M, K, P, starts_with("p")) |> 
  pivot_longer(
    cols = starts_with("p_")
  ) |> 
  filter(!is.na(value))

dfl$Method <- "Anticlustering"
dfl$Method[grepl("osat", dfl$name)] <- "OSAT"
dfl$Method[grepl("p_anticlust_c", dfl$name)] <- "Must-Link Anticlustering"
dfl$Method[grepl("p_ps", dfl$name)] <- "Propensity Score Batch Effect"

## Analyze data by variables

## This Figure is in the paper, which illustrates balance by number of batches (K) and number of variables.

facets <- c(
  `2` = "K = 2",
  `4` = "K = 4",
  `10` = "K = 10"
)

dfl |> 
  group_by(Method, M, K) |> 
  summarise(`p value` = mean(value)) |> 
  ggplot(aes(y = `p value`, x = M, color = Method)) +
  geom_point(aes(color = Method, shape = Method)) + 
  geom_line(aes(color = Method, linetype = Method)) +
  facet_grid(cols = vars(`K`), labeller = as_labeller(facets))+
  theme_bw(base_size = 14) +
  xlab("Number of variables") +
  ylab("Average p value")


## For the supplementary materials, we also do this plot for the other parameters that varied in the simulation:
# Sample size (N), Number of categories per variable (P)

# N (categorized)

dfl$N_category <- santoku::chop(dfl$N, breaks = c(100, 300))

dfl |> 
  group_by(Method, N_category = as.numeric(N_category), K) |> # I hate it so much that ggplot does not draw geom_line for ordered factors
  summarise(`p value` = mean(value)) |> 
  ggplot(aes(y = `p value`, x = N_category, color = Method), group = 1) +
  geom_point(aes(color = Method, shape = Method)) + 
  geom_line(aes(color = Method, linetype = Method)) +
  facet_grid(cols = vars(`K`), labeller = as_labeller(facets))+
  theme_bw(base_size = 14) +
  xlab("Sample size (categorized)") +
  ylab("Average p value") +
  scale_x_continuous(breaks = 1:3, labels = levels(dfl$N_category)) +
  expand_limits(x=c(0.5, 3.5))

# P (number of categories per variable)

dfl |> 
  group_by(Method, P, K) |> 
  summarise(`p value` = mean(value)) |> 
  ggplot(aes(y = `p value`, x = P, color = Method)) +
  geom_point(aes(color = Method, shape = Method)) + 
  geom_line(aes(color = Method, linetype = Method)) +
  facet_grid(cols = vars(`K`), labeller = as_labeller(facets))+
  theme_bw(base_size = 14) +
  xlab("Number of categories per variable") +
  ylab("Average p value")


## Expected number of must-link partners in the simulation:
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
