
library(here)
library(tidyverse)
library(santoku)

tt <- read.csv(here("Simulation_vs_OSAT", "results.csv"), sep = ";")
tt$ID <- 1:nrow(tt)
nrow(tt)
p_values <- tt[, grepl("p_", colnames(tt))]


# global results OSAT vs. anticlust. How often is p value of anticlust better (i.e. larger)?
anticlust_better <- apply(p_values, 1, function(x) x[1:5] < x[6:10])
osat_better <- apply(p_values, 1, function(x) x[1:5] > x[6:10])
same <- apply(p_values, 1, function(x) x[1:5] == x[6:10])

mean(anticlust_better, na.rm = TRUE) |> round(2)
mean(osat_better, na.rm = TRUE) |> round(2)
mean(same, na.rm = TRUE) |> round(2)

# global results PS vs. anticlust. How often is p value of anticlust better (i.e. larger)?
anticlust_better <- apply(p_values, 1, function(x) x[1:5] < x[16:20])
ps_better <- apply(p_values, 1, function(x) x[1:5] > x[16:20])
same <- apply(p_values, 1, function(x) x[1:5] == x[16:20])

mean(anticlust_better, na.rm = TRUE) |> round(2)
mean(ps_better, na.rm = TRUE) |> round(2)
mean(same, na.rm = TRUE) |> round(2)

# PS method seems to be better than OSAT, but still outperformed by anticlust

# Average relative objective of constrained assignment:
mean(tt$diversity_constrained / tt$diversity_unconstrained) |> round(3)


# How often is balance not reduced by constraints:
mean(apply(p_values, 1, function(x) x[6:10] <= x[11:15]), na.rm = TRUE) |> round(2)

# Run times. We also should report those in the paper
colMeans(tt[, grepl("_t", colnames(tt))], na.rm = TRUE) |> round(2)

## Some more sophisticated analyes:

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

# N - OSAT gets better with increasing N (but always worse than anticlustering)
dfl |> 
  group_by(Method, N_category = santoku::chop(N, breaks = c(100, 300))) |> 
  summarise(mean(value))
# M - interesting, OSAT gets worse with increasing M
dfl |> 
  group_by(Method, M) |> 
  summarise(mean(value))

# K (no interaction here, just anticlust always better)
dfl |> 
  group_by(Method, K) |> 
  summarise(mean(value))

# P (OSAT gets worse with increasing P)
dfl |> 
  group_by(Method, P) |> 
  summarise(mean(value))

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


## Expected number of must-link partners:

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
