
# Author: Martin Papenberg
# Year: 2025

# Some analyses on optimal anticlustering with must-link feature. Also see the supplementary materials.
# (Supplementary_Materials/supplemetary_materials.pdf)

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

ti <- read.csv(here("Running_Time_Optimal_Algorithm", "results_runningtime.csv"), sep = ";")

# we need comparator that allows for numeric imprecision of double, to test if heuristic solutions were optimal 
# (as in file running_times.R)

"%==%" <- function(x, y) abs(x - y) < .0000000001
"%>%" <- function(x, y) (x - y) > (.0000000001)

mean(ti$diversity_optimal %>% ti$diversity_heuristic, na.rm = TRUE) 
ti$ID <- 1:nrow(ti)

dfl <- ti |> 
  select(ID, N, K, time_heuristic, time_optimal, time_limit_exceeded) |> 
  pivot_longer(
    cols = c(time_heuristic, time_optimal),
    names_to = "Algorithm", 
    values_to = "Runtime",
    names_prefix = "time_")

dfl$K <- factor(dfl$K)
df <- dfl |> 
  filter(!time_limit_exceeded) |> 
  group_by(N, K, Algorithm) |> 
  summarise(runtime = mean(Runtime))

df |> 
  group_by(K, Algorithm) |> 
  summarise(max_mean_runtime = max(runtime), N_max = max(N)) |> 
  arrange(desc(Algorithm), K)

options(scipen = 999, digits = 2)
df |>
  ggplot(aes(x = N, y = log(runtime))) +
  geom_point(aes(color = K, shape = K)) +
  geom_line(aes(linetype = K, color = K)) +
  theme_bw(base_size = 14) +
  facet_grid(cols = vars(Algorithm), scales = "free") +
  ylab("log(Run time)") + 
  theme(legend.position = "top")+ 
  scale_y_continuous(
    sec.axis = sec_axis(transform = ~ exp(.x),
                        name = "Run time (s)", 
                        breaks = exp(c(-6, -3, 0, 3, 6)))
  ) 

