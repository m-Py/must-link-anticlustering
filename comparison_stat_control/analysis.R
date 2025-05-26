
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)

# read all data files
files <- list.files(here("comparison_stat_control", "data"), full.names = TRUE)

df <- do.call("rbind", lapply(files, read.csv, sep =";"))

head(df)

mean(df[df$treatment_effect > 0, "p_anticlust_control"] <  df[df$treatment_effect > 0, "p_rnd_control"])

# apply(subset(df, select = K:adjust_for_covariate), 2, table)

ldf <- df |>
  select(ID, N, M, K, scale_batch_effect, treatment_effect, adjust_for_covariate, starts_with("p_")) |>
  pivot_longer(
    cols = starts_with("p_"),
    names_prefix = "p_",
    names_to = "Method",
    values_to = "pvalue"
  )

ldf$Control <- as.numeric(!grepl(pattern = "_no_control", ldf$Method))

ldf |>
  filter(treatment_effect > 0) |>
  group_by(Method, K) |>
  summarize(
    Power = mean(pvalue <= .05)
  ) |> arrange(grepl(pattern = "_no_control", Method)) |> 
  pivot_wider(names_from = Method, values_from = Power)


ldf |>
  filter(treatment_effect > 0) |>
  group_by(adjust_for_covariate, K, Method) |>
  summarize(
    Power = mean(pvalue <= .05)
  ) |> arrange(grepl(pattern = "_no_control", Method)) |> 
  pivot_wider(names_from = Method, values_from = Power)


# We have a function that compute confidence interval for power. Use Method 3 from
# Newcombe, R. G. (1998). Two‐sided confidence intervals for the single proportion:
# comparison of seven methods. Statistics in medicine, 17(8), 857-872.

ci_one_prop <- function(ci = 95, r, n, return_string = TRUE, decimals = 3) {
  z <- get.z.score(ci)
  p <- r/n
  lower <- (2 * n * p + z^2 - (z * sqrt(z^2 + 4 * n * p * (1 - p))))/(2 * (n + z^2))
  upper <- (2 * n * p + z^2 + (z * sqrt(z^2 + 4 * n * p * (1 - p))))/(2 * (n + z^2))
  if (!return_string) return(list(p = p, l = lower, u = upper))
  paste0("[", prmisc::decimals_only(lower, decimals), ", ", prmisc::decimals_only(upper, decimals), "]")
}
get.z.score <- function(ci) {
  alpha <- (100 - ci)/100
  z <- qnorm((1 - alpha/2))
  return(z)
}

color_rnd <- "#D55E00"
color_anticlust <- "#0072B2"
color_confound <-  "#009E73"


ldf$K <- factor(ldf$K)
# Power by covariate effect
pd <- position_dodge(.05)
ldf |>
  filter(treatment_effect > 0, !grepl("confound", Method)) |>
  group_by(K, Method, scale_batch_effect) |>
  summarize(
    Power = mean(pvalue <= .05),
    lower = ci_one_prop(r = sum(pvalue <= .05), n = n(), return_string = FALSE)$l,
    upper = ci_one_prop(r = sum(pvalue <= .05), n = n(), return_string = FALSE)$u) |>
  ggplot(aes(x = K, y = Power, colour = Method)) +
  geom_point(aes(color = Method, shape = Method), position = pd, size = 3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, colour = Method, linetype = Method), 
    position = pd, width = .3
  ) +
  facet_grid(cols = vars(scale_batch_effect), scales = "free") +
  xlab("Number of batches") +
  theme_bw(base_size = 16)+
  theme(legend.position = "top", legend.box = "vertical", legend.title = element_blank())
  # scale_color_manual(values = c(color_anticlust, color_minimization, color_rnd))



