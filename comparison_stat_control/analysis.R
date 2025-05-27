
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(santoku)

# read all data files
files <- list.files(here("comparison_stat_control", "data"), full.names = TRUE)

df <- do.call("rbind", lapply(files, read.csv, sep =";"))

head(df)

df1 <- df[df$treatment_effect > 0, ]
df0 <- df[df$treatment_effect == 0, ]

mean(df1$p_anticlust_control < df1$p_rnd_control)

plot_stuff <- function(cor_rnd, p_rnd, cor_anticlust, p_anticlust, cor_confound=NULL, p_confounded=NULL,
                       xlab = "correlation batch effect / treatment effect", ylab = "p value") {
  xlim <- range(cor_confound, cor_anticlust, cor_confound)
  plot(cor_rnd, p_rnd, xlim = xlim, pch = 2, cex = .5, col = "#ABCDEF",
       xlab = xlab, ylab = ylab)
  points(cor_anticlust, p_anticlust, col = "#333333", cex = .6, pch = 4)
  if(!is.null(cor_confound)) {
    points(cor_confound, p_confounded, col = "#f0f", cex = .6, pch = 4)
  }
  legend <-  c("Random Assignment", "Anticlustering")
  if (!is.null(cor_confound)) {
    legend <- c(legend, "Confounded Assignment")
  }
  legend("topright", legend = legend, 
         pch = c(2, 4), col = c("#ABCDEF", "#333333", "#f0f"), cex = .7)
}


# plot_stuff(
#   df1$cor_rnd, df1$p_rnd_no_control, 
#   df1$cor_anticlust, df1$p_anticlust_no_control
# )

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

ldf <- df |>
  select(ID, N, K, scale_batch_effect, treatment_effect, starts_with("p_")) |>
  pivot_longer(
    cols = starts_with("p_"),
    names_prefix = "p_",
    names_to = "Method",
    values_to = "pvalue"
  )

ldf |>
  filter(treatment_effect > 0) |>
  group_by(Method) |>
  summarize(
    Power = mean(pvalue <= .05), 
    `95% CI` = ci_one_prop(r = sum(pvalue <= .05), n = n()), n = papaja::printnum(n(), format = "d")
  ) |>
  arrange(Power) |> 
  mutate(Power = prmisc::decimals_only(Power, 3))


ldf$`Statistical Control` <- "unadjusted analysis"
ldf$`Statistical Control`[!grepl("no_control", ldf$Method)] <- "adjusted analysis"

ldf$Assignment <- "Random Assignment"
ldf$Assignment[grepl("anticlust_", ldf$Method)] <- "Anticlustering"
ldf$Assignment[grepl("confound_", ldf$Method)] <- "Confounded"

ldf$K <- factor(ldf$K)
ldf$Batch_Effects <- "small batch effects"
ldf$Batch_Effects[ldf$scale_batch_effect == 10] <- "large batch effects"

# Power by batch effect and number of batches
pd <- position_dodge(.05)
ldf |>
  filter(treatment_effect > 0, !grepl("confound", Method)) |>
  group_by(K, Method, Assignment, `Statistical Control`, Batch_Effects) |>
  summarize(
    Power = mean(pvalue <= .05),
    lower = ci_one_prop(r = sum(pvalue <= .05), n = n(), return_string = FALSE)$l,
    upper = ci_one_prop(r = sum(pvalue <= .05), n = n(), return_string = FALSE)$u) |>
  ggplot(aes(x = K, y = Power, colour = Method)) +
  geom_point(aes(color = Assignment, shape = `Statistical Control`), position = pd, size = 2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, colour = Assignment, linetype = `Statistical Control`), 
    position = pd, width = .3
  ) +
  facet_grid(cols = vars(Batch_Effects),scales = "free") +
  xlab("Number of batches") +
  theme_bw(base_size = 16)+
  theme(legend.position = "top", legend.box = "vertical", legend.title = element_blank())
  # scale_color_manual(values = c(color_anticlust, color_minimization, color_rnd))

# Interesting: With considerable batch effects, anticlustering is even advantageous for unadjusted analysis


# Power by N
ldf$N_category <- santoku::chop(ldf$N, breaks = c(100, 300))
ldf |>
  filter(treatment_effect > 0, !grepl("confound", Method)) |>
  group_by(N_category, Method, Assignment, `Statistical Control`, Batch_Effects) |>
  summarize(
    Power = mean(pvalue <= .05),
    lower = ci_one_prop(r = sum(pvalue <= .05), n = n(), return_string = FALSE)$l,
    upper = ci_one_prop(r = sum(pvalue <= .05), n = n(), return_string = FALSE)$u) |>
  ggplot(aes(x = N_category, y = Power, colour = Method)) +
  geom_point(aes(color = Assignment, shape = `Statistical Control`), position = pd, size = 2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, colour = Assignment, linetype = `Statistical Control`), 
    position = pd, width = .3
  ) +
  facet_grid(cols = vars(Batch_Effects), scales = "free") +
  xlab("Number of samples (categorized)") +
  theme_bw(base_size = 16)+
  theme(legend.position = "top", legend.box = "vertical", legend.title = element_blank())
# scale_color_manual(values = c(color_anticlust, color_minimization, color_rnd))


# Power by N and K and batch effects
ldf$N_category <- ordered(ldf$N_category, levels = levels(ldf$N_category), labels = paste0("N = ", levels(ldf$N_category)))
ldf |>
  filter(treatment_effect > 0, !grepl("confound", Method)) |>
  group_by(N_category, K, Method, Assignment, `Statistical Control`, Batch_Effects) |>
  summarize(
    Power = mean(pvalue <= .05),
    lower = ci_one_prop(r = sum(pvalue <= .05), n = n(), return_string = FALSE)$l,
    upper = ci_one_prop(r = sum(pvalue <= .05), n = n(), return_string = FALSE)$u) |>
  ggplot(aes(x = K, y = Power, colour = Method)) +
  geom_point(aes(color = Assignment, shape = `Statistical Control`), position = pd, size = 2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, colour = Assignment, linetype = `Statistical Control`), 
    position = pd, width = .3
  ) +
  facet_grid(cols = vars(Batch_Effects), vars(N_category), scales = "free") +
  xlab("Number of covariates") +
  theme_bw(base_size = 16)+
  theme(legend.position = "top", legend.box = "vertical", legend.title = element_blank())
# scale_color_manual(values = c(color_anticlust, color_minimization, color_rnd))



