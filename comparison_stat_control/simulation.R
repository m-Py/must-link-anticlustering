# Small scale simulation that compares anticlustering-based assignment versus random assignment and 
# confounded assignment. For each assignment technique, we use statistical control of batches, or not.

# Outcome (O) is modelled as:
# O = b0*C + b1*X + b2*T + e 
#  b0 = confounding effect of batch, C = Confounding Batch
#  b1 = Effect of Covariates, X = covariates
#  b2 = Effect of treatment, T = treatment (1 / 0)
#  e = residual error

# Simulation uses the following design: Some factors are crossed, other variables randomly vary

# Factors: 
# - Effect size treatment: 0, 0.5, 1, 1.5, 2
# - Effect size interaction covariate / treatment: 0, 1, 2 ("none", "small", "medium", "large")
# - "Effect size" (scaling factor) batches: "small" (1), "substantial" (10)
# - Number of batches: 2, 5, 10, 20

# Fixed: 
# - Residual error SD = 2
# - Covariate Effect = 1
# - Number of covariates = 1
# - Number of levels in covariate = 2
# - Probability that treatment is received = 0.5

# Randomly varies: 
# - sample_sizes <- 40:400; N <- sample(sample_sizes[sample_sizes %% K == 0], 1) # sample size

# For statistical control, we use 2-way ANOVA (see NYGAARD et al., 2016; 
# "Methods that remove batch effects while retaining group differences 
# may lead to exaggerated confidence in downstream analyses"),
# Here, it is equivalent to simply using lm() (because we only use additive effects).

library(anticlust) # for anticlustering
library(prmisc) # for formatting a p value. 
library(here)
library(parallel)
# functions that implement the data generating process + statistical control
source(here("comparison_stat_control", "simulation_function.R")) 

# only create cluster object once at the start to not confuse the computer
if (!"cl" %in% ls()) {
  cl <- makeCluster(getOption("cl.cores", min(12, detectCores() - 2))) # use 12 if possible but always leave 2 cores alone;
}

conditions <- expand.grid(
  scale_batch_effect = c(1, 10),
  covariate_effect = c(1),
  treatment_effect = c(0, 0.5, 1, 1.5),
  interaction_effect = c(0, 2.5),
  K = c(2, 5, 10, 20), 
  prob_treatment = c(.50)
)
nrow(conditions)

nsim <- 4000 # number of repetitions per simulation conditions

start <- Sys.time()
for (i in 1:nrow(conditions)) {
  results <- parSapply(
    cl, 
    1:nsim, 
    simulate_parallel,
    K = conditions[i, "K"],   # number of batches.
    covariate_effect = conditions[i, "covariate_effect"], # relative effect of covariate on outcome
    scale_batch_effect = conditions[i, "scale_batch_effect"], # relative effect of batches. Must be somewhat large to see effect of stat. control. 
    treatment_effect = conditions[i, "treatment_effect"],     # effect size of a treatment (0 = null effect; check for alpha errors)
    interaction_effect = conditions[i, "interaction_effect"],
    #(note that these effect sizes are not really comparable in their magnitude; 
    # they are differently related to the regression that creates the data)
    # scale_batch_effect is useful to obtain a discrepancy in power between adjusted and unadjusted analysis
    SD_residual = 2, # residual error SD
    prob_treatment = conditions[i, "prob_treatment"] # case / control has the same probability with prob_treatment = .5
  )
  results <- as.data.frame(t(results))
  
  write.table(
    results,
    file = here("comparison_stat_control", "data", 
                paste0(format(Sys.time(), "%Y-%m-%d"), "_DAT_", paste0(sample(LETTERS, 5), sample(0:9, 5), collapse = ""), ".csv")),
    row.names = FALSE,
    quote = FALSE,
    sep = ";"
  )
}
Sys.time() - start 
