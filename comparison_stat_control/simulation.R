# Small scale simulation that compares anticlustering-based assignment versus random assignment and 
# confounded assignment. For each assignment technique, we use statistical control of batches. 

# Outcome (O) is modelled as:
# O = b0*C + b1*X + b2*T + e 
#  b0 = confounding effect of batch, C = Confounding Batch
#  b1 = Effect of Covariates, X = covariates
#  b2 = Effect of treatment, T = treatment (1 / 0)
#  e = residual error

# For statistical control, we use 2-way ANOVA (see NYGAARD et al., 2016; 
# "Methods that remove batch effects while retaining group differences 
# may lead to exaggerated confidence in downstream analyses"),
# as implemented in the afex package 

library(anticlust) # for anticlustering
library(afex) # for 2-way ANOVA
library(prmisc) # for formatting a p value. 
library(here)

source(here("comparison_stat_control", "functions.R")) # functions that implement the data generating process + statistical control

# set.seed(123)
nsim <- 1000
start <- Sys.time()
results <- replicate(
  nsim, 
  simulate(
    N = 200,  # number of samples
    K = 20,   # number of batches. Seems that anticlustering improves inferences for smaller batches (i.e., more batches with constant N)
    M = 4,    # number of covariates
    P = 4,    # number of class per covariate
    scale_batch_effect = 10,  # relative effect of batches. Must be somewhat large (e.g., >= 10). 
    treatment_effect = 1,     # effect size of a treatment (0 = null effect; check for alpha errors)
    #(note that these effect sizes are not really comparable in their magnitude; 
    # they are differently related to the regression that creates the data)
    # scale_batch_effect is useful to obtain a discrepancy in power between adjusted and unadjusted analysis
    SD_residual = 2, # residual error SD
    prob_treatment = .5 # case / control has the same probability with prob_treatment = .5
  )
)
Sys.time() - start

# overall power per method:
sort(rowMeans(results < .05)) |> round(2)

cat("For adjusted analysis: p value of anticlustering was lower (=better) in ", 
    mean(results["p_anticlust_control", ] < results["p_rnd_control", ]) * 100,
    "% of all cases,", 
    format_p(prop.test(sum(results["p_anticlust_control", ] < results["p_rnd_control", ]), ncol(results), p = .50)$p.value),
    "\n"
)

cat("For unadjusted analysis: p value of anticlustering was lower (=better) in ", 
    mean(results["p_anticlust_no_control", ] < results["p_rnd_no_control", ]) * 100,
    "% of all cases,", 
    format_p(prop.test(sum(results["p_anticlust_no_control", ] < results["p_rnd_no_control", ]), ncol(results), p = .50)$p.value),
    "\n"
)
