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
# may lead to exaggerated confidence in downstream analyses")
# as implemented in the afex package 

library(anticlust) # for anticlustering
library(afex) # for 2-way ANOVA
library(prmisc) # for formatting a p value. 

source("functions.R") # functions that implement the data generating process + statistical control

# set.seed(123)
nsim <- 1000
start <- Sys.time()
results <- replicate(
  nsim, 
  simulate(
    N = 200,  # number of samples
    K = 20,   # number of batches
    M = 10,    # number of covariates
    P = 10,   # number of class per covariate
    scale_batch_effect = 10, # relative effect of batches
    treatment_effect = 1.5,   # effect size of a treatment
    SD_residual = 2
  )
)
Sys.time() - start

# replicate(1000, foo(N = 200, M = 5, P = 10, scale_batch_effect = 20, confound = FALSE))
# -> finds significant, but marginal improvement of anticlustering

sort(rowMeans(results < .05)) |> round(2)
# Uncontroled confounded analysis has very high "power". It also produces alpha errors because it confounds batch
# effects (which are present) for treatment effects. This is very bad.
cat("p value of anticlustering was lower (=better) in ", 
    mean(results["p_anticlust_control", ] < results["p_rnd_control", ]) * 100,
    "% of all cases,", 
    format_p(prop.test(sum(results["p_anticlust_control", ] < results["p_rnd_control", ]), ncol(results), p = .50)$p.value)
)

t.test(log(results["p_anticlust_control", ]), log(results["p_rnd_control", ]), paired = TRUE) # anticlust significantly improves inference over mere statistical control + random assignment
