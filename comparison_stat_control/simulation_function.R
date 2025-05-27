
# K = 20
# scale_batch_effect = 10
# covariate_effect = 1
# SD_residual = 2
# treatment_effect = 1.2
# prob_treatment = .5
# adjust_for_covariate = FALSE

simulate_parallel <- function(X, 
    K = 20, 
    scale_batch_effect = 10, 
    covariate_effect = 1,
    SD_residual = 2, 
    treatment_effect = 1.2, 
    prob_treatment = .5, 
    adjust_for_covariate = FALSE) {

  # Set variables that randomly vary
  sample_sizes <- 40:400
  N <- sample(sample_sizes[sample_sizes %% K == 0], 1) 

  # LOAD LIBRARIES AND FUNCTIONS
  library(anticlust)
  ## DEFINE HELPER FUNCTIONS HERE BECAUSE OTHERWISE IN PARALLEL THE CORES DO NOT HAVE THEM
  
  cor_batch_effect_treatment_effect <- function(batches, treatment, b0, b2) {
    batch_effect <- c(categories_to_binary(batches) %*% b0)
    treatment_effect <- treatment * b2
    cor(batch_effect, treatment_effect)
  }
  
  get_p_value_treatment <- function(N, outcome, treatment, batches, covariate, statistical_adjustment, adjust_for_covariate)  {

    stopifnot(length(covariate) == N)
    stopifnot(length(outcome) == N)
    stopifnot(length(treatment) == N)
    stopifnot(nrow(batches) == N)
    stopifnot(length(statistical_adjustment) == 1)
    stopifnot(length(adjust_for_covariate) == 1)
    
    dat <- data.frame(
      id = 1:N, 
      outcome = outcome,
      treatment = treatment,
      batch = as.factor(batches),
      covariate = covariate
    )
    
    model <- "outcome ~ treatment"
    
    if (statistical_adjustment) {
      model <- paste(model, "+ batch")
    } 
    
    if (adjust_for_covariate) {
      model <- paste(model, "+ covariate")
    }
    coef(summary(lm(as.formula(model), data = dat)))["treatment", "Pr(>|t|)"]
  }
  
  # get outcome data depending on assignment scheme; without residual error though
  # (residual is added later on output of this function for different assignments,
  # so it is the same for all assignments to reduce error variance between
  # assignment schemes in the simulation)
  get_batch_data <- function(N, batches, treatment, covariate, b0, b1, b2) {
    stopifnot(length(batches) == N)
    stopifnot(length(b0) == length(unique(batches)))
    covariate_effect <- c(as.matrix(covariate) %*% b1)
    batch_effect <- c(categories_to_binary(batches) %*% b0) # need as vector
    treatment_effect <- treatment * b2 # treatment = 0/1 coded
    batch_effect + covariate_effect + treatment_effect
  }
  
  covariate <- sample(0:1, size = N, replace = TRUE)
  b1 <- covariate_effect # effect of covariate on outcome
  b2 <- treatment_effect #effect of treatment on outcome
  
  batches_rnd <- sample(rep_len(1:K, N))
  treatment <- sample(0:1, size = N, replace = TRUE, prob = c(1-prob_treatment, prob_treatment))

  batches_anticlust <- fast_anticlustering(cbind(covariate, treatment), K = K) # balance treatment + covariates
  b0 <- sample(((1:K) / K) * scale_batch_effect)
  residual <- rnorm(N, sd = SD_residual)
  
  batches_confounded <- balanced_clustering(treatment, K = K)
  
  outcome_rnd <- get_batch_data(
    N, batches_rnd, 
    treatment, 
    covariate, 
    b0, b1, b2
  ) + residual
  outcome_anticlust <- get_batch_data(
    N, batches_anticlust, 
    treatment, 
    covariate, 
    b0, b1, b2
  ) + residual
  outcome_confound <- get_batch_data(
    N, batches_confounded, 
    treatment, 
    covariate, 
    b0, b1, b2
  ) + residual

  c(
    ID = as.numeric(paste0(sample(0:9, size = 15, replace = TRUE), collapse = "")),
    N = N, K = K, 
    scale_batch_effect = scale_batch_effect, 
    covariate_effect = covariate_effect,
    SD_residual = SD_residual, 
    treatment_effect = treatment_effect, 
    prob_treatment = prob_treatment, 
    adjust_for_covariate = adjust_for_covariate,
    p_rnd_no_control = get_p_value_treatment(N, outcome_rnd, treatment, batches_rnd, covariate, FALSE, FALSE),
    p_rnd_control = get_p_value_treatment(N, outcome_rnd, treatment, batches_rnd, covariate, TRUE, adjust_for_covariate),
    p_anticlust_no_control = get_p_value_treatment(N, outcome_anticlust, treatment, batches_anticlust, covariate, FALSE, FALSE),
    p_anticlust_control = get_p_value_treatment(N, outcome_anticlust, treatment, batches_anticlust, covariate, TRUE, adjust_for_covariate),
    p_confound_no_control = get_p_value_treatment(N, outcome_confound, treatment, batches_confounded, covariate, FALSE, FALSE),
    p_confound_control = get_p_value_treatment(N, outcome_confound, treatment, batches_confounded, covariate, TRUE, adjust_for_covariate),
    cor_rnd = ifelse(treatment_effect > 0, cor_batch_effect_treatment_effect(batches_rnd, treatment, b0, b2), NA),
    cor_anticlust = ifelse(treatment_effect > 0, cor_batch_effect_treatment_effect(batches_anticlust, treatment, b0, b2), NA),
    cor_confound = ifelse(treatment_effect > 0, cor_batch_effect_treatment_effect(batches_confounded, treatment, b0, b2), NA),
    balance_rnd =  chisq.test(batches_rnd, treatment)$p.value,
    balance_anticlust = chisq.test(batches_anticlust, treatment)$p.value,
    balance_confound = chisq.test(batches_confounded, treatment)$p.value
  )
}
