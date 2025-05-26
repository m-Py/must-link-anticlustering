 
N = 200; K = 2; M = 1; P = 2;
scale_batch_effect = 10; SD_residual = 2;
treatment_effect = 1; prob_treatment = .5; adjustment_method = "afex";

simulate <- function(X, 
    N = 200, K = 20, M = 10, P = 5, 
    scale_batch_effect = 10, SD_residual = 2, 
    treatment_effect = 1.2, prob_treatment = .5, 
    adjustment_method = c("lm", "afex")) {
  
  adjustment_method <- ifelse(length(adjustment_method) == 1, adjustment_method, adjustment_method[1])
  covariates <- generate_categorical_data(N, M, P = P)
  covariates_binary <- categories_to_binary(covariates)
  b1 <- rnorm(ncol(covariates_binary)) # effect of covariate on outcome
  b2 <- treatment_effect #effect of treatment on outcome
  
  batches_rnd <- sample(rep_len(1:K, N))
  treatment <- sample(0:1, size = N, replace = TRUE, prob = c(1-prob_treatment, prob_treatment))

  batches_anticlust <- fast_anticlustering(cbind(covariates_binary, treatment), K = K) # balance treatment + covariates
  b0 <- sample(((1:K) / K) * scale_batch_effect)
  residual <- rnorm(N, sd = SD_residual)
  
  outcome_rnd <- get_batch_data(
    N, batches_rnd, 
    treatment, 
    covariates_binary, 
    b0, b1, b2
  ) + residual
  outcome_anticlust <- get_batch_data(
    N, batches_anticlust, 
    treatment, 
    covariates_binary, 
    b0, b1, b2
  ) + residual

  c(
    p_rnd_no_control = get_p_value_treatment(N, outcome_rnd, treatment, batches_rnd, FALSE, adjustment_method),
    p_rnd_control = get_p_value_treatment(N, outcome_rnd, treatment, batches_rnd, TRUE, adjustment_method),
    p_anticlust_no_control = get_p_value_treatment(N, outcome_anticlust, treatment, batches_anticlust, FALSE, adjustment_method),
    p_anticlust_control = get_p_value_treatment(N, outcome_anticlust, treatment, batches_anticlust, TRUE, adjustment_method),
    cor_rnd = cor_batch_effect_treatment_effect(batches_rnd, treatment, b0, b2),
    cor_anticlust = cor_batch_effect_treatment_effect(batches_anticlust, treatment, b0, b2)
  )
}
