
# N = 200; K = 20; M = 10; P = 5; 
# scale_batch_effect = 10; SD_residual = 2;
# treatment_effect = 1.2

simulate <- function(
    N = 200, K = 20, M = 10, P = 5, 
    scale_batch_effect = 10, SD_residual = 2, 
    treatment_effect = 1.2) {
  covariates <- generate_categorical_data(N, M, P = P)
  covariates_binary <- categories_to_binary(covariates)
  b1 <- rnorm(ncol(covariates_binary)) # effect of covariate on outcome
  b2 <- treatment_effect #effect of treatment on outcome
  
  batches_rnd <- sample(rep_len(1:K, N))
  treatment_rnd <- sample(0:1, size = N, replace = TRUE) 
  treatment_confound <- batches_rnd + rnorm(N, sd = 10) < (K/2) 
  
  batches_anticlust <- fast_anticlustering(cbind(covariates_binary, treatment_rnd), K = K) # balance treatment + covariates
  b0 <- ((1:K) / K) * scale_batch_effect # batch effect on outcome
  
  residual <- rnorm(N, sd = SD_residual)
  
  outcome_rnd <- get_batch_data(
    N, batches_rnd, 
    treatment_rnd, 
    covariates_binary, 
    b0, b1, b2
  ) + residual
  outcome_anticlust <- get_batch_data(
    N, batches_anticlust, 
    treatment_rnd, 
    covariates_binary, 
    b0, b1, b2
  ) + residual
  outcome_confound <- get_batch_data(
    N, batches_rnd, 
    treatment_confound, 
    covariates_binary, 
    b0, b1, b2
  ) + residual

  
  c(
    p_rnd_no_control = get_p_value_treatment(N, outcome_rnd, treatment_rnd, batches_rnd, FALSE),
    p_rnd_control = get_p_value_treatment(N, outcome_rnd, treatment_rnd, batches_rnd, TRUE),
    p_confound_no_control = get_p_value_treatment(N, outcome_confound, treatment_confound, batches_rnd, FALSE),
    p_confound_control = get_p_value_treatment(N, outcome_confound, treatment_confound, batches_rnd, TRUE),
    p_anticlust_no_control = get_p_value_treatment(N, outcome_anticlust, treatment_rnd, batches_anticlust, FALSE),
    p_anticlust_control = get_p_value_treatment(N, outcome_anticlust, treatment_rnd, batches_anticlust, TRUE)
  )
}

get_p_value_treatment <- function(N, outcome, treatment, batches, statistical_adjustment)  {
  
  dat <- data.frame(
    id = 1:N, 
    outcome = outcome,
    treatment = treatment,
    batch = batches
  )
  if (statistical_adjustment) {
    analysis <- suppressMessages(aov_ez(
      id = "id",
      dv = "outcome", 
      between = "treatment",
      covariate = "batch",
      data = dat
    ))
  } else {
    analysis <- suppressMessages(aov_ez(
      id = "id",
      dv = "outcome", 
      between = "treatment",
      data = dat
    ))
  }
  analysis$anova_table["treatment", "Pr(>F)"] 
}

# get outcome data depending on assignment scheme; without residual error though (is added later, the same for all assignments to reduce error variance between assignment schemes in the simulation)
get_batch_data <- function(N, batches, treatment, covariates, b0, b1, b2) {
  stopifnot(length(batches) == N)
  stopifnot(length(b0) == length(unique(batches)))
  covariate_effect <- c(covariates %*% b1) # get as vector
  batch_effect <- c(categories_to_binary(batches) %*% b0)
  treatment_effect <- treatment * b2 # treatment = 0/1 coded
  batch_effect + covariate_effect + treatment_effect
}

generate_categorical_data <- function(N, M, P, prob = NULL) { # prob given to sample
  if (is.null(M)) return(NULL)
  stopifnot(is.numeric(N) && is.numeric(M) && is.numeric(P))
  stopifnot(length(N) == 1 && length(M) == 1 && length(P) == 1)
  if (is.null(P)) stop("Number of classes must be passed via argument P for categorical variables")
  data <- data.frame(matrix(sample(P, replace = TRUE, size = N * M, prob = prob), ncol = M))
  as.data.frame(lapply(data, as.factor))
}
