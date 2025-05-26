
cor_batch_effect_treatment_effect <- function(batches, treatment, b0, b2) {
  batch_effect <- c(categories_to_binary(batches) %*% b0)
  treatment_effect <- treatment * b2
  cor(batch_effect, treatment_effect)
}

get_p_value_treatment <- function(N, outcome, treatment, batches, statistical_adjustment, adjustment_method)  {
  
  dat <- data.frame(
    id = 1:N, 
    outcome = outcome,
    treatment = treatment,
    batch = as.factor(batches)
  )
  
  if (statistical_adjustment && adjustment_method == "afex") {
    pvalue <- suppressMessages(aov_ez(
      id = "id",
      dv = "outcome", 
      between = "treatment",
      covariate = "batch",
      data = dat
    ))$anova_table["treatment", "Pr(>F)"] 
  } else if (!statistical_adjustment && adjustment_method == "afex") {
    pvalue <- suppressMessages(aov_ez(
      id = "id",
      dv = "outcome", 
      between = "treatment",
      data = dat
    ))$anova_table["treatment", "Pr(>F)"] 
  } else if (statistical_adjustment && adjustment_method == "lm") {
    pvalue <- coef(summary(lm(outcome ~ treatment + batch, data = dat)))["treatment", "Pr(>|t|)"]
  } else if (!statistical_adjustment && adjustment_method == "lm") {
    pvalue <- coef(summary(lm(outcome ~ treatment, data = dat)))["treatment", "Pr(>|t|)"]
  }
  pvalue
}




# get outcome data depending on assignment scheme; without residual error though
# (residual is added later on output of this function for different assignments,
# so it is the same for all assignments to reduce error variance between
# assignment schemes in the simulation)
get_batch_data <- function(N, batches, treatment, covariates, b0, b1, b2, b3) {
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
