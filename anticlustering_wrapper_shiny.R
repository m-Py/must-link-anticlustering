# anticlustering wrapper function for shiny

library(anticlust)
## Error handling must be done (at least partly) outside of this function
# This function *always* implements standardization of the input variables
anticlustering_shiny <- function(
    numeric_vars = NULL, categorical_vars = NULL, 
    must_link_constraints = NULL, K
) {

  if (!is.null(must_link_constraints)) {
    return(must_link_anticlustering_shiny(numeric_vars, categorical_vars, must_link_constraints, K))
  }
  
  
  N <- max(nrow(numeric_vars), nrow(categorical_vars)) # deal with NULL input
  if (N > 1000) { # "Large" data sets: 
    if (!is.null(must_link_constraints)) stop("deal with must link constraints and large N")
    return(fast_anticlustering_shiny(numeric_vars, categorical_vars, K))
  }
  
  if (N <= 1000) {
    repetitions <- 1
    reps <- "1 repetition"
  }
  if (N <= 200) {
    repetitions <- 10
    reps <- "10 repetitions"
  }
  if (N <= 100) {
    repetitions <- 100
    reps <- "100 repetitions"
  }
  if (!is.null(categorical_vars)) {
    categorical_vars <- categories_to_binary(categorical_vars)
  }
  input <- cbind(numeric_vars, categorical_vars)
  message("N = ", N, ". Using anticlustering with diversity criterion, method is 'local-maximum' with ", reps, ".\n")
  anticlustering(input, K = K, method = "local-maximum", repetitions = repetitions, standardize = TRUE)
}


must_link_anticlustering_shiny <- function(numeric_vars, categorical_vars, must_link_constraints, K) {
  if (!is.null(categorical_vars)) {
    categorical_vars <- categories_to_binary(categorical_vars)
  }
  input <- cbind(numeric_vars, categorical_vars)
  N <- nrow(input) 
  repetitions <- 2
  if (N > 15000) {
    stop("Sorry, we do not offer the must-link for data sets with more than 15000 elements. Try out the R package anticlust directly.")
  }
  if (N <= 1000) {
    repetitions <- 1
  }
  if (N <= 200) {
    repetitions <- 10
  }
  if (N <= 100) {
    repetitions <- 100
  }
  reps <- paste0(repetitions, " repetitions (", repetitions/2, "x Phase 1, ", repetitions/2, "x Phase 2)")
  
  message("N = ", N, ". Using anticlustering with diversity criterion, method is '2PML' (for must-link constraints) with ", reps, ".\n")
  anticlustering(input, K = K, method = "2PML", repetitions = repetitions, standardize = TRUE, must_link = must_link_constraints)
  
}

fast_anticlustering_shiny <- function(numeric_vars, categorical_vars, K) {
  if (!is.null(categorical_vars)) {
    categorical_vars <- categories_to_binary(categorical_vars)
  }
  if (!is.null(numeric_vars)) {
    numeric_vars <- kplus_moment_variables(numeric_vars, 2, FALSE)
    criterion <- "k-plus"
  } else {
    criterion <- "k-means"
  }
  input <- cbind(numeric_vars, categorical_vars)
  input <- scale(input) # use standardization
  if (N > 10000) { # "very large data"
    nn_method <- ifelse(N > 100000, "random", "RANN")
    message(
      "N = ", N, ". Using ", criterion, 
      " anticlustering with fast exchange method (100 exchange partners, selection of exchange partners using method ", 
      nn_method, ").\n")
    return(
      fast_anticlustering(
        input, K = K,
        exchange_partners = generate_exchange_partners(100, N = N, features = numeric_vars, method = nn_method)
      )
    )
  } else {
    message("N = ", N, ". Using ", criterion, " anticlustering with exchange method.")
    return(fast_anticlustering(input, K = K))
  }
}
