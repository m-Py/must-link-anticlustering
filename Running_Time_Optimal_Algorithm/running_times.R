# Running times 

library(anticlust)

nrep <- 5
results_file <- "results_runningtime.csv"

solver <- c("gurobi") # you can use lpSolve, glpk or symphony as open source options; we will not get as far as with gurobi though

## Set min and max for sample sizes by K
Ns <- list()
Ns[["K=2"]] <- c(min = 20, max = 1000)
Ns[["K=3"]] <- c(min = 21, max = 999)
Ns[["K=4"]] <- c(min = 20, max = 1000)
Ns[["K=5"]] <- c(min = 20, max = 1000)
Ns[["K=10"]] <- c(min = 20, max = 1000)

TIME_LIMIT <- 1800

# this became somewhat complicated so it is a function now:
get_next_N <- function(N, K, N_max) {
  if (N < 100) {
    step <- K
  } else if (N < 500) {
    step <- 10 * K
  } else {
    step <- 20 * K 
  } #  lateron, increase faster
  if (N < N_max) {
    N <- N + step
  }
  if (N > N_max) {
    N <- N_max # always use last iteration with N_max
  } else if (N == N_max) {
    N <- N + 1 # just to break out of while loop after last iteration
  }
  N
}

generate_categorical_data <- function(N, M, P) {
  data <- data.frame(matrix(sample(P, replace = TRUE, size = N*M), ncol = M))
  as.data.frame(lapply(data, as.factor))
}

edit_must_link_distances <- function(distances, must_link) {
  INFINITY <- sum(distances) + 1
  distances <- as.matrix(distances)
  same <- as.numeric(names(table(must_link)[table(must_link) > 1]))
  for (i in same) {
    distances[must_link == i, must_link == i] <- INFINITY
  }
  distances
}

must_link_constraints_valid <- function(cl, must_link) {
  same <- as.numeric(names(table(must_link)[table(must_link) > 1]))
  all_good <- c()
  for (i in same) {
    all_good <- c(all_good, all(cl[must_link == i] == cl[must_link == i][1]))
  }
  all(all_good)
}

## I need functions for testing equality /smaller than with some allowance for numeric imprecision...
"%==%" <- function(x, y) abs(x - y) < .0000000001
"%>%" <- function(x, y) (x - y) > (.0000000001)



for (K in c(2:5, 10)) {
  N <- Ns[[paste0("K=", K)]]["min"]
  time_limit_exceeded <- FALSE # for ending the while loop
  N_max <- Ns[[paste0("K=", K)]]["max"]
  while (!time_limit_exceeded && N <= N_max) {
    message("N = ", N, ", K = ", K, ". Time limit: ", TIME_LIMIT, "s.")
    for (j in 1:nrep) {
      P <- sample(2:5, size = 1)
      M <- sample(2:5, size = 1)
      data <- generate_categorical_data(N, M, P)
      distances <- as.matrix(dist(categories_to_binary(generate_categorical_data(N, M, P))))
      
      must_link <- sample(N, replace = TRUE)
      
      MUST_LINK_CONSTRAINS_FULFILLED <- FALSE
      while (!MUST_LINK_CONSTRAINS_FULFILLED) {
        start <- Sys.time()
        heuristic <- tryCatch(
          anticlustering(distances, K = K, must_link = must_link, method = "2PML", repetitions = 1000),
          error = function(e) e
        )
        time1 <- as.numeric(difftime(Sys.time(), start, units = "s")) |> round(2)
        if ("simpleError" %in% class(heuristic)) { 
          # re must-link constraints if must link constraints cannot be fulfilled
          must_link <- sample(N, replace = TRUE)
        } else {
          MUST_LINK_CONSTRAINS_FULFILLED <- TRUE
        }
      }
      
      ## OPTIMAL ANTICLUSTERING GIVEN CONSTRAINTS
      start <- Sys.time()
      edited_distances <- edit_must_link_distances(distances, must_link)
      opt <- tryCatch(
        optimal_anticlustering(edited_distances, K, objective = "diversity", solver = solver, time_limit = TIME_LIMIT),
        error = function(e) e
      )
      time2 <- as.numeric(difftime(Sys.time(), start, units = "s")) |> round(2)
      if ("simpleError" %in% class(opt)) { 
        time_limit_exceeded <- TRUE
        message("Time limit exceeded!")
      }

      #################
      
      if (!time_limit_exceeded) {
        stopifnot(must_link_constraints_valid(heuristic, must_link))
        stopifnot(must_link_constraints_valid(opt, must_link))
        
        diversity_heuristic <- diversity_objective(distances, heuristic)
        diversity_optimal   <- diversity_objective(distances, opt)
        if (diversity_optimal %>% diversity_heuristic) {
          cat("Heuristic did not find optimal solution.\n")
        } else if (diversity_heuristic %==% diversity_optimal) {
          cat("Heuristic found optimal solution.\n")
        } else {
          write.table(
            distances, 
            paste0("./logs/DISTANCS_N=", N, "_K=", K, "_ERROR_ILP.csv"),
            quote = FALSE,
            row.names = FALSE,
            sep = ";"
          )
          write.table(
            edited_distances, 
            paste0("./logs/EDITED_DISTANCS_N=", N, "_K=", K, "_ERROR_ILP.csv"),
            quote = FALSE,
            row.names = FALSE,
            sep = ";"
          )
          cat("Something went wrong, ILP did not find optimal solution. Data set was stored.") #FYI this does not happen
        }
      } else {
        diversity_heuristic <- NA
        diversity_optimal <- NA
      }

      
      ## Write results to file 
      df <- data.frame(
        N = N,
        K = K,
        diversity_heuristic = diversity_heuristic,
        diversity_optimal = diversity_optimal,
        time_limit = TIME_LIMIT,
        time_limit_exceeded = time_limit_exceeded,
        time_heuristic = time1, 
        time_optimal = time2,
        solver = solver
      )
      results_file_exists <- FALSE
      if (file.exists(results_file)) {
        results_file_exists <- TRUE
      }
      write.table(
        df, 
        file = results_file, 
        append = results_file_exists, 
        col.names = !results_file_exists,
        row.names = FALSE, 
        quote = FALSE,
        sep = ";"
      )
      if (time_limit_exceeded) {
        break
      }
    }
    N <- get_next_N(N, K, N_max)
  }
}
