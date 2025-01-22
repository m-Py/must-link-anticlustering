
# Author: Martin Papenberg
# Year: 2025


## For Phase 2 of the 2PML algorithm, we generate all combinations of numbers 
# (representing the sizes of cliques / singletons) sizes that in total add
# up to the size of a clique for which we search exchange partners. 
# Generating these combination in fact constitutes the NP hard subset sum problem.
# https://en.wikipedia.org/wiki/Subset_sum_problem
# So we can't do this exhaustively for arbitrary sizes of cliques. The code
# below was used to solve the problem for clique sizes up to 10. In anticlust,
# I pasted the results into the code, rather than recomputing them each time.
# I suppose there are faster ways to solve this problem than with brute force (in R).
# It crashed my R session with a clique size of 11. I guess we could do better
# with a smarter algorithm / better implementation.
# For n_clique = 10, the subset, which has to be searched is
# [1]  1  1  1  1  1  1  1  1  1  1  2  2  2  2  2  3  3  3  4  4  5  5  6  7  8  9 10
# (given via function `get_set_for_subset_problem()`)
# and it already has size of 27, making complete enumeration (which is done here) hard. 
# For n_clique = 11, it is `get_set_for_subset_problem(11)`
#  [1]  1  1  1  1  1  1  1  1  1  1  1  2  2  2  2  2  3  3  3  4  4  5  5  6  7  8  9 10 11 
# which has length 29, crashing my personal computer. 

library(here) # for writing to file in the "correct" directory, not used for computation

valid_sums_clique <- function(n_clique) { # returns a list of combinations that sum to n_clique
  set <- get_set_for_subset_problem(n_clique)
  subset_sums <- subset_sum(set, n_clique)
  subset_sums[sapply(subset_sums, function(x) !all(x == 1))] # exclude the entry with only Ones
}

get_set_for_subset_problem <- function(n_clique) {
  init <- rep(1:n_clique, n_clique:1)
  init_list <- mapply(FUN = rep, 1:n_clique, n_clique:1)
  cumsums <- lapply(1:n_clique, function(x) cumsum(init[init == x]))
  valid_cumsums <- lapply(cumsums, function(x) x[x <= n_clique])
  indices <- lapply(valid_cumsums, function(x) 1:length(x))
  set <- unlist(mapply("[", init_list, indices)) # set of numbers that can add up to n_clique
  set <- set[!is.na(set)] # some cleanup
  set
}

# all combinations of a vector that have sum = n_clique
subset_sum <- function(x, n_clique) {
  l <- list()
  for (i in 1:(length(x)-1)) {
    l[[i]] <- combn(x, i+1)
  }
  all_sums <- lapply(l, function(x) colSums(x))
  indices_good <- lapply(l, function(x) which(colSums(x) == n_clique))
  nl <- list()
  for (i in seq_along(l)) {
    values <- l[[i]][, indices_good[[i]], drop = FALSE]
    nl[[i]] <- values
  }
  nl <- unlist(lapply(nl, function(x) as.list(as.data.frame(x))), recursive = FALSE) #just some cleanup
  unname(nl[!duplicated(nl)])
}


## Generate the sums and write to file:
subset_sums <- list()

for (i in 3:10) {
  subset_sums[[as.character(i)]] <- valid_sums_clique(i)
}

write(deparse(dput(subset_sums)), file = here("2PML", "sums.txt")) # I include this output in the code of anticlust
