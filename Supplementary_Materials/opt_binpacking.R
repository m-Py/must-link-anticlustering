
K <- 3L
N <- 9L

c <- as.integer(N/3)
m <- K

n_batches <- 3
n <- 5
capacity <- 5
weights <- sample(3, size = n, replace = TRUE)

get_col_names <- function(n_batches, n) {
  combs <- expand.grid(1:n_batches, 1:n)
  paste0("b", combs[,1], "_item_", combs[, 2], "_")
}

variables <- get_col_names(n_batches, n)

constraints2 <- matrix(0, ncol = length(variables), nrow = n)
colnames(constraints2) <- variables
for (i in 1:nrow(constraints2)) {
  constraints2[i, grepl(paste0("item_", i, "_"), variables)] <- 1
}



constraints1 <- matrix(0, ncol = length(variables), nrow = n_batches)
colnames(constraints1) <- variables
for (i in 1:nrow(constraints1)) {
  constraints1[i, grepl(paste0("b", i, "_"), variables)] <- weights
}


dirs <- c(rep("==", nrow(constraints2)), rep("<=", n_batches))
rhs <- c(rep(1, nrow(constraints2)), rep(capacity, n_batches))

tt <- lpSolve::lp(
  "max", # min
  rep(weights, each = n_batches), 
  rbind(constraints2, constraints1), 
  dirs,
  rhs,
  all.bin = TRUE
)

tt$status

variables[tt$solution == 1]
weights
