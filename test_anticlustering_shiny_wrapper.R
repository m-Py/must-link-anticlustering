

library(anticlust)
library(tableone)
source("anticlustering_wrapper_shiny.R")

generate_categorical_data <- function(N, M, P) {
  data <- data.frame(matrix(sample(P, replace = TRUE, size = N*M), ncol = M))
  as.data.frame(lapply(data, as.factor))
}

# Small data 

N <- 100
M_numeric <- 3  # number of numeric features
M_categorical <- 3
numeric_vars <- matrix(rnorm(N*M_numeric), ncol = M_numeric)
categorical_vars <- generate_categorical_data(N, M_categorical, 3)

data <- data.frame(numeric_vars, categorical_vars)
head(data)

data$Groups <- anticlustering_shiny(
  numeric_vars = numeric_vars,
  categorical_vars = categorical_vars, 
  K = 5
)

print(tableone::CreateTableOne(data = data, vars = colnames(data)[1:(M_numeric+M_categorical)], strata = "Groups"), smd = TRUE)


# Rather Small data 

N <- 200
M_numeric <- 3  # number of numeric features
M_categorical <- 3
numeric_vars <- matrix(rnorm(N*M_numeric), ncol = M_numeric)
categorical_vars <- generate_categorical_data(N, M_categorical, 3)

data <- data.frame(numeric_vars, categorical_vars)
head(data)

data$Groups <- anticlustering_shiny(
  numeric_vars = numeric_vars,
  categorical_vars = categorical_vars, 
  K = 5
)

print(tableone::CreateTableOne(data = data, vars = colnames(data)[1:(M_numeric+M_categorical)], strata = "Groups"), smd = TRUE)


# Smallish data 

N <- 600
M_numeric <- 3  # number of numeric features
M_categorical <- 3
numeric_vars <- matrix(rnorm(N*M_numeric), ncol = M_numeric)
categorical_vars <- generate_categorical_data(N, M_categorical, 3)

data <- data.frame(numeric_vars, categorical_vars)
head(data)

data$Groups <- anticlustering_shiny(
  numeric_vars = numeric_vars,
  categorical_vars = categorical_vars, 
  K = 5
)

print(tableone::CreateTableOne(data = data, vars = colnames(data)[1:(M_numeric+M_categorical)], strata = "Groups"), smd = TRUE)


# "large" data:

N <- 2000
M_numeric <- 3  # number of numeric features
M_categorical <- 3
numeric_vars <- matrix(rnorm(N*M_numeric), ncol = M_numeric)
categorical_vars <- generate_categorical_data(N, M_categorical, 3)

data <- data.frame(numeric_vars, categorical_vars)
head(data)

data$Groups <- anticlustering_shiny(
  numeric_vars = numeric_vars,
  categorical_vars = categorical_vars, 
  K = 5
)

print(tableone::CreateTableOne(data = data, vars = colnames(data)[1:(M_numeric+M_categorical)], strata = "Groups"), smd = TRUE)


# Very large data:

N <- 12000
M_numeric <- 3  # number of numeric features
M_categorical <- 3
numeric_vars <- matrix(rnorm(N*M_numeric), ncol = M_numeric)
categorical_vars <- generate_categorical_data(N, M_categorical, 3)

data <- data.frame(numeric_vars, categorical_vars)
head(data)

data$Groups <- anticlustering_shiny(
  numeric_vars = numeric_vars,
  categorical_vars = categorical_vars, 
  K = 5
)

print(tableone::CreateTableOne(data = data, vars = colnames(data)[1:(M_numeric+M_categorical)], strata = "Groups"), smd = TRUE)



