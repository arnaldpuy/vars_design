



















# FUNCTION TO CREATE STAR-VARS ------------------------------------------------

vars_matrices_NEW <- function(N, params, h) {
  out <- center <- sections <- A <- B <- AB <- X <- out <- list()
  mat <- randtoolbox::sobol(n = N, dim = length(params))
  for(i in 1:nrow(mat)) {
    center[[i]] <- mat[i, ]
    sections[[i]] <- sapply(mat[i, ], function(x) {
      all <- seq(x %% h, 1, h)
      non.zeros <- all[all!= 0] # Remove zeroes
    })
    A[[i]] <- matrix(mat[i, ], nrow = nrow(sections[[i]]), 
                     ncol = ncol(mat), byrow = TRUE)
    X[[i]] <- rbind(A[[i]], sections[[i]])
    for(j in 1:ncol(A[[i]])) {
      AB[[i]] <- A[[i]]
      AB[[i]][, j] <- sections[[i]][, j]
      X[[i]] <- rbind(X[[i]], AB[[i]])
    }
    AB[[i]] <- X[[i]][(2 * nrow(sections[[i]]) + 1):nrow(X[[i]]), ]
  }
  star.vars <- do.call(rbind, AB)
  tmp <- lapply(1:N, function(x) 
    rowSums(star.vars == mat[x, ][col(star.vars)]) == ncol(star.vars)) 
  indices.star.output <- unlist(lapply(1:N, function(x) which(tmp[[x]] == TRUE)[1]))
  return(list(star.vars, indices.star.output))
}

pairs_tmp <- function(x, h) {
  da <- list()
  for(i in 1:((1/h) - 1)) {
    da[[i]] <- c(x[i], x[i+1])
  }
  out <- do.call(rbind, da)
  return(out)
}

pairs_fun <- function(x, h) {
  da <- list()
  for(j in 1:ncol(x)) {
    da[[j]] <- pairs_tmp(x[, j], h = h)
  }
  return(da)
}

vars_ti_NEW <- function(Y, N, indices.var, params, h) {
  mat <- matrix(Y, ncol = N)
  indices <- CutBySize(nrow(mat), nb = length(params))
  out <- list()
  for(i in 1:nrow(indices)) {
    out[[i]] <- mat[indices[i, "lower"]:indices[i, "upper"], ]
  }
 d <- lapply(out, function(x) pairs_fun(x, h = h))
 variogr <- unlist(lapply(d, function(x) lapply(x, function(y) 
   mean(0.5 * (y[, 1] - y[, 2]) ^ 2))) %>%
     lapply(., function(x) do.call(rbind, x)) %>%
     lapply(., mean))
covariogr <- unlist(lapply(d, function(x)
  lapply(x, function(y) cov(y[, 1], y[, 2]))) %>%
    lapply(., function(x) Rfast::colmeans(do.call(rbind, x))))
 VY <- var(Y[indices.var])
 output <- (variogr + covariogr) / VY
 return(output)
}




variogr <- unlist(lapply(out, function(x) lapply(x, function(y) 
  mean(0.5 * (y[, 1] - y[, 2]) ^ 2))) %>%
    lapply(., function(x) do.call(rbind, x)) %>%
    lapply(., mean))


N <- 100
h <- 0.1
params <- paste("X", 1:3, sep = "")
mat <- vars_matrices_NEW(N = N, params = params, h = h)
nrow(mat[[1]])
N * (length(params) * ((1 / h) - 1) + 1) + N * length(params) - N
mat2 <- vars_matrices(N = N, params = params, h = h)
N * (length(params) * ((1 / h) - 1) + 1)
nrow(mat2)

(N * length(params)) / h






Y <- sensobol::ishigami_Fun(mat[[1]])
output <- vars_ti_NEW(Y = Y, N = N, params = params, h = h, 
                      indices.var = mat[[2]])
round(output, digits = 3)
  
params <- paste("X", 1:8, sep = "")
mat <- vars_matrices_NEW(N = N, params = params, h = h)
Y <- sensobol::sobol_Fun(mat[[1]])
output <- vars_ti_NEW(Y = Y, N = N, params = params, h = h, 
                      indices.var = mat[[2]])
round(output, digits = 3)

params <- paste("X", 1:20, sep = "")
mat <- vars_matrices_NEW(N = N, params = params, h = h)
Y <- sensitivity::morris.fun(mat[[1]])
output <- vars_ti_NEW(Y = Y, N = N, params = params, h = h, 
                      indices.var = mat[[2]])
round(output, digits = 3)

 


# DEFINE THE SETTINGS FOR A STAR-VARS SAMPLE MATRIX ---------------------------

N <- 200 # Star centers
h <- 0.2 # h step
params <- paste("X", 1:6, sep = "")

mt <- vars_matrices(N = N, params = params, h = h)

mt[, 1] <- -sin(pi * mt[, 1]) - 0.3 * sin(3.33 * pi * mt[, 1])
mt[, 2] <- -0.76 * sin(pi * (mt[, 2] - 0.2)) - 0.315
mt[, 3] <- -0.12 * sin(1.05 * pi * (mt[, 3] - 0.2)) - 0.02 * sin(95.24 * pi * mt[, 3]) - 0.96 
mt[, 4] <- -0.12 * sin(1.05 * pi * (mt[, 4] - 0.2)) - 0.96
mt[, 5] <- -0.05 * sin(pi * (mt[, 5] - 0.2)) - 1.02
mt[, 6] <-  -1.08

Y <- Rfast::rowsums(mt)
output <- vars_ti(Y = Y, N = N, params = params, h = h)
output



da <- vars_matrices_NEW(N = N, params = params, h = h)
mt <- da[[1]]

mt[, 1] <- -sin(pi * mt[, 1]) - 0.3 * sin(3.33 * pi * mt[, 1])
mt[, 2] <- -0.76 * sin(pi * (mt[, 2] - 0.2)) - 0.315
mt[, 3] <- -0.12 * sin(1.05 * pi * (mt[, 3] - 0.2)) - 0.02 * sin(95.24 * pi * mt[, 3]) - 0.96 
mt[, 4] <- -0.12 * sin(1.05 * pi * (mt[, 4] - 0.2)) - 0.96
mt[, 5] <- -0.05 * sin(pi * (mt[, 5] - 0.2)) - 1.02
mt[, 6] <-  -1.08

Y <- Rfast::rowsums(mt)
output <- vars_ti_NEW(Y = Y, N = N, params = params, h = h, indices.var = mt[[2]])
output


