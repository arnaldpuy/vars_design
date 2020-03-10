
library(tidyverse)

# FUNCTION TO CREATE VARS MATRICES --------------------------------------------

vars_matrices <- function(N, params, h) {
  out <- center <- sections <- A <- B <- AB <- X <- out <- list()
  mat <- randtoolbox::sobol(n = N, dim = length(params))
  for(i in 1:nrow(mat)) {
    center[[i]] <- mat[i, ]
    sections[[i]] <- sapply(center[[i]], function(x) {
      all <- seq(x %% h, 1, h)
      non.zeros <- all[all!= 0] # Remove zeroes
      })
    B[[i]] <- sapply(1:ncol(mat), function(x) sections[[i]][, x][!sections[[i]][, x] %in% center[[i]][x]])
    A[[i]] <- matrix(center[[i]], nrow = nrow(B[[i]]), ncol = length(center[[i]]), byrow = TRUE)
    X[[i]] <- rbind(A[[i]], B[[i]])
    for(j in 1:ncol(A[[i]])) {
      AB[[i]] <- A[[i]]
      AB[[i]][, j] <- B[[i]][, j]
      X[[i]] <- rbind(X[[i]], AB[[i]])
    }
    AB[[i]] <- X[[i]][(2 * nrow(B[[i]]) + 1):nrow(X[[i]]), ]
    out[[i]] <- rbind(unname(center[[i]]), AB[[i]])
  }
  return(do.call(rbind, out))
}

# Function to cut by size
CutBySize <- function(m, block.size, nb = ceiling(m / block.size)) {
  int <- m / nb
  upper <- round(1:nb * int)
  lower <- c(1, upper[-nb] + 1)
  size <- c(upper[1], diff(upper))
  cbind(lower, upper, size)
}

# Function to compute variogram
variogram <- function(x) 1/ (2 * nrow(x)) * sum(x[, 1] - x[, 2]) ^ 2

# Function to compute Total variance
variance_f <- function(x) {
  f0 <- 1 / (2 * nrow(x)) * sum(x[, 1] + x[, 2])
  VY <- 1 / (2 * nrow(x) - 1) * sum((x[, 1] - f0) ^ 2 + (x[, 2] - f0) ^ 2)
  return(VY)
}

# SETTINGS --------------------------------------------------------------------
N <- 2
params <- paste("X", 1:8, sep = "")
h <- 0.1

# Create STAR-VARS
mat <- vars_matrices(N = N, params = params, h = h)

# Check that number of rows of mat is identical to what is expected according to Razavi & Gupta
nrow(mat) == N * (length(params) * ((1 / h) -1) + 1)

# Compute model output
Y <- sensobol::sobol_Fun(mat)

# Index for star centers
n.cross.points <- length(params) * ((1 / h) - 1) + 1

index.centers <- seq(1, nrow(mat), n.cross.points)
mat.nocenters <- matrix(Y[-index.centers], ncol = N)
mat.centers <- matrix(Y[index.centers], 
                      nrow = nrow(mat.nocenters) + length(params),
                      ncol = N, 
                      byrow = TRUE)
location.centers <- seq(1, nrow(mat.nocenters), 1 / h)
mat.centers[-location.centers, ] <- mat.nocenters

indices <- CutBySize(nrow(mat.centers), nb = length(params))

out <- list()
for(i in 1:nrow(indices)) {
  out[[i]] <- mat.centers[indices[i, "lower"]:indices[i, "upper"], ]
}

tmp <- lapply(1:length(params), function(x) 
  lapply(1:ncol(out[[x]]), function(y) 
    t(combn(out[[x]][, y], 2))))

prove2 <- lapply(tmp, function(x) do.call(rbind, x))


# Compute variogram
variogr <- unlist(lapply(prove2, function(x) 1/ (2 * nrow(x)) * sum(x[, 1] - x[, 2]) ^ 2))

# Covariogram
covariog <- unlist(lapply(prove2, function(x) 1 / (2 * nrow(x)) * cov(x[, 1], x[, 2])))

# Total variance
VY <- variance_f(do.call(rbind, prove2))

# Final value
(variogr + covariog) / VY



# COMPARISON WITH JANSEN ------------------------
N <- 1000
params <- paste("X", 1:3, sep = "")
mat <- sensobol::sobol_matrices(N = N, params = params)
Y <- sensobol::ishigami_Fun(mat)

k <- length(params)
m <- matrix(Y, nrow = N)

Y_A <- m[, 1]
Y_B <- m[, 2]
Y_AB <- m[, -c(1, 2)]
f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
VY <- 1 / (2 * N - 1) * sum((Y_A - f0) ^ 2 + (Y_B - f0) ^ 2)

Vi <- 1 / (2 * N) * Rfast::colsums((Y_A - Y_AB) ^ 2) 

Vi/VY

lapply(prove2, function(x) (1 / 2 * var(x[, 1] - x[, 2]) ^ 2)  + cov(x[, 1], x[, 2]))









f0 <- (1 / N) * sum(m[, 1] * m[, 2])
VY <- 1 / (2 * N - 1) * sum(m[, 1] ^ 2 + m[, 2] ^ 2) - f0




cross.section.mat <- matrix(Y[-index.centers], ncol = N)

indices <- CutBySize(nrow(cross.section.mat), nb = length(params))

# List of Matrices: each slot is a model input. The first row in each
# slot is the star centers. Each column is a star.
out <- list()
for(i in 1:nrow(indices)) {
  out[[i]] <- cross.section.mat[indices[i, "lower"]:indices[i, "upper"], ]
  out[[i]] <- rbind(Y[index.centers], out[[i]])
}












mat.without.centers <- matrix(Y[-index.centers], ncol = N)






mat.without.centers <- matrix(Y[-index.centers], ncol = N)

location.centers <- seq(1, nrow(mat.without.centers), 1 / h)

new.mat <- matrix(Y[index.centers], 
                  nrow = nrow(mat.without.centers) + length(params), 
                  ncol = N, 
                  byrow = TRUE)

new.mat[-index.centers, ] <- mat.without.centers












mat.without <- matrix(Y[-index.centers], ncol = N)

new.mat <- matrix(0, nrow = nrow(mat.without) + length(index.centers), 
                  ncol = ncol(mat.without))

new_mat[-index.centers,] <- Y[index.centers]   
new_mat


matrix(Y, ncol = N)

Y.without.centers <- Y[-index.centers]

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

lapply(chunk2(Y.without.centers, n = N), 
       function(x) chunk2(x, n = length(params)))






mat <- matrix(1:18,6)
vec <- c(2, 5, 6)

# New matrix 'new_mat' with all zeros, 
# No. of rows = original matrix rows + number new rows to be added
new_mat <- matrix(Y[index.centers], 
                  nrow = nrow(mat.without) + length(index.centers), 
                  ncol = N, 
                  byrow = TRUE)

new_mat[-index.centers, ] <- mat.without











# 'new_mat' rows getting filled with `mat` values
new_mat[-vec,] <- mat   
new_mat


# Index for star centers
index.centers <- seq(1, nrow(mat), n.cross.points)

mat.without <- matrix(Y[-index.centers], ncol = N)

new.mat <- matrix(0, 
                  nrow = nrow(mat.without) + length(index.centers), 
                  ncol = ncol(mat.without))

new.mat[-index.centers] <- mat.without








# Matrix with only the cross-section points 
cross.section.mat <- matrix(Y[-index.centers], ncol = N)

indices <- CutBySize(nrow(cross.section.mat), nb = length(params))

# List of Matrices: each slot is a model input. The first row in each
# slot is the star centers. Each column is a star.
out <- list()
for(i in 1:nrow(indices)) {
  out[[i]] <- cross.section.mat[indices[i, "lower"]:indices[i, "upper"], ]
  out[[i]] <- rbind(Y[index.centers], out[[i]])
}










# Each slot is a parameter, all stars are row-binded
prove <- lapply(1:length(params), function(x) 
  lapply(1:ncol(out[[x]]), function(y) 
  t(combn(out[[x]][, y], 2))))

prove2 <- lapply(prove, function(x) do.call(rbind, x))

sum(sapply(1:10, function(x) N * (( 1 / h) - x)))

# Compute variogram
variogr <- unlist(lapply(prove2, variogram))

# Covariogram
covariog <- unlist(lapply(prove2, function(x) 1 / (2 * nrow(x)) * cov(x[, 1], x[, 2])))

# Total variance
VY <- unlist(lapply(prove2, variance_f))

# Final value
(variogr - covariog) / VY

























pairs.points <- t(combn(1:(((1 / h) - 1) + 1), 2))

h.steps <- abs(pairs.points[, 1] / 10 - pairs.points[, 2] / 10)


cbind(pairs.points, h.steps)






# Compute variogram
vars <- unlist(lapply(prove2, variogram))

# compute covariogram
covars <- unlist(lapply(prove2, function(x) cov(x[, 1], x[, 2])))

# compute variance
VY <- unlist(lapply(prove2, function(x) var(x[, 1], x[, 2])))











(vars-covars) / VY



do.call(rbind, prove)


var







prove <- list()
for(j in 1:ncol(out[[1]])) {
  prove[[j]] <- do.call(rbind, combn(out[[1]][, j], 2, simplify = FALSE))
}



prove <- list()
for(j in 1:ncol(out[[1]])) {
  prove[[j]] <- do.call(rbind, combn(out[[1]][, j], 2, simplify = FALSE))
}

sum(sapply(1:10, function(x) N * (( 1 / h) - x)))


input3 <- do.call(rbind, prove)

((1 / (2 * nrow(input3)) * sum(input3[, 1] - input3[, 2]) ^ 2) - cov(x = input3[, 1], y = input3[, 2])) / 
  var(x = input3[, 1], y = input3[, 2])







lapply(lapply(out, c), function(x) combn(x, 2, simplify = FALSE))



do.call(rbind, combn(out[[1]], 2, simplify = FALSE))












matrix(Y[-index.centers], ncol = N)






mapply(rbind, t(Y[index.centers]), Y[-index.centers])


A <- randtoolbox::sobol(n = N, dim = length(params))
 cross.sections <- B <- A <- X <- AB <- list()
for(i in 1:nrow(mat)) {
  cross.sections[[i]] <- sapply(mat[i, ], function(x) seq(x %% h, 1, h))
  B[[i]] <- sapply(1:ncol(mat), function(x) 
    cross.sections[[i]][, x][!cross.sections[[i]][, x] %in% mat[i, ][x]])
  A[[i]] <- matrix(mat[i, ], nrow = nrow(B[[i]]), ncol = length(mat[i, ]), byrow = TRUE)
  X[[i]] <- rbind(A[[i]], B[[i]])
  for(j in 1:ncol(A[[i]])) {
    AB[[i]] <- A[[i]]
    AB[[i]][, j] <- B[[i]][, j]
    X[[i]] <- rbind(X[[i]], AB[[i]])
  }
}

cross.sections

all(unlist(lapply(cross.sections, nrow)))
all(unlist(lapply(B, nrow)))
all(unlist(lapply(A, nrow)))
all(unlist(lapply(X, nrow)))




mat <- vars_matrices(N = N, params = params, h = h)





nrow(mat)
N * (length(params) * ((1 / h) - 1) + 1)



# PROVE-----------------------
N <- 8
params <- paste("X", 1:3, sep = "")
h <- 0.1
mat <- randtoolbox::sobol(n = N, dim = length(params))
center <- mat[N, ]


sections <- sapply(center, function(x) { 
  all <- seq(x %% h, 1, h)
  non.zeros <- all[all!=0]
})


sapply(center, function(x) seq(x %% h, 1, h))


sapply(1:ncol(mat), function(x) sections[, x][!sections[, x] %in% center[x]])


























Y <- sensobol::ishigami_Fun(mat)


rows.star <- length(params) * ((1 / h) - 1) + 1




Y.mat <- matrix(Y, ncol = N)

star.centers <- Y.mat[1, ]
Y.mat.without.centers <- Y.mat[-1, ]



indices <- CutBySize(nrow(Y.mat.without.centers), nb = length(params))




lapply(list(1:10, c(1, 11:19)), function(x) combn(x, 2, simplify = FALSE))




CutBySize <- function(m, block.size, nb = ceiling(m / block.size)) {
  int <- m / nb
  upper <- round(1:nb * int)
  lower <- c(1, upper[-nb] + 1)
  size <- c(upper[1], diff(upper))
  cbind(lower, upper, size)
}

da <- CutBySize(nrow(Y.mat[-1, ]), nb = length(params))

out <- list()
for(i in 1:nrow(da)) {
  out[[i]] <- Y.mat[da[i, "lower"]:da[i, "upper"], ]
}






apply(da, 1, function(x) da["lower":"upper", ])




rows.star <- length(params) * ((1 / h) - 1) + 1

cross.points <- (1 / 0.1) - 1












A <- randtoolbox::sobol(N, length(params))


pairs <- matrix(Y[1:10], ncol = 3)

pairs2 <- do.call(rbind, combn(Y[1:10], 2, simplify = FALSE))

variogram <- 1 / (2 * nrow(pairs2)) * sum(pairs2[, 1] - pairs2[, 2]) ^ 2
covariogram <- cov(pairs2[, 1], pairs2[, 2])


(variogram - covariogram) / var(pairs2[, 1], pairs2[, 2])

combn(1:10)


var(pairs2[, 1], pairs2[, 2])

sum(sapply(1:10, function(x) N * (( 1 / h) - x)))



cov(mat[, 1], mat[, 2])


matrix(Y[1:10], nrow = length(params) * (( 1 / h) - 1) + 1)










N * ((1 / h) - 1)
N * ((1 / h) - 2)
N * ((1 / h) - 3)

da <- sapply(1:3, function(x) N * ((1 / h) - x))

sum(da)


combn(mat[, 1], 2)

N * ((1 / h) - 2)


vars_matrices <- function(N, params, h) {
  out <- center <- sections <- A <- B <- AB <- X <- out <- list()
  mat <- randtoolbox::sobol(n = N, dim = length(params))
  for(i in 1:nrow(mat)) {
    center[[i]] <- mat[i, ]
    sections[[i]] <- sapply(center[[i]], function(x) {
      all <- seq(x %% h, 1, h)
      non.zeros <- all[all!= 0] # Remove zeroes
    })
    B[[i]] <- sapply(1:ncol(mat), function(x) sections[[i]][, x][!sections[[i]][, x] %in% center[[i]][x]])
    A[[i]] <- matrix(center[[i]], nrow = nrow(B[[i]]), ncol = length(center[[i]]), byrow = TRUE)
    X[[i]] <- rbind(A[[i]], B[[i]])
    for(j in 1:ncol(A[[i]])) {
      AB[[i]] <- A[[i]]
      AB[[i]][, j] <- B[[i]][, j]
      X[[i]] <- rbind(X[[i]], AB[[i]])
    }
    AB[[i]] <- X[[i]][(2 * nrow(B[[i]]) + 1):nrow(X[[i]]), ]
    out[[i]] <- rbind(mat, AB[[i]])
  }
  return(do.call(rbind, out))
}














nrow(out)
N * (length(params) * (( 1 / h) -1) + 1)




length(params) * (( 1 / h) - 1)
















N <- 10
params <- paste("X", 1:3, sep = "")
mat <- randtoolbox::sobol(N, length(params))

h <- 0.3
center <- mat[1, ]

seq(h %% 0.1, 1, 0.1)

sections <- sapply(center, function(x) { 
  all <- seq(x %% h, 1, h)
  non.zeros <- all[all!=0]
})

B <- sapply(1:ncol(mat), function(x) sections[, x][!sections[, x] %in% center[x]])

A <- matrix(center, nrow = nrow(B), ncol = length(center), byrow = TRUE)



length(params) * (length(params) - 1) / 2

X <- rbind(A, B)
for(j in 1:ncol(A)) {
  AB <- A
  AB[, j] <- B[, j]
  X <- rbind(X, AB)
}

AB <- X[(2 * nrow(B) + 1):nrow(X), ]
out <- rbind(unname(center), AB)

length(params) * ((1 / h) - 1)






N * (length(params) * ((1 / h) - 1) + 1)

N <- 70
params <- paste("X", 1:120, sep ="")
h <- 0.2

mat <- randtoolbox::sobol(N, length(params))

vars_matrices(N = N, params = params, h = h)







N * (length(params) * ((1 / h) - 1) + 1)





da <- vars_matrices(N = N, params = params, h = h)
nrow(da)

N * (length(params) * ((1 / h)) + 1)

}
out <- center <- sections <- A <- B <- AB <- X <- out <- list()
for(i in 1:nrow(mat)) {
  center[[i]] <- mat[i, ]
  sections[[i]] <- sapply(center[[i]], function(x) seq(x %% h, 1, h))
  B[[i]] <- sapply(1:ncol(mat), function(x) sections[[i]][, x][!sections[[i]][, x] %in% center[[i]][x]])
  A[[i]] <- matrix(center[[i]], nrow = nrow(B[[i]]), ncol = length(center[[i]]), byrow = TRUE)
  X[[i]] <- rbind(A[[i]], B[[i]])
  for(j in 1:ncol(A[[i]])) {
    AB[[i]] <- A[[i]]
    AB[[i]][, j] <- B[[i]][, j]
    X[[i]] <- rbind(X[[i]], AB[[i]])
  }
  AB[[i]] <- X[[i]][(2 * nrow(B[[i]]) + 1):nrow(X[[i]]), ]
  out[[i]] <- rbind(unname(center[[i]]), AB[[i]])
  return(do.call(rbind, out))
}


