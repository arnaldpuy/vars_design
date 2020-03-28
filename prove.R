
devtools::install_github("arnaldpuy/sensobol", build_vignettes = TRUE)
library(sensobol)








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

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

# SETTINGS --------------------------------------------------------------------
N <- 20
params <- paste("X", 1:3, sep = "")
h <- 0.1

# Create STAR-VARS
mat <- vars_matrices(N = N, params = params, h = h)

# Check that number of rows of mat is identical to what is expected according to Razavi & Gupta
nrow(mat) == N * (length(params) * ((1 / h) -1) + 1)

# Compute model output
Y <- sensobol::ishigami_Fun(mat)

# Calculate number of points in star
n.points.star <- length(params) * ((1 / h) - 1) + 1

# Locate star centers
location.centers <- seq(1, length(Y), n.points.star)

# Extract star centers
star.centers <- Y[location.centers]

da <- matrix(Y[-location.centers], ncol = N)

# Split matrix
indices <- CutBySize(nrow(da), nb = length(params))

out <- list()
for(i in 1:nrow(indices)) {
  out[[i]] <- da[indices[i, "lower"]:indices[i, "upper"], ]
  out[[i]] <- rbind(out[[i]], Y[location.centers])
}

pairs.points <- lapply(1:length(params), function(x) 
  lapply(1:ncol(out[[x]]), function(y) 
    t(combn(out[[x]][, y], 2))))
pairs.points <- lapply(pairs.points, function(x) do.call(rbind, x))

# Compute VY
VY <- var(Y[location.centers])

d <- unlist(lapply(pairs.points, function(x) 
  (1 / (2 * nrow(x)) * sum(x[, 1] - x[, 2]) ^ 2)))

d / VY









# Remove star centers from output vector and bind with star centers
mt2 <- cbind(matrix(Y[-location.centers]), rep(star.centers, each = n.points.star-1))

# Split matrix
indices <- CutBySize(nrow(mt2), nb = length(params))

out <- list()
for(i in 1:nrow(indices)) {
  out[[i]] <- mt2[indices[i, "lower"]:indices[i, "upper"], ]
}




















n.points.star <- length(params) * ((1 / h) - 1) + 1
n.cross.points <- (1 / h) - 1
index.centers <- seq(1, length(Y), n.points.star)

matrix(Y, ncol = N)








# Index for star centers
n.cross.points <- length(params) * ((1 / h) - 1) + 1
index.centers <- seq(1, length(Y), n.cross.points)
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

VY <- var(Y[index.centers])
prove2 <- lapply(tmp, function(x) do.call(rbind, x))

unlist(lapply(prove2, function(x) 
  (1/ (2 * nrow(x)) * sum(x[, 1] - x[, 2]) ^ 2) + cov(x[, 1], x[, 2]))) / VY










split(mat.centers, rep(1:ceiling((((1 / h) - 1) + 1)/length(params)), each=((1 / h) - 1) + 1, length.out=length(params)))








index.centers <- seq(1, length(Y), n.cross.points)
mat.nocenters <- matrix(Y[-index.centers], ncol = N)
mat.centers <- matrix(Y[index.centers], 
                      nrow = nrow(mat.nocenters) + length(params),
                      ncol = N, 
                      byrow = TRUE)
location.centers <- seq(1, nrow(mat.nocenters), 1 / h)
mat.centers[-location.centers, ] <- mat.nocenters




matrix(Y, ncol = N)



index.centers <- seq(1, length(Y), n.cross.points)
mat.nocenters <- matrix(Y[-index.centers], ncol = N)
mat.centers <- matrix(Y[index.centers], 
                      nrow = nrow(mat.nocenters) + length(params),
                      ncol = N, 
                      byrow = TRUE)
location.centers <- seq(1, nrow(mat.nocenters), 1 / h)
mat.centers[-location.centers, ] <- mat.nocenters










matrix(Y[-index.centers], nrow = N * length((params)))

matrix(a, ncol = 3, byrow = FALSE)

















matrix(Y[-index.centers], ncol = length(params))



length(Y[-index.centers]) / 3


do.call(cbind, chunk2(Y, N))




cbind(mat, Y)



Y[index.centers]

do.call(cbind, chunk2(Y[-index.centers], length(params)))







N <- 60
params <- paste("X", 1:3, sep = "")
h <- 0.1

# Create STAR-VARS
mat <- vars_matrices(N = N, params = params, h = h)

# Check that number of rows of mat is identical to what is expected according to Razavi & Gupta
nrow(mat) == N * (length(params) * ((1 / h) -1) + 1)

# Compute model output
Y <- sensobol::ishigami_Fun(mat)
n.points.star <- length(params) * ((1 / h) - 1) + 1
n.cross.points <- (1 / h) - 1
index.centers <- seq(1, length(Y), n.points.star)



VY <- var(Y[index.centers])

mm <- matrix(Y[-index.centers])
p <- cbind(mm, rep(Y[index.centers], each = nrow(mm) / N))
dt <- data.table::data.table(p)[
  , parameters:= rep(params, each = (1 / h) - 1, times = N)] 



da <- dt[, ((1 / (2 * .N) * sum(V1 - V2) ^ 2) + cov(V1, V2)), parameters]

da$V1 / VY



n <- 1000
params <- paste("X", 1:3, sep = "")
A <- sensobol::sobol_matrices(N = n, params = params)
Y <- sensobol::ishigami_Fun(A)
ind <- sensobol::sobol_indices(Y = Y, N = n, params = params)





indices <- CutBySize(nrow(p), nb = ((1 / h) - 1))
out <- list()
for(i in 1:nrow(indices)) {
  out[[i]] <- p[indices[i, "lower"]:indices[i, "upper"], ]
}








split(p, rep(1:ncol(p), each = nrow(mm) / N))



CutBySize(nrow(p), nb = length(params))




54/3





n.cross.points <- length(params) * ((1 / h) - 1) + 1
index.centers <- seq(1, length(Y), n.cross.points)
mat.nocenters <- matrix(Y[-index.centers], ncol = N)



Y[index.centers]





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
pairs.points <- lapply(1:length(params), function(x) 
  lapply(1:ncol(out[[x]]), function(y) 
    t(combn(out[[x]][, y], 2))))
pairs.points <- lapply(pairs.points, function(x) do.call(rbind, x))



lapply(pairs.points, function(x) 1 / (2 * nrow(x)) * colsums(x[, 2] - x[, 1]) ^ 2 )


colsum

























# SETTINGS --------------------------------------------------------------------
N <- 2
params <- paste("X", 1:3, sep = "")
h <- 0.2

mat.prove <- randtoolbox::sobol(n = N, dim = length(params))
# Create STAR-VARS
mat <- vars_matrices(N = N, params = params, h = h)

# Check that number of rows of mat is identical to what is expected according to Razavi & Gupta
nrow(mat) == N * (length(params) * ((1 / h) -1) + 1)

# Compute model output
Y <- sensobol::ishigami_Fun(mat)

n.cross.points <- length(params) * ((1 / h) - 1) + 1
index.centers <- seq(1, length(Y), n.cross.points)
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
pairs.points <- lapply(1:length(params), function(x) 
  lapply(1:ncol(out[[x]]), function(y) 
    t(combn(out[[x]][, y], 2))))
pairs.points <- lapply(pairs.points, function(x) do.call(rbind, x))
# Compute variogram
variogr <- unlist(lapply(pairs.points, function(x) 1 / 2 * (var(x[, 1] - x[, 2]) ^ 2) + cov(x[, 1], x[, 2])))
# Covariogram
covariog <- unlist(lapply(pairs.points, function(x) 1 / (2 * nrow(x)) * cov(x[, 1], x[, 2])))
# Total variance
VY <- variance_f(do.call(rbind, pairs.points))
final <- (variogr + covariog) / VY




STi <- 1 / (2 * N) * Rfast::colsums((Y_A - Y_AB) ^ 2) / VY





cols_f <- function(x) {
  pairs <- c(x[1], rep(x[2:(length(x) - 1)], each = 2), tail(x, 1))
  out <- matrix(pairs, nrow = length(pairs) / 2, byrow = TRUE)
  return(out)
}

pairs.points <- lapply(1:length(params), function(x) 
  lapply(1:ncol(out[[x]]), function(y) 
    cols_f(out[[x]][, y])))

pairs.points <- lapply(pairs.points, function(x) do.call(rbind, x))

# Compute variogram
variogr <- unlist(lapply(pairs.points, function(x) 1 / 2 * var(x[, 1] - x[, 2]) ^  2))

pp <- pairs.points[[1]]

cov(pp[, 1], pp[, 2])

# Covariogram
covariog <- unlist(lapply(pairs.points, function(x) cov(x[, 1], x[, 2])))
# Total variance
VY <- variance_f(do.call(rbind, pairs.points))

final <- (variogr + covariog) / VY
final



out <- list()
for(j in 1:ncol(prove)) {
  out[[j]] <- cols_f(prove[, j])
}

da <- c(prove[1], rep(prove[2:(length(prove) - 1)], each = 2), tail(prove, 1))
matrix(da, nrow = length(da) / 2, byrow = TRUE)














di <- rep(prove[-c(prove[1], tail(prove, 1))], each = 2)



da[-c(1, tail(da, 1))]

tail(prove, 1)










pairs.points <- lapply(1:length(params), function(x) 
  lapply(1:ncol(out[[x]]), function(y) 
    t(combn(out[[x]][, y], 2))))
pairs.points <- lapply(pairs.points, function(x) do.call(rbind, x))





vars_f <- function(x)  1 / 2 * mean((x[, 1] - x[, 2]) ^ 2) + (1 / 2 * cov(x[, 1], x[, 2]))


variogr <- unlist(lapply(pairs.points, vars_f))
variogr / VY

# Compute variogram
variogr <- unlist(lapply(pairs.points, function(x) 1/ (2 * nrow(x)) * sum(x[, 1] - x[, 2]) ^ 2))
# Covariogram
covariog <- unlist(lapply(pairs.points, function(x) cov(x[, 1], x[, 2])))
# Total variance
VY <- variance_f(do.call(rbind, pairs.points))

final <- (variogr + covariog) / VY
final



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
VY

Vi <- 1 / (2 * N) * Rfast::colsums((Y_A - Y_AB) ^ 2) 
Vi
Vi/VY





























































# SETTINGS --------------------------------------------------------------------
N <- 50
params <- paste("X", 1:3, sep = "")
h <- 0.1

# Create STAR-VARS
mat <- vars_matrices(N = N, params = params, h = h)

# Check that number of rows of mat is identical to what is expected according to Razavi & Gupta
nrow(mat) == N * (length(params) * ((1 / h) -1) + 1)




# Compute model output
Y <- sensobol::ishigami_Fun(mat)

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
VY <- var(Y[index.centers])

# Compute variogram
variogr <- unlist(lapply(prove2, function(x) 1/ (2 * nrow(x)) * sum(x[, 1] - x[, 2]) ^ 2))





# Covariogram
covariog <- unlist(lapply(prove2, function(x) 1 / (2 * nrow(x)) * cov(x[, 1], x[, 2])))

# Final value
(variogr + covariog) / VY













