
N <- 10
params <- paste("X", 1:3, sep = "")
mat <- randtoolbox::sobol(N, length(params))

h <- 0.1
center <- mat[1, ]

seq(0.5 %% 0.1, 1, 0.1)

sections <- sapply(center, function(x) { 
  all <- seq(x %% h, 1, h)
  non.zeros <- all[all!=0]
})

B <- sapply(1:ncol(mat), function(x) sections[, x][!sections[, x] %in% center[x]])

A <- matrix(center, nrow = nrow(B), ncol = length(center), byrow = TRUE)


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
  return(out)
}

N <- 10
params <- paste("X", 1:3, sep ="")
h <- 0.1

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


