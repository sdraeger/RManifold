library(class)
library(tidyverse)
library(rgl)
library(RDRToolbox)
library(cccd)
library(Rfast)
library(matlib)
library(RANN)
library(SNFtool)
library(igraph)
library(PRIMME)
library(pracma)
library(Rtsne)


isomap <- function(X, k, d) {
  # Compute k neighbor graph
  graph = nng(data, k=k)
  # Compute distance matrix
  ds = dist(data)
  dsm = as.matrix(ds)
  # Compute shortest pairwise paths between vertices
  shortest_paths = floyd(dsm)
  G = shortest_paths ^ 2
  G = G * -0.5
  # Apply MDS on distance matrix
  fit = cmdscale(G, eig=TRUE, k=d)
  fit
}

pairwiseDists = function(data){
  # num_samples = nrow(data)
  #
  # d = t(rowSums(data^2))
  # d = d[rep(1, each = num_samples),]
  # d = d + t(d) - 2*data%*%t(data)
  # diag(d) = 0
  # d = sqrt(d)

  #return(d)
  dist(data)
}

lle <- function(X, dim, k) {
  k = round(k)
  dim = round(dim)

  num_samples = nrow(data)
  num_features = ncol(data)

  ## compute pairwise distances
  d = pairwiseDists(data)

  ## determine the k nearest neighbors
  sort_idx = apply(d,2,order)
  neighbors = sort_idx[2:(k+1),]

  W = matrix(0,k,num_samples)
  for(i in 1:num_samples){
    ## compute covariance matrix of the neighbors of the ith sample (center and square samples)
    N = t(X[neighbors[,i],]) - X[i,]
    Cov = t(N) %*% N

    ## weights W are solution of C * W = 1
    Cov = Cov + reg * diag(dim(Cov)[1]) * sum(diag(Cov))
    W[,i] = solve(Cov, rep(1, k))
    ## rescale weights that rows sum to one
    W[,i] = W[,i] / sum(W[,i])
  }

  ## build the cost matrix M = (I-W)^T (I-W)
  M = diag(1,num_samples)
  for(i in 1:num_samples){
    w = W[,i]
    n = neighbors[,i]
    M[i,n] = M[i,n] - t(w)
    M[n,i] = M[n,i] - w
    M[n,n] = M[n,n] + w %*% t(w)
  }

  ## low dimensional embedding Y given by the eigenvectors belonging to the (dim+1) smallest eigenvalues (first eigenvalue is trivial)
  eig_M = eigen(M)
  sweep(eig_M$vectors, 2, sqrt(colSums(eig_M$vectors^2)), "/")    ## normalize eingenvectors
  Y = eig_M$vectors[,(num_samples-1):(num_samples-dim)] * sqrt(num_samples)

  Y
}

mds <- function(X, d) {
  N = dim(X)[1]
  D = as.matrix(pairwiseDists(X))
  m = rep(1, N) / N
  C = diag(N) - matrix(1, N, N)
  S = -0.5 * (C %*% (D %*% t(C)))
  S_eig = eigen(S)
  S_eigvals = S_eig$values
  S_eigvecs = S_eig$vectors
  S_eigvals[S_eigvals < 0] = 0
  Lambda = diag(S_eigvals)
  U = S_eigvecs
  scores = sqrt(inv(diag(m))) %*% (U %*% sqrt(Lambda))
  scores
}

hessian_lle <- function(X, d, k) {
  indata <- data
  n <- nrow(indata)
  hs <- d * (d + 1) / 2
  W <- Matrix::sparseMatrix(i = numeric(0),
                            j = numeric(0),
                            x = numeric(0),
                            dims = c(n, hs * n))
  ii <- jj <- ww <- list()
  ## Identify neighbors
  nnidx <- RANN::nn2(data = indata, query = indata, k = k + 1,
                     treetype = "kd", "standard", eps = 0)$nn.idx
  for (i in seq_len(n)) {
    cat(i, "/", n, "\r", sep = "")
    ## get neighborhood
    Nui <- indata[nnidx[i, ], , drop = FALSE]

    ## Form tangent coordinates:
    Nui <- sweep(Nui, 2, colMeans(Nui), "-")
    tc <- svd(Nui, nu = d, nv = 0)$u

    ## Develop Hessian Estimator
    Xi <- cbind(
      1, tc, tc ^ 2,
      apply(combn(seq_len(d), 2), 2,
            function(x) tc[, x[1]] * tc[, x[2]])
    )
    tHi <- qr.Q(qr(Xi))[, -(1:(d + 1)),
                        drop = FALSE]

    ## Add quadratic form to hessian
    ii[[i]] <- rep(nnidx[i, ], hs)
    jj[[i]] <- rep((i - 1) * hs + (1:hs), each = ncol(nnidx))
    ww[[i]] <- as.vector(tHi)
  }
  H <- as(Matrix::tcrossprod(Matrix::spMatrix(
    i = unlist(ii, FALSE, FALSE),
    j = unlist(jj, FALSE, FALSE),
    x = unlist(ww, FALSE, FALSE),
    nrow = n, ncol = n * hs)
  ), "dgCMatrix")

  ## Find null space
  ## eigs and eigs_sym converges much more reliably and faster
  ## with sigma = -eps than with which = "L*"
  outdata <- RSpectra::eigs_sym(H, k = d + 1, sigma = -1e-5)
  outdata <- outdata$vectors[, order(outdata$values)[-1], drop = FALSE]
  colnames(outdata) <- paste0("HLLE", seq_len(ncol(outdata)))

  outdata
}

spectral_embedding <- function(X, d) {
  adj_mat = as.matrix(pairwiseDists(X))
  G = graph_from_adjacency_matrix(adj_mat)
  Y = embed_adjacency_matrix(G, d)$Y
  Y
}

ltsa <- function(X, d, k, n_coord = 2) {
  N = dim(X)[1]
  nnidx <- RANN::nn2(data = data, query = data, k = k + 1,
                     treetype = "kd", "standard", eps = 0)$nn.idx
  M = matrix(0, N, N)

  for (i in 1:N) {
    nbrs_len = length(nnidx[i,]) - 1
    neighs = nnidx[i,2:(nbrs_len+1)]
    Xi = X[neighs,]
    Xi = Xi - rowMeans(Xi)

    if (k > dim(X)[2]) {
      v = svd(Xi, nrow(Xi), ncol(Xi))$u
    } else {
      Ci = Xi %*% t(Xi)
      v = eigen(Ci)$vectors
    }

    Gi = matrix(0, k, n_coord + 1)
    Gi[,2:(dim(Gi)[2])] = v[,1:n_coord]
    Gi[,1] = 1.0 / sqrt(k)

    GiGiT = Gi %*% t(Gi)
    M[neighs, neighs] = M[neighs, neighs] - GiGiT + 1
  }

  outdata <- RSpectra::eigs_sym(M, k = d + 1, sigma = -1e-5)
  outdata <- outdata$vectors[, order(outdata$values)[-1], drop = FALSE]
  outdata
}

tsne <- function(...) {
  Rtsne::Rtsne(...)
}

mlle <- function(X, d, k, n_coord = 2) {
  N = dim(X)[1]
  d_in = dim(X)[2]
  nnidx <- RANN::nn2(data = X, query = data, k = k + 1,
                     treetype = "kd", "standard", eps = 0)$nn.idx
  nnidx = nnidx[,2:(k+1)]
  V = array(0, dim=c(N, k, k))
  nev = min(d_in, k)
  evals = matrix(0, N, nev)

  if (k > d_in) {
    for (i in 1:N) {
      Nui <- X[nnidx[i,], , drop = FALSE]
      Nui <- sweep(Nui, 2, colMeans(Nui), "-")
      tc <- svd(Nui, k, k)
      V[i,,] = tc$u
      evals[i,] = tc$d
    }
    evals = evals ^ 2
  } else {
    for (i in 1:N) {
      Nui <- X[nnidx[i,], , drop = FALSE]
      Nui <- sweep(Nui, 2, colMeans(Nui), "-")
      C_nbrs = Nui %*% t(Nui)
      edecomp = eigen(C_nbrs)
      evecs = edecomp$vectors
      evals[i,] = edecomp$values
      V[i,,] = evecs
    }
  }

  W = matrix(0, k, N)
  for(i in 1:N){
    ## compute covariance matrix of the neighbors of the ith sample (center and square samples)
    neigh = t(X[nnidx[i,],]) - X[i,]
    Cov = t(neigh) %*% neigh

    ## weights W are solution of C * W = 1
    Cov = Cov + reg * diag(dim(Cov)[1]) * sum(diag(Cov))
    W[,i] = solve(Cov, rep(1, k))
    ## rescale weights that rows sum to one
    W[,i] = W[,i] / sum(W[,i])
  }

  rho = rowSums(evals[, n_coord:dim(evals)[2]])
  eta = median(rho)

  s_range = rep(0, N)
  evals_cumsum = t(apply(evals, 1, cumsum))
  eta_range = matrix(0, dim(evals_cumsum)[1], dim(evals_cumsum)[2] - 1)
  for (i in 1:dim(eta_range)[1]) {
    for (j in 1:(dim(eta_range)[2] - 1)) {
      eta_range[i, j] = evals_cumsum[i, dim(evals_cumsum)[2]] / evals_cumsum[i, j] - 1
    }
  }

  for (i in 1:N) {
    point <- which(order(c(eta, eta_range[i,])) == 1)
    s_range = append(s_range, eta, point - 1)
  }

  s_range = s_range + k - nev
  outdata <- RSpectra::eigs_sym(s_range, k = d + 1, sigma = -1e-5)
  outdata <- outdata$vectors[, order(outdata$values)[-1], drop = FALSE]
}

mlle(data, 2, 5)

tsne(data)

ltemb = ltsa(data, 2, 5)
df = data.frame(ltemb)
colnames(df) = c('x', 'y')

ggplot(df, aes(x=x, y=y)) +
  geom_point()

Semb = spectral_embedding(data, 2)

m = mds(data, 2)
l = lle(data, 2, 5)

df = data.frame(l)
colnames(df) = c('x', 'y')

df = data.frame(Semb)
colnames(df) = c('x', 'y')

ggplot(df, aes(x=x, y=y)) +
  geom_point()

# k = 5
# d = 2
# fit = isomap(data, k, d)
# emb_df = data.frame(fit$points)
# colnames(emb_df) = c('x', 'y')
# ggplot(emb_df, aes(x=x, y=y)) + geom_point()
