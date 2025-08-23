#' @importFrom RSpectra eigs_sym
#' @importFrom diffpriv DPParamsEps DPMechLaplace releaseResponse
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats glm binomial lm
#' @importFrom randnet LSM.PGD
#' @importFrom igraph graph.adjacency degree count_triangles eigen_centrality harmonic_centrality
#' @importFrom HCD gen.A.from.P
#' @importFrom transport wasserstein1d

f <- function(X) X

ase <- function(A, K) {
  A2 <- t(A) %*% A
  eig <- eigs_sym(A2, k = K)
  Xhat <- t(t(eig$vectors[, 1:K, drop = FALSE]) * eig$values[1:K]^(1/4))

  return(Xhat)
}

Add.Laplace <- function(X, eps = 1) {
  p <- ncol(X)
  n <- nrow(X)
  pparams <- DPParamsEps(epsilon = eps / p)
  X_lap <- X

  for (j in 1:p) {
    mechanism <- DPMechLaplace(target = f, sensitivity = max(abs(X[, j])), dims = n)
    r_lap <- releaseResponse(mechanism, privacyParams = pparams, X = X[, j])
    X_lap[, j] <- r_lap$response
  }

  return(X_lap)
}

#' Generate Latent Space Model Network
#'
#' @title Generate Latent Space Model Network
#' @description Generates a random network following LSM (Latent Space Model) with specified parameters.
#'   LSM assumes that each node has a latent position in a k-dimensional space, and edge probabilities
#'   depend on the inner products between these latent positions.
#' @param n Integer. Number of nodes in the network.
#' @param k Integer. Dimension of the latent space.
#' @param K Integer. Number of communities/groups.
#' @param avg.d Numeric. Target average node degree. If NULL, no degree adjustment is performed. Default is "NULL".
#' @details 
#'   The Latent Space Model generates networks based on the following process:
#'   \enumerate{
#'     \item Node-specific intercept parameters: \eqn{\alpha_i \sim \mathsf{Unif}(1, 3)}, then transformed as \eqn{\alpha_i = -\alpha_i/2}
#'     \item Community assignments: Each node is randomly assigned to one of K communities with equal probability
#'     \item Community centers: \eqn{\mu_g} sampled uniformly from \eqn{[-1, 1]^k} hypercube
#'     \item Latent positions: \eqn{Z_i \sim \mathcal{N}_{[-2, 2]}(\mu_{g_i}, I_k)} where \eqn{\mu_{g_i}} is the community center for node i's community
#'     \item Edge probabilities: \deqn{P_{ij} = \text{logit}^{-1}(\alpha_i + \alpha_j + Z_i^\top Z_j) = \frac{1}{1 + \exp(-(\alpha_i + \alpha_j + Z_i^\top Z_j))}}
#'     \item Adjacency matrix: \eqn{A_{ij} \sim \mathsf{Bern}(P_{ij})} for \eqn{i < j}, with \eqn{A_{ii} = 0} and \eqn{A_{ji} = A_{ij}}
#'   }
#'   If \code{avg.d} is specified, the edge probabilities are scaled to achieve the target average node degree.
#'   
#'   For more details, see the "simulation studies" section of Ma, Ma, and Yuan (2020), noting that our implementation has slight differences.
#' @return A list containing:
#' \itemize{
#'   \item A: Adjacency matrix of the generated network
#'   \item g: Graph object of the generated network
#'   \item P: Probability matrix of the generated network
#'   \item alpha: Node-specific intercept parameters
#'   \item Z: Latent positions in k-dimensional space
#'   \item idx: Community assignments for each node
#' }
#' @references 
#'   P. D. Hoff, A. E. Raftery, and M. S. Handcock. Latent space approaches to social network analysis. 
#'   Journal of the American Statistical Association, 97(460):1090–1098, 2002.
#'   
#'   Z. Ma, Z. Ma, and H. Yuan. Universal latent space model fitting for large networks with edge covariates.
#'   Journal of Machine Learning Research, 21(4):1–67, 2020.
#'   
#'   S. Liu, X. Bi, and T. Li. GRAND: Graph Release with Assured Node Differential Privacy. 
#'   arXiv preprint arXiv:2507.00402, 2025.
#' @export
#' @examples
#' # Generate a network with 400 nodes, 2D latent space, 3 communities
#' network <- LSM.Gen(n = 400, k = 2, K = 3)
#' # Generate with target average degree of 40
#' network2 <- LSM.Gen(n = 400, k = 2, K = 3, avg.d = 40)
LSM.Gen <- function(n, k, K, avg.d = NULL) {
  alpha <- runif(n, 1, 3)
  alpha <- -alpha / 2

  mu <- matrix(runif(K * k, -1, 1), K, k)
  idx <- sample(1:K, n, replace = TRUE)
  mu <- mu[idx, ]

  Z <- mu + matrix(rtruncnorm(n * k, -2, 2), n, k)
  J <- diag(n) - rep(1, n) %*% t(rep(1, n)) / n
  Z <- J %*% Z
  G <- Z %*% t(Z)
  Z <- Z / sqrt(sqrt(sum(G^2)) / n)

  theta <- alpha %*% t(rep(1, n)) + rep(1, n) %*% t(alpha) + Z %*% t(Z)
  P <- 1 / (1 + exp(-theta))

  if (!is.null(avg.d)) {
    for (i in 1:10) {
      ratio <- mean(rowSums(P)) / avg.d
      beta <- -log(ratio)
      alpha <- alpha + beta / 2
      theta <- alpha %*% t(rep(1, n)) + rep(1, n) %*% t(alpha) + Z %*% t(Z)
      P <- 1 / (1 + exp(-theta))
    }
  }

  upper.index <- which(upper.tri(P))
  upper.p <- P[upper.index]
  upper.u <- runif(length(upper.p))
  upper.A <- rep(0, length(upper.p))
  upper.A[upper.u < upper.p] <- 1
  A <- matrix(0, n, n)
  A[upper.index] <- upper.A
  A <- A + t(A)
  diag(A) <- 0
  
  g <- graph.adjacency(A, "undirected")

  return(list(A = A, g = g, P = P, alpha = alpha, Z = Z, idx = idx))
}

GRAND.one.node <- function(A, given.index, new.index, given.Z, given.alpha = NULL, model = c("LSM", "RDPG")) {
  model <- match.arg(model)

  Y <- A[new.index, given.index]

  if (model == "LSM") {
    if (is.null(given.alpha)) {
      stop("given.alpha must be provided for LSM model.")
    }

    fit <- glm(Y ~ given.Z, family = binomial(), offset = given.alpha, maxit = 1000)
  } else {
    fit <- lm(Y ~ given.Z - 1)
  }

  return(fit)
}

GRAND.estimate <- function(A, K, holdout.index, release.index, model = c("LSM", "RDPG"), niter = 500, verbose = TRUE) {
  model <- match.arg(model)

  m <- length(holdout.index)
  n <- length(release.index)
  A22 <- A[holdout.index, holdout.index]

  if (model == "LSM") {
    fitting.count <- 0

    while (fitting.count < 5) {
      fit.holdout <- LSM.PGD(A22, k = K, niter = niter)
      fitting.count <- fitting.count + 1

      if (sum(is.na(fit.holdout$Z)) == 0) {
        break
      }
    }

    if (sum(colSums(abs(fit.holdout$Z)) < 1e-4) > 0) {
      stop(paste0("Potential deficiency of ", sum(colSums(abs(fit.holdout$Z)) < 1e-4), "."))
    }

    if (verbose) cat("PGD used", length(fit.holdout$obj), "iterations.\n")

    Z.holdout <- fit.holdout$Z
    alpha.holdout <- fit.holdout$alpha
  } else {
    Z.holdout <- ase(A22, K)
    alpha.holdout <- NULL
  }

  Z.hat <- matrix(0, nrow = m + n, ncol = K)
  alpha.hat <- rep(0, m + n)
  Z.hat[holdout.index, ] <- Z.holdout
  if (!is.null(alpha.holdout)) {
    alpha.hat[holdout.index] <- alpha.holdout
  }

  for (j in 1:n) {
    fit <- GRAND.one.node(A = A,
                          given.index = holdout.index,
                          new.index = release.index[j],
                          given.Z = Z.holdout,
                          given.alpha = alpha.holdout,
                          model = model)

    if (model == "LSM") {
      Z.hat[release.index[j], ] <- fit$coefficients[2:(K + 1)]
      alpha.hat[release.index[j]] <- fit$coefficients[1]
    } else {
      Z.hat[release.index[j], ] <- fit$coefficients[1:K]
    }
  }

  if (model == "LSM") {
    return(list(Z = Z.hat, alpha = alpha.hat))
  } else {
    return(list(Z = Z.hat))
  }
}

#' GRAND Privatization of Network Data
#'
#' @title GRAND Privatization of Network Data
#' @description Applies the GRAND (Graph Release with Assured Node Differential privacy) method
#' to privatize network data using differential privacy. The method estimates latent positions
#' from the network and applies multivariate differential privacy to protect sensitive information.
#' @param A Matrix. Adjacency matrix of the input network.
#' @param K Integer. Dimension of the latent space for network embedding.
#' @param idx Integer vector. Indices of nodes to be privatized.
#' @param eps Numeric or numeric vector. Privacy budget parameter(s) for differential privacy. Default is 1.
#' @param model Character. Model type, either "LSM" (Latent Space Model) or "RDPG" (Random Dot Product Graph). Default is "LSM".
#' @param niter Integer. Number of iterations for the optimization algorithm. Default is 500.
#' @param rho Numeric. Parameter controlling the neighborhood size for conditional distributions. Default is 0.05.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#' @details 
#'   The GRAND privatization algorithm consists of the following steps:
#'   \enumerate{
#'     \item \strong{Network partitioning}: Split nodes into release set (idx) and holdout set
#'     \item \strong{Latent position estimation}: Use LSM.PGD (Projected Gradient Descent) or ASE (Adjacency Spectral Embedding) to estimate latent positions
#'     \item \strong{Differential privacy}: Apply multivariate DIP (Distribution-Invariant differential Privacy) mechanism to protect latent positions
#'     \item \strong{Network reconstruction}: Generate privatized networks from perturbed latent positions
#'   }
#'   The method also computes standard Laplace mechanism results for comparison.
#'
#'   For more details, see the "proposed method" section of Liu, Bi, and Li (2025).
#' @return A list containing:
#' \itemize{
#'   \item non.private.result: Results without privacy, including the original and estimated data
#'   \item GRAND.result: List with one element per epsilon value. Each element contains GRAND privatization results for that specific epsilon
#'   \item Laplace.result: List with one element per epsilon value. Each element contains baseline Laplace mechanism results for that specific epsilon
#'   \item eps: Vector of privacy budget parameter(s) used
#' }
#' @references 
#'   P. D. Hoff, A. E. Raftery, and M. S. Handcock. Latent space approaches to social network analysis. 
#'   Journal of the American Statistical Association, 97(460):1090–1098, 2002.
#'   
#'   S. J. Young and E. R. Scheinerman. Random dot product graph models for social networks. 
#'   In International Workshop on Algorithms and Models for the Web-Graph, pages 138–149. Springer, 2007.
#'   
#'   Z. Ma, Z. Ma, and H. Yuan. Universal latent space model fitting for large networks with edge covariates. 
#'   Journal of Machine Learning Research, 21(4):1–67, 2020.
#'
#'   A. Athreya, D. E. Fishkind, M. Tang, C. E. Priebe, Y. Park, J. T. Vogelstein, K. Levin, V. Lyzinski, Y. Qin, and D. L. Sussman. Statistical inference on random dot product graphs: a survey. 
#'   Journal of Machine Learning Research, 18(226):1–92, 2018.
#'
#'   P. Rubin-Delanchy, J. Cape, M. Tang, and C. E. Priebe. A statistical interpretation of spectral embedding: The generalised random dot product graph. 
#'   Journal of the Royal Statistical Society Series B: Statistical Methodology, 84(4):1446–1473, 2022.
#'
#'   X. Bi and X. Shen. Distribution-invariant differential privacy. 
#'   Journal of Econometrics, 235(2):444–453, 2023.
#'
#'   S. Liu, X. Bi, and T. Li. GRAND: Graph Release with Assured Node Differential Privacy. 
#'   arXiv preprint arXiv:2507.00402, 2025.
#' @export
#' @examples
#' # Generate a sample network
#' network <- LSM.Gen(n = 400, k = 2, K = 3, avg.d = 40)
#' # Privatize the first 200 nodes with epsilon = 1, 2, 5, 10
#' result <- GRAND.privatize(A = network$A, K = 2, idx = 1:200, eps = c(1, 2, 5, 10), model = "LSM")
GRAND.privatize <- function(A, K, idx, eps = 1, model = c("LSM", "RDPG"), niter = 500, rho = 0.05, verbose = TRUE) {
  model <- match.arg(model)

  n <- nrow(A)
  n1 <- length(idx)
  n2 <- n - n1
  holdout.idx <- setdiff(1:n, idx)
  A11 <- A[idx, idx]
  g1 <- graph.adjacency(A11, "undirected")
  A22 <- A[holdout.idx, holdout.idx]
  g2 <- graph.adjacency(A22, "undirected")

  fit <- GRAND.estimate(A = A,
                        K = K,
                        holdout.index = holdout.idx,
                        release.index = idx,
                        model = model,
                        niter = niter,
                        verbose = verbose)
  X.trans <- if (model == "LSM") cbind(fit$alpha, fit$Z) else fit$Z

  X1.hat <- X.trans[idx, , drop = FALSE]
  X2.hat <- X.trans[holdout.idx, , drop = FALSE]

  if (model == "LSM") {
    alpha1.hat <- X1.hat[, 1, drop = FALSE]
    Z1.hat <- X1.hat[, -1, drop = FALSE]
    theta1.hat <- alpha1.hat %*% t(rep(1, n1)) + rep(1, n1) %*% t(alpha1.hat) + Z1.hat %*% t(Z1.hat)
    P1.hat <- 1 / (1 + exp(-theta1.hat))
  } else {
    Z1.hat <- X1.hat
    theta1.hat <- Z1.hat %*% t(Z1.hat)
    P1.hat <- pmin(pmax(theta1.hat, 1e-5), 1 - 1e-5)
  }
  A1.hat <- gen.A.from.P(P1.hat)
  g1.hat <- graph.adjacency(A1.hat, "undirected")

  if (model == "LSM") {
    alpha2.hat <- X2.hat[, 1, drop = FALSE]
    Z2.hat <- X2.hat[, -1, drop = FALSE]
    theta2.hat <- alpha2.hat %*% t(rep(1, n2)) + rep(1, n2) %*% t(alpha2.hat) + Z2.hat %*% t(Z2.hat)
    P2.hat <- 1 / (1 + exp(-theta2.hat))
  } else {
    Z2.hat <- X2.hat
    theta2.hat <- Z2.hat %*% t(Z2.hat)
    P2.hat <- pmin(pmax(theta2.hat, 1e-5), 1 - 1e-5)
  }
  A2.hat <- gen.A.from.P(P2.hat)
  g2.hat <- graph.adjacency(A2.hat, "undirected")

  non.private.result <- list(A1 = A11, A2 = A22,
                             A1.hat = A1.hat, A2.hat = A2.hat,
                             P1.hat = P1.hat, P2.hat = P2.hat,
                             g1 = g1, g2 = g2,
                             g1.hat = g1.hat, g2.hat = g2.hat,
                             X1.hat = X1.hat, X2.hat = X2.hat)
  L <- length(eps)
  GRAND.result <- list()

  for (jj in 1:L) {
    if (verbose) cat(paste0("Calling GRAND with \u03B5=", eps[jj], ".\n"))
    X1.dip <- DIP.multivariate(X1.hat, eps[jj], X2.hat, rho)

    if (model == "LSM") {
      alpha1.dip <- X1.dip[, 1, drop = FALSE]
      Z1.dip <- X1.dip[, -1, drop = FALSE]
      theta1.dip <- alpha1.dip %*% t(rep(1, n1)) + rep(1, n1) %*% t(alpha1.dip) + Z1.dip %*% t(Z1.dip)
      P1.dip <- 1 / (1 + exp(-theta1.dip))
    } else {
      Z1.dip <- X1.dip
      theta1.dip <- Z1.dip %*% t(Z1.dip)
      P1.dip <- pmin(pmax(theta1.dip, 1e-5), 1 - 1e-5)
    }
    A1.dip <- gen.A.from.P(P1.dip)
    g1.dip <- graph.adjacency(A1.dip, "undirected")
    GRAND.result[[jj]] <- list(A1.grand = A1.dip,
                               P1.grand = P1.dip,
                               g1.grand = g1.dip,
                               X1.grand = X1.dip)
  }
  if (verbose) cat("Finish GRAND.\n")

  Laplace.result <- list()
  for (jj in 1:L) {
    if (verbose) cat(paste0("Calling Laplace with \u03B5=", eps[jj], ".\n"))
    X1.Lap <- Add.Laplace(X = X1.hat, eps = eps[jj])

    if (model == "LSM") {
      alpha1.Lap <- X1.Lap[, 1, drop = FALSE]
      Z1.Lap <- X1.Lap[, -1, drop = FALSE]
      theta1.Lap <- alpha1.Lap %*% t(rep(1, n1)) + rep(1, n1) %*% t(alpha1.Lap) + Z1.Lap %*% t(Z1.Lap)
      P1.Lap <- 1 / (1 + exp(-theta1.Lap))
    } else {
      Z1.Lap <- X1.Lap
      theta1.Lap <- Z1.Lap %*% t(Z1.Lap)
      P1.Lap <- pmin(pmax(theta1.Lap, 1e-5), 1 - 1e-5)
    }
    A1.Lap <- gen.A.from.P(P1.Lap)
    g1.Lap <- graph.adjacency(A1.Lap, "undirected")
    Laplace.result[[jj]] <- list(A1.Lap = A1.Lap,
                                 P1.Lap = P1.Lap,
                                 g1.Lap = g1.Lap,
                                 X1.Lap = X1.Lap)
  }
  if (verbose) cat("Finish Laplace.\n")

  return(list(non.private.result = non.private.result, GRAND.result = GRAND.result, Laplace.result = Laplace.result, eps = eps))
}

GRAND.evaluate.degree <- function(result) {
  degree.true <- log(1 + degree(result$non.private.result$g1))
  degree.hat <- log(1 + degree(result$non.private.result$g1.hat))
  degree.true2 <- log(1 + degree(result$non.private.result$g2))
  degree.hat2 <- log(1 + degree(result$non.private.result$g2.hat))
  degree.grand <- lapply(result$GRAND.result, function(x) log(1 + degree(x$g1.grand)))
  degree.lap <- lapply(result$Laplace.result, function(x) log(1 + degree(x$g1.Lap)))
  degree.mat <- data.frame(stat = rep("Node Degree", length(result$eps)),
                           eps = result$eps,
                           Hat = rep(wasserstein1d(degree.true, degree.hat), length(result$eps)),
                           Hat2 = rep(wasserstein1d(degree.true2, degree.hat2), length(result$eps)),
                           GRAND = unlist(lapply(degree.grand, function(x) wasserstein1d(degree.true, x))),
                           Laplace = unlist(lapply(degree.lap, function(x) wasserstein1d(degree.true, x))))

  return(degree.mat)
}

GRAND.evaluate.vshape <- function(result) {
  vs.true <- log(1 + get.v(result$non.private.result$g1))
  vs.hat <- log(1 + get.v(result$non.private.result$g1.hat))
  vs.true2 <- log(1 + get.v(result$non.private.result$g2))
  vs.hat2 <- log(1 + get.v(result$non.private.result$g2.hat))
  vs.grand <- lapply(result$GRAND.result, function(x) log(1 + get.v(x$g1.grand)))
  vs.lap <- lapply(result$Laplace.result, function(x) log(1 + get.v(x$g1.Lap)))
  vs.mat <- data.frame(stat = rep("V-Shape Count", length(result$eps)),
                       eps = result$eps,
                       Hat = rep(wasserstein1d(vs.true, vs.hat), length(result$eps)),
                       Hat2 = rep(wasserstein1d(vs.true2, vs.hat2), length(result$eps)),
                       GRAND = unlist(lapply(vs.grand, function(x) wasserstein1d(vs.true, x))),
                       Laplace = unlist(lapply(vs.lap, function(x) wasserstein1d(vs.true, x))))


  return(vs.mat)
}

GRAND.evaluate.triangle <- function(result) {
  tri.true <- log(1 + count_triangles(result$non.private.result$g1))
  tri.hat <- log(1 + count_triangles(result$non.private.result$g1.hat))
  tri.true2 <- log(1 + count_triangles(result$non.private.result$g2))
  tri.hat2 <- log(1 + count_triangles(result$non.private.result$g2.hat))
  tri.grand <- lapply(result$GRAND.result, function(x) log(1 + count_triangles(x$g1.grand)))
  tri.lap <- lapply(result$Laplace.result, function(x) log(1 + count_triangles(x$g1.Lap)))
  tri.mat <- data.frame(stat = rep("Triangle Count", length(result$eps)),
                        eps = result$eps,
                        Hat = rep(wasserstein1d(tri.true, tri.hat), length(result$eps)),
                        Hat2 = rep(wasserstein1d(tri.true2, tri.hat2), length(result$eps)),
                        GRAND = unlist(lapply(tri.grand, function(x) wasserstein1d(tri.true, x))),
                        Laplace = unlist(lapply(tri.lap, function(x) wasserstein1d(tri.true, x))))

  return(tri.mat)
}

GRAND.evaluate.eigen <- function(result) {
  eigen.true <- eigen_centrality(result$non.private.result$g1)$vector
  eigen.hat <- eigen_centrality(result$non.private.result$g1.hat)$vector
  eigen.true2 <- eigen_centrality(result$non.private.result$g2)$vector
  eigen.hat2 <- eigen_centrality(result$non.private.result$g2.hat)$vector
  eigen.grand <- lapply(result$GRAND.result, function(x) eigen_centrality(x$g1.grand)$vector)
  eigen.lap <- lapply(result$Laplace.result, function(x) eigen_centrality(x$g1.Lap)$vector)
  eigen.mat <- data.frame(stat = rep("Eigen Centrality", length(result$eps)),
                          eps = result$eps,
                          Hat = rep(wasserstein1d(eigen.true, eigen.hat), length(result$eps)),
                          Hat2 = rep(wasserstein1d(eigen.true2, eigen.hat2), length(result$eps)),
                          GRAND = unlist(lapply(eigen.grand, function(x) wasserstein1d(eigen.true, x))),
                          Laplace = unlist(lapply(eigen.lap, function(x) wasserstein1d(eigen.true, x))))


  return(eigen.mat)
}

GRAND.evaluate.harmonic <- function(result) {
  harmonic.true <- harmonic_centrality(result$non.private.result$g1)
  harmonic.hat <- harmonic_centrality(result$non.private.result$g1.hat)
  harmonic.true2 <- harmonic_centrality(result$non.private.result$g2)
  harmonic.hat2 <- harmonic_centrality(result$non.private.result$g2.hat)
  harmonic.grand <- lapply(result$GRAND.result, function(x) harmonic_centrality(x$g1.grand))
  harmonic.lap <- lapply(result$Laplace.result, function(x) harmonic_centrality(x$g1.Lap))
  harmonic.mat <- data.frame(stat = rep("Harmonic Centrality", length(result$eps)),
                             eps = result$eps,
                             Hat = rep(wasserstein1d(harmonic.true, harmonic.hat), length(result$eps)),
                             Hat2 = rep(wasserstein1d(harmonic.true2, harmonic.hat2), length(result$eps)),
                             GRAND = unlist(lapply(harmonic.grand, function(x) wasserstein1d(harmonic.true, x))),
                             Laplace = unlist(lapply(harmonic.lap, function(x) wasserstein1d(harmonic.true, x))))


  return(harmonic.mat)
}

get.v <- function(g) {
  deg_vec <- degree(g)
  twostar_counts <- choose(deg_vec, 2)

  return(twostar_counts)
}

#' GRAND Evaluation of Network Data
#'
#' @title GRAND Evaluation of Network Data
#' @description Evaluates the quality of GRAND privatization results by comparing
#' various network statistics between the original and privatized networks using
#' Wasserstein distance.
#' @param result List. Output from GRAND.privatize function containing privatization results.
#' @param statistics Character vector. Network statistics to evaluate. Options include:
#' "degree" (Node Degree), "vshape" (V-Shape Count), "triangle" (Triangle Count), "eigen" (Eigen Centrality), "harmonic" (Harmonic Centrality). Default is all statistics.
#' @details 
#'   This function evaluates privatization quality by comparing network statistics between 
#'   original and privatized networks using Wasserstein-1 distance. The evaluation covers:
#'   \itemize{
#'     \item \strong{Hat}: Performance on release set (without privacy)
#'     \item \strong{Hat2}: Performance on holdout set (without privacy)
#'     \item \strong{GRAND}: Performance of GRAND privatization method
#'     \item \strong{Laplace}: Performance of standard Laplace mechanism
#'   }
#'   Lower Wasserstein distances indicate better utility preservation.
#'
#'   For more details, see the "simulation experiments" section of Liu, Bi, and Li (2025).
#' @return A data frame containing evaluation results with columns:
#' \itemize{
#'   \item stat: Type of network statistic(s) evaluated
#'   \item eps: Privacy budget parameter(s) used
#'   \item Hat: Wasserstein distance for release set estimation
#'   \item Hat2: Wasserstein distance for holdout set estimation
#'   \item GRAND: Wasserstein distance for GRAND privatization method
#'   \item Laplace: Wasserstein distance for standard Laplace mechanism
#' }
#' @references 
#'   S. Liu, X. Bi, and T. Li. GRAND: Graph Release with Assured Node Differential Privacy. 
#'   arXiv preprint arXiv:2507.00402, 2025.
#' @export
#' @examples
#' # Generate and privatize a network
#' network <- LSM.Gen(n = 400, k = 2, K = 3, avg.d = 40)
#' result <- GRAND.privatize(A = network$A, K = 2, idx = 1:200, eps = c(1, 2, 5, 10), model = "LSM")
#' # Evaluate results for all statistics
#' evaluation <- GRAND.evaluate(result)
#' # Evaluate only degree and triangle statistics
#' evaluation_subset <- GRAND.evaluate(result, statistics = c("degree", "triangle"))
GRAND.evaluate <- function(result, statistics = c("degree", "vshape", "triangle", "eigen", "harmonic")) {
  statistic_funcs <- list(degree = GRAND.evaluate.degree,
                          vshape = GRAND.evaluate.vshape,
                          triangle = GRAND.evaluate.triangle,
                          eigen = GRAND.evaluate.eigen,
                          harmonic = GRAND.evaluate.harmonic)

  selected_funcs <- statistic_funcs[statistics]
  results <- lapply(selected_funcs, function(f) f(result))
  output <- do.call(rbind, results)
  rownames(output) <- NULL

  return(output)
}
