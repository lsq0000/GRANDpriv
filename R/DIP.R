#' @importFrom rmutil rlaplace
#' @importFrom EnvStats pemp qemp
#' @importFrom stats runif dist quantile

G <- function(x, b = 1) {
  ex <- exp(x / b)
  e1 <- exp(-1 / b)
  y <- 0.5 * b * ex * (1 - e1) * (x < 0) +
    (x + 0.5 * b / ex - 0.5 * b * ex * e1) * (x >= 0 & x <= 1) +
    (1 - 0.5 * b / e1 / ex * (1 - e1)) * (x > 1)

  return(y)
}

DIP.univariate <- function(Y, ep, distn, value) {
  n <- length(Y)
  eta <- rlaplace(n, 0, 1 / ep)

  if (distn == "empirical") {
    if (length(value) == 2) {
      if (value[[1]] == "continuous") {
        FY <- pemp(Y, obs = value[[2]])
        Z <- FY + eta
        GZ <- G(Z, 1 / ep)
        Yptb <- qemp(GZ, obs = value[[2]])
      } else if (value[[1]] == "discrete") {
        U <- runif(n)
        V <- Y - U
        FV <- pemp(V, obs = value[[2]])
        Z <- FV + eta
        GZ <- G(Z, 1 / ep)
        Yptb_raw <- qemp(GZ, obs = value[[2]])
        Yptb <- ceiling(Yptb_raw)
      } else {
        stop("value[[1]] for empirical distributions has to be either `continuous` or `discrete`.")
      }
    } else {
      stop("Parameter vector length has to be 2.")
    }
  } else {
    stop("The input distribution is not supported.")
  }

  return(Yptb)
}

DIP.conditional <- function(zp, zlp, zlo, ep, X, rho = 0.05) {
  x.full <- as.matrix(cbind(zlo, zp))
  n <- nrow(x.full)
  p <- ncol(x.full)

  zlp <- as.matrix(zlp)
  zlo <- as.matrix(zlo)
  xo <- X[, 1:p]

  A <- zlo
  B <- as.matrix(xo[, 1:(p - 1)])

  dist_matrix <- as.matrix(dist(rbind(A, B), method = "manhattan"))
  n_A <- nrow(A)
  n_B <- nrow(B)
  dst <- dist_matrix[1:n_A, (n_A + 1):(n_A + n_B)]
  ind <- t(apply(dst, 1, function(x) order(x)[1:ceiling(n_B * rho)]))

  xo_ind <- matrix(xo[c(ind), p], dim(ind)[1], dim(ind)[2])
  xo_ind_minus_zp <- (xo_ind - zp %*% t(rep(1, dim(ind)[2])) <= 0)
  Fx2i <- apply(xo_ind_minus_zp, 1, mean)

  eta <- rlaplace(n_A, 0, 1 / ep)
  z <- Fx2i + eta
  gz <- G(z, 1 / ep)

  C <- zlp

  dist_matrix2 <- as.matrix(dist(rbind(C, B), method = "manhattan"))
  n_C <- nrow(C)
  dst2 <- dist_matrix2[1:n_C, (n_C + 1):(n_C + n_B)]
  ind2 <- t(apply(dst2, 1, function(x) order(x)[1:ceiling(n_B * rho)]))

  xo_ind2 <- matrix(xo[c(ind2), p], dim(ind2)[1], dim(ind2)[2])
  xp.zip <- mapply(quantile, split(xo_ind2, row(xo_ind2)), gz)
  xt2 <- as.vector(xp.zip)

  return(xt2)
}

DIP.multivariate <- function(Z, ep, X, rho = 0.05) {
  p <- dim(Z)[2]
  s <- dim(Z)[1]
  n <- dim(X)[1]

  disc_ind <- sapply(1:p, function(j) length(unique(Z[, j])) < s)
  p1 <- sum(disc_ind)

  XC <- X
  XC[, disc_ind] <- X[, disc_ind] - matrix(runif(n * p1), n, p1)
  ZC <- Z
  ZC[, disc_ind] <- Z[, disc_ind] - matrix(runif(s * p1), s, p1)

  ZCdip <- ZC
  ZCdip[, 1] <- DIP.univariate(ZC[, 1], ep / p, "empirical", list("continuous", XC[, 1]))

  for (j in 2:p) {
    ZCdip[, j] <- DIP.conditional(ZC[, j], ZCdip[, 1:(j - 1)], ZC[, 1:(j - 1)], ep / p, XC, rho)
  }

  Zdip <- ZCdip
  Zdip[, disc_ind] <- ceiling(ZCdip[, disc_ind])

  for (j in which(disc_ind)) {
    Zdip[, j] <- pmax(Zdip[, j], min(X[, j]))
  }

  return(Zdip)
}
