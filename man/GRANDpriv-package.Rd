\name{GRANDpriv-package}
\alias{GRANDpriv-package}
\alias{GRANDpriv}
\alias{_PACKAGE}
\docType{package}
\title{
\packageTitle{GRANDpriv}
}
\description{
\strong{GRANDpriv} (\strong{G}raph \strong{R}elease with \strong{A}ssured \strong{N}ode \strong{D}ifferential \strong{priv}acy) is an R package that implements a novel method for privatizing network data using differential privacy. The package provides functions for generating synthetic networks based on LSM (Latent Space Model), applying differential privacy to network latent positions to achieve overall network privatization, and evaluating the utility of privatized networks through various network statistics. The privatize and evaluate functions support both LSM and RDPG (Random Dot Product Graph). For generating RDPG networks, users are encouraged to use the \strong{randnet} package.
}
\details{
The package implements a two-step approach:
\enumerate{
  \item \strong{Latent Position Estimation}: Estimates latent positions from the network structure using either LSM.PGD (Projected Gradient Descent, for LSM) or ASE (Adjacency Spectral Embedding, for RDPG)
  \item \strong{Multivariate Differential Privacy}: Applies DIP (Distribution-Invariant differential Privacy) mechanism to protect latent positions while preserving network utility
}

\strong{Main Functions:}
\itemize{
  \item \code{\link{LSM.Gen}}: Generate Latent Space Model Network
  \item \code{\link{GRAND.privatize}}: GRAND Privatization of Network Data
  \item \code{\link{GRAND.evaluate}}: GRAND Evaluation of Network Data
}

\strong{Key Features:}
\itemize{
  \item Network Generation: Generate synthetic networks using Latent Space Models
  \item Differential Privacy: Apply node-level differential privacy to network data
  \item Multiple Models: Support for both LSM and RDPG models
  \item Comprehensive Evaluation: Evaluate utility through multiple network statistics
  \item Flexible Privacy Budgets: Support for multiple privacy levels in a single run
}

The package is designed for researchers and practitioners working with sensitive network data who need to balance privacy protection with data utility.
}
\author{
\packageAuthor{GRANDpriv}

Maintainer: Suqing Liu <\email{liusuqing0123@uchicago.edu}>
}
\references{
P. D. Hoff, A. E. Raftery, and M. S. Handcock. Latent space approaches to social network analysis.
Journal of the American Statistical Association, 97(460):1090–1098, 2002.

S. J. Young and E. R. Scheinerman. Random dot product graph models for social networks.
In International Workshop on Algorithms and Models for the Web-Graph, pages 138–149. Springer, 2007.

Z. Ma, Z. Ma, and H. Yuan. Universal latent space model fitting for large networks with edge covariates.
Journal of Machine Learning Research, 21(4):1–67, 2020.

A. Athreya, D. E. Fishkind, M. Tang, C. E. Priebe, Y. Park, J. T. Vogelstein, K. Levin, V. Lyzinski, Y. Qin, and D. L. Sussman. Statistical inference on random dot product graphs: a survey.
Journal of Machine Learning Research, 18(226):1–92, 2018.

P. Rubin-Delanchy, J. Cape, M. Tang, and C. E. Priebe. A statistical interpretation of spectral embedding: The generalised random dot product graph.
Journal of the Royal Statistical Society Series B: Statistical Methodology, 84(4):1446–1473, 2022.

X. Bi and X. Shen. Distribution-invariant differential privacy.
Journal of Econometrics, 235(2):444–453, 2023.

S. Liu, X. Bi, and T. Li. GRAND: Graph Release with Assured Node Differential Privacy.
arXiv preprint arXiv:2507.00402, 2025.
}
\keyword{package}
