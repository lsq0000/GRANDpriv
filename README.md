# GRANDpriv: Graph Release with Assured Node Differential Privacy

![Version](https://img.shields.io/badge/version-v0.1.3-blue.svg) [![CRAN](https://img.shields.io/badge/CRAN-available-brightgreen.svg)](https://CRAN.R-project.org/package=GRANDpriv) [![arXiv](https://img.shields.io/badge/arXiv-2507.00402-b31b1b.svg)](https://arxiv.org/abs/2507.00402)

## Overview

**GRANDpriv** (**G**raph **R**elease with **A**ssured **N**ode **D**ifferential **priv**acy) is an R package that implements a novel method for privatizing network data using differential privacy. The package provides functions for generating synthetic networks based on LSM (Latent Space Model), applying differential privacy to network latent positions to achieve overall network privatization, and evaluating the utility of privatized networks through various network statistics. The privatize and evaluate functions support both LSM and RDPG (Random Dot Product Graph). For generating RDPG networks, users are encouraged to use the [`randnet`](https://CRAN.R-project.org/package=randnet) package. For more details, see the "proposed method" section of [Liu, Bi, and Li (2025)](https://arxiv.org/abs/2507.00402).

## Installation

### From CRAN (Stable Version)

```r
# Install GRANDpriv from CRAN
install.packages("GRANDpriv")
```

<!-- *Note: CRAN version is currently under review and not yet available. Please use the development version from GitHub.* -->

### From GitHub (Development Version)

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install GRANDpriv from GitHub
devtools::install_github("lsq0000/GRANDpriv")
```

### Dependencies

The package requires the following R packages:
- EnvStats
- rmutil  
- RSpectra
- diffpriv
- truncnorm
- randnet
- igraph
- HCD
- transport

## Quick Start

```r
library(GRANDpriv)

# Generate a sample network using Latent Space Model
# Note: Use larger networks (n >= 1000) for better stability
network <- LSM.Gen(n = 400, k = 2, K = 3, avg.d = 40)

# Privatize the first 200 nodes with different privacy budgets
result <- GRAND.privatize(
  A = network$A, 
  K = 2, 
  idx = 1:200, 
  eps = c(1, 2, 5, 10), 
  model = "LSM"
)

# Evaluate the privatization results
evaluation <- GRAND.evaluate(result)
print(evaluation)

# Evaluate specific statistics only
evaluation_subset <- GRAND.evaluate(result, statistics = c("degree", "triangle"))
print(evaluation_subset)
```

## Main Functions

### 1️⃣ `LSM.Gen(n, k, K, avg.d = NULL)`
Generates a random network following LSM (Latent Space Model) with specified parameters.

- `n`: ***Integer***. Number of nodes in the network.
- `k`: ***Integer***. Dimension of the latent space.
- `K`: ***Integer***. Number of communities/groups.
- `avg.d`: ***Numeric***. Target average node degree. If NULL, no degree adjustment is performed.

**Returns**: A list containing:
- `A`: Adjacency matrix of the generated network
- `g`: Graph object of the generated network
- `P`: Probability matrix of the generated network
- `alpha`: Node-specific intercept parameters
- `Z`: Latent positions in k-dimensional space
- `idx`: Community assignments for each node

### 2️⃣ `GRAND.privatize(A, K, idx, eps = 1, model = "LSM", niter = 500, rho = 0.05, verbose = TRUE)`
Applies the GRAND (Graph Release with Assured Node Differential privacy) method to privatize network data using differential privacy.

- `A`: ***Matrix***. Adjacency matrix of the input network.
- `K`: ***Integer***. Dimension of the latent space for network embedding.
- `idx`: ***Integer vector***. Indices of nodes to be privatized.
- `eps`: ***Numeric or numeric vector***. Privacy budget parameter(s) for differential privacy. Default is 1.
- `model`: ***Character***. Model type, either "LSM" (Latent Space Model) or "RDPG" (Random Dot Product Graph). Default is "LSM".
- `niter`: ***Integer***. Number of iterations for the optimization algorithm. Default is 500.
- `rho`: ***Numeric***. Parameter controlling the neighborhood size for conditional distributions. Default is 0.05.
- `verbose`: ***Logical***. Whether to print progress messages. Default is TRUE.

**Returns**: A list containing:
- `non.private.result`: Results without privacy (original and estimated data)
- `GRAND.result`: List with one element per epsilon value. Each element contains privatization results for that specific epsilon
- `Laplace.result`: List with one element per epsilon value. Each element contains baseline Laplace mechanism results for that specific epsilon
- `eps`: Vector of privacy budget parameter(s) used

### 3️⃣ `GRAND.evaluate(result, statistics = c("degree", "vshape", "triangle", "eigen", "harmonic"))`
Evaluates the quality of GRAND privatization results by comparing various network statistics between the original and privatized networks using Wasserstein distance.

- `result`: ***List***. Output from GRAND.privatize function containing privatization results.
- `statistics`: ***Character vector***. Network statistics to evaluate. Options include: "degree" (Node Degree), "vshape" (V-Shape Count), "triangle" (Triangle Count), "eigen" (Eigen Centrality), "harmonic" (Harmonic Centrality). Default is all statistics.

**Returns**: A data frame containing evaluation results with columns:
- `stat`: Type of network statistic(s) evaluated
- `eps`: Privacy budget parameter(s) used
- `Hat`: Wasserstein distance for release set estimation
- `Hat2`: Wasserstein distance for holdout set estimation
- `GRAND`: Wasserstein distance for GRAND privatization method
- `Laplace`: Wasserstein distance for standard Laplace mechanism

## Features

- **Network Generation**: Generate synthetic networks using Latent Space Model
- **Differential Privacy**: Apply node-level differential privacy to network data
- **Multiple Models**: Support for both LSM and RDPG models
- **Comprehensive Evaluation**: Evaluate utility through multiple network statistics
- **Flexible Privacy Budgets**: Support for multiple privacy levels in a single run

## Methodology

GRAND uses a two-step approach:
1. **Latent Position Estimation**: Estimates latent positions from the network structure using either LSM.PGD (Projected Gradient Descent, for LSM) or ASE (Adjacency Spectral Embedding, for RDPG)
2. **Multivariate Differential Privacy**: Applies DIP (Distribution-Invariant differential Privacy) mechanism to protect latent positions while preserving network utility

## License

GPL (>= 3)

## Citation

If you use **GRANDpriv** in your research, please cite:

```bibtex
@misc{liu2025grand,
  title={GRAND: Graph Release with Assured Node Differential Privacy},
  author={Suqing Liu and Xuan Bi and Tianxi Li},
  year={2025},
  eprint={2507.00402},
  archivePrefix={arXiv},
  primaryClass={stat.ML},
  url={https://arxiv.org/abs/2507.00402}
}
```

Or in text format:
```
S. Liu, X. Bi, and T. Li. GRAND: Graph Release with Assured Node Differential Privacy. arXiv preprint arXiv:2507.00402, 2025.
```

## Issues and Contributions

Please report issues or contribute to the package through GitHub.
