#' SLICE estimator
#'
#' This function implements sparse + low-rank inverse covariance estimation (SLICE).
#'
#' @export
#' @importFrom clime clime
#' @importFrom glasso glasso
#' @importFrom huge huge
#' @importFrom Matrix chol
#' @importFrom Matrix chol2inv
#' @importFrom RSpectra eigs
#' @importFrom RSpectra svds
#'
#' @param Sigma A matrix. The input covariance matrix.
#' @param rho A numeric. Regularization parameter for sparse estimator.
#' @param r An integer. Rank for latent component.
#' @param Sest A character string. Type of sparse estimator to use, default =
#' \code{"glasso"}(Friedman et al, 2008), other choices include \code{"gscad"}
#' (Fan et al., 2009), \code{"clime"}(Cai et al., 2011), and
#' \code{"huge_glasso"} (Zhao et al., 2012).
#' @param tol A numeric. Tolerance for algorithm, default = 1e-3.
#' @param maxiter An integer. Maximum number of iterations, default = 100.
#'
#' @details Given sample covariance \eqn{\boldsymbol{\tilde{\Sigma}}} the
#' objective, for the L1 penalized variant, is to find \eqn{\boldsymbol{
#' \hat{S}}} and \eqn{\boldsymbol{\hat{L}}} which minimize the following function:
#' \deqn{
#' \underbrace{- \mathcal{L}(\boldsymbol{\hat{S}};(\boldsymbol{\tilde{\Sigma}}
#' ^{-1} - \boldsymbol{\hat{L}})^{-1}) + \rho \|\boldsymbol{\hat{S}}\|_1}_{
#' \text{penalized negative log likelihood}} + \underbrace{\|\boldsymbol{\tilde
#' {\Sigma}}(\boldsymbol{\hat{S}} + \boldsymbol{\hat{L}}) - \boldsymbol{I}
#' \|_F^2}_{\text{covariance fidelity}} \\
#' \text{s.t. } \mathcal{R}(\boldsymbol{\hat{L}}) = r, \ \text{where } 0 < r < p
#' }
#' where \eqn{\rho} and \eqn{r} are regularization parameters for the sparse
#' and latent components, respectively.
#'
#' @return An S3 class \code{slice} object with:
#' \item{S}{A matrix corresoponding to the estimated sparse component.}
#' \item{L}{A matrix corresoponding to the estimated latent component.}
#' \item{rho}{A numeric of the regularization
#' parameter used for the sparse component.}
#' \item{r}{An integer of the rank used for the latent component.}
#' \item{misc}{contains additional outputs
#' related to the convergence of the algorithm.}
#'
#' @seealso \code{\link{cv.slice}}
#'
#' @references
#' Cai, T., Liu, W., and Luo, X. A constrained l1 minimization
#' approach to sparse precision matrix estimation. \emph{Journal
#' of the American Statistical Association}, 106(494):594–607, 2011.
#'
#' Fan, J., Feng, Y., and Wu, Y. Network exploration via the
#' adaptive lasso and scad penalties. \emph{The annals of
#' applied statistics}, 3(2):521, 2009.
#'
#' Friedman, J., Hastie, T., and Tibshirani, R. Sparse inverse
#' covariance estimation with the graphical lasso.
#' \emph{Biostatistics}, 9(3):432–441, 2008.
#'
#' Zhao, T., Liu, H., Roeder, K., Lafferty, J., and
#' Wasserman, L. The huge package for high dimensional
#' undirected graph estimation in r. \emph{The Journal
#' of Machine Learning Research}, 13(1):1059–1062, 2012.
#'
#' @examples
#' set.seed(123)
#' p <- 100 # Number of nodes
#'
#' S <- outer(1:p, 1:p,
#'         function(i, j) 1 * exp(-0.5 * abs(i - j))) # Exponential decay
#' S[S < 0.01] <- 0
#' perm <- sample(p) # Permute
#' S <- S[perm, perm]
#'
#' r <- 4 # Rank of latent
#' probs <- runif(r)
#' probs <- probs / sum(probs)
#' Z <- matrix(0, p, r)
#' indices <- sample(1:r, p, replace = TRUE, prob = probs)
#' Z[cbind(1:p, indices)] <- 1
#' L <- Z %*% t(Z)
#'
#' Sigma <- solve(S + L) # Define Sigma
#'
#' out <- slice(Sigma, 0.01, r) # Run SLICE
#'
#' out <- slice(Sigma, 0.01, r, Sest = "gscad") # Run SLICE with SCAD
slice <- function(Sigma, rho, r, Sest = "glasso",
                  tol = 1e-3, maxiter = 100){
  ### Checks ###
  if(!Sest %in% c("glasso", "clime", "gscad", "huge_glasso")){
    stop(paste(Sest, "is not a valid sparse model"))
  }
  if(!isSymmetric(Sigma)){
    stop("Sigma must be symmetric")
  }
  ### Checks ###

  p <- ncol(Sigma) # Make Sigma PD
  Sigma <- makePD(Sigma)
  invSigma <- Matrix::chol2inv(Matrix::chol(Sigma))

  L <- 0 # zero initialization L
  E <- invSigma - L # Expectation
  S <- 0 # Empty S

  deltaS <- deltaL <- deltalogL <- c()
  for(i in 1:maxiter){
    if(!isPD(E)){
      E <- makePD(E) # Make expectation PD
    }

    Sold <- S # Sparse step
    if(Sest == "glasso"){
      S <- glasso(Matrix::chol2inv(Matrix::chol(E)), rho,
                  thr = tol, maxit = maxiter)$wi
    } else if(Sest == "clime"){
      S <- clime(Matrix::chol2inv(Matrix::chol(E)), rho,
                 sigma = TRUE, linsolver = "simplex")$Omegalist[[1]]
      S[abs(S) < tol] <- 0
    } else if(Sest == "gscad"){
      S <- gscad(Matrix::chol2inv(Matrix::chol(E)), rho)
    } else if(Sest == "huge_glasso"){
      S <- huge(Matrix::chol2inv(Matrix::chol(E)), rho,
                method = "glasso", verbose = FALSE)$icov[[1]]
    }
    S <- (S + t(S))/2

    Lold <- L # Latent step
    tsvdL <- svds(invSigma - S, r)
    L <- tsvdL$v %*% diag(tsvdL$d) %*% t(tsvdL$v)
    L <- (L + t(L))/2

    E <- invSigma - L # New expectation
    E <- (E + t(E))/2

    deltaS <- c(deltaS, sqrt(sum((S - Sold)^2))) # Convergence check
    deltaL <- c(deltaL, sqrt(sum((L - Lold)^2)))
    if(i == 1){
      deltalogL <- c(deltalogL,
                     suppressWarnings(abs(logL(Sigma, S + L) -
                                            logL(Sigma, diag(p)*tol))))
    } else{
      deltalogL <- c(deltalogL,
                     suppressWarnings(abs(logL(Sigma, S + L) -
                                            logL(Sigma, Sold + Lold))))
    }
    if((deltaS[i] < tol) && (deltaL[i] < tol) |
       ifelse(is.na(deltalogL[i] < tol), FALSE, deltalogL[i] < tol)){
      break
    }
  }
  # Return result
  result <- list(S = S, L = L, rho = rho, r = r,
              misc = list(converged = ifelse(i < maxiter, TRUE, FALSE), iters = i,
              deltaS = deltaS, deltaL = deltaL, deltalogL = deltalogL))
  class(result) <- "slice"
  return(result)
}
