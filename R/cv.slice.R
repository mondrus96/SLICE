#' Cross validation for model selection in SLICE
#'
#' This function implements cross validation for SLICE. For full 
#' details, please see the original publication (Ondrus et al, 2024).
#'
#' @export
#'
#' @param X A matrix. The input data matrix.
#' @param folds A numeric. The number of folds to split the data into.
#' @param rhos A vector of numerics. Regularization parameter for sparse estimator.
#' @param rs A vector of integers. Ranks for latent component.
#' @param Sest A character string. Type of sparse estimator to use, default =
#' \code{"glasso"}(Friedman et al, 2008), other choices include \code{"gscad"}
#' (Fan et al., 2009), \code{"clime"}(Cai et al., 2011), and
#' \code{"huge_glasso"} (Zhao et al., 2012).
#' @param tol A numeric. Tolerance for algorithm, default = 1e-3.
#' @param maxiter An integer. Maximum number of iterations, default = 100.
#' @param verbose A logical. Whether to print progress, default = TRUE
#'
#' @details This function implements a grid search over all combinations of
#' \eqn{\rho} and \eqn{r}, and finds the best combination as determined by
#' log likelihood.
#'
#' @return An S3 class \code{cv.slice} object with:
#' \item{cvmat}{A matrix of log likelihood values for each combination of
#' rho and r.}
#' \item{maxlogL}{The maximum log likehood value.}
#' \item{rho}{A numeric of the regularization parameter corresponding to the
#' highest likelihood.}
#' \item{r}{An integer of the rank corresponding to the highest likelihood.}
#'
#' @seealso \code{\link{slice}}
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
#' Ondrus, M., & Cribben, I. A direct method for the estimation 
#' of the sparse and latent variable components of a Gaussian 
#' graphical model. In 2024 10th International Conference on 
#' Control, Decision and Information Technologies (CoDIT), 
#' 522-527, 2024.
#'
#' Zhao, T., Liu, H., Roeder, K., Lafferty, J., and
#' Wasserman, L. The huge package for high dimensional
#' undirected graph estimation in R. \emph{The Journal
#' of Machine Learning Research}, 13(1):1059–1062, 2012.
#'
#' @examples
#' sim_out <- sim_slice_data(r = 4, p = 50, n = 500, seed = 123)
#'
#' X <- sim_out$X
#'
#' out <- cv.slice(X, folds = 3)
#'
#' # Access selected parameters
#' out$rho
#' out$r
cv.slice = function(X, folds = 3, rhos = logseq(1e-5, 0.1, 5), rs = 2:6,
                    Sest = "glasso", tol = 1e-3, maxiter = 100, verbose = TRUE){
  n <- nrow(X) # number of samples
  cvmat <- matrix(NA, length(rs), length(rhos))
  rownames(cvmat) <- rs; colnames(cvmat) <- rhos

  # Go over grid of rhos and rs
  for(i in 1:length(rs)){
    if(verbose){
      print(paste0("rank: ", rs[i]))
    }
    for(j in 1:length(rhos)){
      if(verbose){
        print(paste0("rho: ", rhos[j]))
      }
      ind <- sample(1:folds, n, replace = TRUE) # Define indices
      mulogL <- c()
      for(k in 1:folds){
        train <- X[ind != k,]; test <- X[ind == k,] # Train and test sets

        out <- slice(stats::cov(train), rhos[j], rs[i],
                     Sest, tol = tol, maxiter = maxiter) # Run method
        S <- out$S; L <- out$L

        likl <- logL(stats::cov(test), S + L) # Append to mulogL
        mulogL <- c(mulogL, likl)
      }
      cvmat[i, j] <- mean(mulogL)
    }
  }
  best <- which(cvmat == max(cvmat, na.rm=TRUE), arr.ind = TRUE)
  if(nrow(best) > 2){
    best <- best[1,]
  }
  result <- list(cvmat = cvmat, maxlogL = max(cvmat),
              rho = rhos[best[2]], r = rs[best[1]])
  class(result) <- "cv.slice"
  return(result)
}
