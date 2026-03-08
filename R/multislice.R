#' multiSLICE estimator
#'
#' This function implements multilayer sparse + low-rank inverse covariance 
#' estimation (multiSLICE) for multilayer graphical models. For full details, 
#' please see the original publication (Ondrus et al, 2025).
#' 
#' @export
#'
#' @param Sigmas A list. The set of input covariance matrices, where each
#' element in the list is a matrix.
#' @param rhos Either a numeric (same penalty across all layers) or a vector
#' of numeric values (layer-wise penalty). Regularization parameter(s) for 
#' sparse estimators.
#' @param r An integer. Rank for latent component.
#' @param Sest A character string. Type of sparse estimator to use, default =
#' \code{"glasso"}(Friedman et al, 2008), other choices include \code{"gscad"}
#' (Fan et al., 2009), \code{"clime"}(Cai et al., 2011), and
#' \code{"huge_glasso"} (Zhao et al., 2012).
#' @param tol A numeric. Tolerance for algorithm, default = 1e-3.
#' @param maxiter An integer. Maximum number of iterations, default = 100.
#'
#' @details Given sample covariances \eqn{\{\boldsymbol{\tilde{\Sigma}}\}_{
#' \boldsymbol{\alpha}=1}^l}the objective, for the L1 penalized variant, is 
#' to find \eqn{\{\boldsymbol{\hat{S}}\}_{\boldsymbol{\alpha}=1}^l} and \eqn{
#' \boldsymbol{\hat{L}}} which minimize the following function:
#' \deqn{
#' \underbrace{\sum_{\boldsymbol{\alpha}=1}^{l}}_{\text{layers}} \left(
#' \underbrace{- \mathcal{L}(\boldsymbol{\hat{S}}_{\boldsymbol{\alpha}};
#' (\boldsymbol{\tilde{\Sigma}}_{\boldsymbol{\alpha}}^{-1} -
#' \boldsymbol{\hat{L}}_{\boldsymbol{\alpha}})^{-1}) + \rho
#' \|\boldsymbol{\hat{S}}_{\boldsymbol{\alpha}}\|_1}_{\text{penalized negative
#' log likelihood}} + \underbrace{\|\boldsymbol{\tilde{\Sigma}}_{\boldsymbol{
#' \alpha}}(\boldsymbol{\hat{S}}_{\boldsymbol{\alpha}} + \boldsymbol{\hat{L}}_{
#' \boldsymbol{\alpha}}) - \boldsymbol{I}\|_F^2}_{\text{covariance fidelity}}
#' \right) \\
#' \text{s.t. } \mathcal{R}(\boldsymbol{\hat{L}}) = r, \ \text{where } 0 < r < p
#' }
#' where \eqn{\rho} and \eqn{r} are regularization parameters for the sparse
#' and latent components, respectively.
#'
#' @return An S3 class \code{multislice} object with:
#' \item{Ss}{A list of matrices corresoponding to the estimated sparse components.}
#' \item{L}{A matrix corresoponding to the estimated latent component.}
#' \item{rhos}{A vector of the regularization parameter(s) used for the 
#' sparse components.}
#' \item{r}{An integer of the rank used for the latent component.}
#' \item{misc}{contains additional outputs
#' related to the convergence of the algorithm.}
#'
#' @seealso \code{\link{cv.multislice}}
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
#' Ondrus, M., Cribben, I., & Feng, Y. A Latent Multilayer 
#' Graphical Model For Complex, Interdependent Systems. \emph{The 
#' Thirty-ninth Annual Conference on Neural Information 
#' Processing Systems}, 2025.
#'
#' Zhao, T., Liu, H., Roeder, K., Lafferty, J., and
#' Wasserman, L. The huge package for high dimensional
#' undirected graph estimation in R. \emph{The Journal
#' of Machine Learning Research}, 13(1):1059–1062, 2012.
#'
#' @examples
#' sim_out <- sim_multislice_data(r = 2, p = 50, l = 2, n = 500, seed = 123)
#'
#' Sigmas <- sim_out$Sigmas
#' r <- 2
#'
#' out <- multislice(Sigmas, 0.01, r) # Run multiSLICE
#'
#' out <- multislice(Sigmas, 0.01, r, Sest = "gscad") # Run multiSLICE with SCAD
multislice <- function(Sigmas, rhos, r, Sest = "glasso",
                  tol = 1e-3, maxiter = 100){
  ### Checks ###
  if(!Sest %in% c("glasso", "clime", "gscad", "huge_glasso")){
    stop(paste(Sest, "is not a valid sparse model"))
  }
  if (!all(sapply(Sigmas, isSymmetric))) {
    stop("All Sigma matrices must be symmetric")
  }
  ### Checks ###
  
  if(length(rhos) == 1){
    rhos = rep(rhos, length(Sigmas))
  }
  
  # Apply SLICE
  labs <- lapply(Sigmas, nrow)
  labs <- Map(seq, cumsum(c(1, utils::head(labs, -1))), cumsum(labs))
  slices <- mapply(slice, Sigmas, rhos, MoreArgs = list(r = r,
                                                        Sest = Sest,
                                                        tol = tol,
                                                        maxiter = maxiter), SIMPLIFY = FALSE) # Apply independent slice models
  Ls <- lapply(slices, "[[", "L"); Ss <- lapply(slices, "[[", "S") # Get Ls and Ss
  Ss <- Map(function(matrix, labels){ # For relabeling Ss with actual labs
    colnames(matrix) <- labels
    rownames(matrix) <- labels
    matrix
  }, Ss, labs)
  
  L <- list2mat(Ls, labs) # Bring together Ls
  L <- matcomp(L, r, labs)
  
  # Return result
  slices_misc <- lapply(slices, "[[", "misc")
  result <- list(S = Ss, L = L, rhos = rhos, r = r,
                 misc = list(converged =
                               sapply(slices_misc, "[[", "converged",
                                        simplify = FALSE),
                             iters = sapply(slices_misc, "[[", "iters", 
                                              simplify = FALSE),
                             deltaS = sapply(slices_misc, "[[", "deltaS",
                                             simplify = FALSE),
                             deltaL = sapply(slices_misc, "[[", "deltaL",
                                             simplify = FALSE)))
  return(result)
}