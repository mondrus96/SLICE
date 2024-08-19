#' Log likelihood function
#'
#' This function is used to calculate the log likelihood.
#' @param Sigma A matrix. The masked latent matrix.
#' @param invSigmahat A matrix. The estimated precision matrix.
#' @return The log likelihood.
#' @keywords internal
logL = function(Sigma, invSigmahat){
  if(!isPD(invSigmahat)){
    invSigmahat <- makePD(invSigmahat) # Make PD if not already
  }
  return(log(det(invSigmahat)) - sum(diag(Sigma %*% (invSigmahat))))
}

#' Make a matrix positive definite (PD)
#'
#' This function adds a small value to diagonal to force a matrix to be PD.
#' @param mat A matrix. The input which is potentially non-PD.
#' @return The PD version of the matrix.
#' @keywords internal
makePD = function(mat){
  p = ncol(mat)
  eigvals = suppressWarnings(eigs(mat, ncol(mat), opts = list(retvec = FALSE))$values)
  perturb = max(max(eigvals) - p*min(eigvals), 0)/(p-1)
  mat = mat+diag(p)*perturb
  return(mat)
}

#' Check if a matrix is PD
#'
#' Use Cholesky decomposition to determine if a matrix is PD (faster than full eigendecomp).
#' @param mat A matrix. The input which is potentially non-PD.
#' @return TRUE/FALSE whether the matrix is PD.
#' @keywords internal
isPD = function(mat){
  tryCatch({
    chol(mat)
    return(TRUE)
  }, error = function(e){
    return(FALSE)
  })
}

#' Generate a log sequence.
#'
#' Returns logarthmically spaced sequence of values.
#' @param beg A numeric. The beginning of the sequence
#' @param end A numeric. The end of the sequence
#' @param len A numeric. The number of values to return.
#' @return A vector of logarithmically spaced values.
#' @keywords internal
logseq <- function(beg, end, len) {
  log_beg <- log10(beg); log_end <- log10(end)
  log_seq <- seq(log_beg, log_end, length.out = len)
  return(10^log_seq)
}
