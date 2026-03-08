#' Log likelihood function
#'
#' This function is used to calculate the log likelihood.
#' 
#' @param Sigma A matrix. The masked latent matrix.
#' @param invSigmahat A matrix. The estimated precision matrix.
#' @return The log likelihood.
#' @keywords internal
#' @noRd
logL = function(Sigma, invSigmahat){
  if(!isPD(invSigmahat)){
    invSigmahat <- makePD(invSigmahat) # Make PD if not already
  }
  return(log(det(invSigmahat)) - sum(diag(Sigma %*% (invSigmahat))))
}

#' Make a matrix positive definite (PD)
#'
#' This function adds a small value to diagonal to force a matrix to be PD.
#' 
#' @param mat A matrix. The input which is potentially non-PD.
#' @return The PD version of the matrix.
#' @keywords internal
#' @noRd
makePD = function(mat){
  p = ncol(mat)
  eigvals = suppressWarnings(RSpectra::eigs(mat, ncol(mat), opts = list(retvec = FALSE))$values)
  perturb = max(max(eigvals) - p*min(eigvals), 0)/(p-1)
  mat = mat+diag(p)*perturb
  return(mat)
}

#' Check if a matrix is PD
#'
#' Use Cholesky decomposition to determine if a matrix 
#' is PD (faster than full eigendecomp).
#' 
#' @param mat A matrix. The input which is potentially non-PD.
#' @return TRUE/FALSE whether the matrix is PD.
#' @keywords internal
#' @noRd
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
#' 
#' @param beg A numeric. The beginning of the sequence.
#' @param end A numeric. The end of the sequence.
#' @param len A numeric. The number of values to return.
#' @return A vector of logarithmically spaced values.
#' @keywords internal
#' @noRd
logseq <- function(beg, end, len) {
  log_beg <- log10(beg); log_end <- log10(end)
  log_seq <- seq(log_beg, log_end, length.out = len)
  return(10^log_seq)
}

#' Generate index sequences for layers.
#'
#' Returns a list of consecutive index sequences of length \code{p} for each
#' layer.
#' 
#' @param p A numeric. The length of each sequence.
#' @param layers A numeric. The number of layers (sequences) to generate.
#' @return A list where each element is a vector of consecutive indices
#' corresponding to a layer.
#' @keywords internal
#' @noRd
makeseq <- function(p, layers){
  seqs <- vector("list", layers)  # Initialize the list to store sequences
  start <- 1  # Starting index for the first sequence
  
  for (i in 1:layers) {
    end <- start + p - 1  # Calculate the ending index
    seqs[[i]] <- start:end  # Assign the sequence to the list
    start <- end + 1  # Update start for the next sequence
  }
  
  return(seqs)
}