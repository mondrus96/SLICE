#' Matrix completion.
#'
#' Returns completed matrix from missing entries.
#' 
#' @param L A matrix. The matrix to be completed, with missing entries as NA.
#' @param r A numeric. The rank of the joint low-rank space.
#' @param labs A list. The node labels corresponding to each layer.
#' @return A matrix with completed entries.
#' @keywords internal
#' @noRd
matcomp <- function(L, r, labs) {
  # One‐pass truncated SVD per block to build factor H
  p <- nrow(L)
  H <- matrix(0, p, r)
  for (idx in labs) {
    X  <- L[idx, idx]
    sv <- RSpectra::svds(X, r)
    H[idx, ] <- sv$u[, 1:r] %*% diag(sqrt(sv$d[1:r]), r, r)
  }
  
  # Reconstruct low‐rank completion
  return(H %*% t(H))
}

#' Project a joint matrix to list of matrices.
#'
#' Returns a list of matrices given an input joint matrix and 
#' layer-wise labels.
#' @param mat A matrix. The joint matrix to be projected into 
#' a layer-wise list of matrices.
#' @param labs A list. The node labels corresponding to each layer.
#' @return A list of matrices.
#' @keywords internal
#' @noRd
mat2list <- function(mat, labs) {
  uniqlabs <- unique(unlist(labs)) # Get unique labels
  # List to store decomposed matrices
  matlist <- vector("list", length(labs))
  
  # Function to extract relevant matrix part
  extractMat <- function(lab) {
    inds <- match(lab, uniqlabs) # Find indices of labels in uniqlabs
    submat <- mat[inds, inds] # Extract submatrix for given labels
    submat
  }
  
  # Apply extractMat to each set of labels in labs
  matlist <- lapply(labs, extractMat)
  
  return(matlist)
}

#' Project a list of matrices to a joint matrix.
#'
#' Returns a joint matrix given an input list of matrices and 
#' layer-wise labels.
#' 
#' @param matlist A list. The layer-wise list of matrices to be 
#' projected into a joint matrix.
#' @param labs A list. The node labels corresponding to each layer.
#' @return A matrx.
#' @keywords internal
#' @noRd
list2mat <- function(matlist, labs){
  uniqlabs <- unique(unlist(labs)) # Get unique labels
  
  mat <- countmat <- matrix(0, nrow = length(uniqlabs), ncol = length(uniqlabs)) # Create matrices
  labmats <- vector("list", length(matlist)) # Create an empty list to collect index matrices
  
  # Fill labmats with the corresponding index matrices
  for(i in seq_along(matlist)){
    indmat <- matrix(FALSE, nrow = length(uniqlabs), ncol = length(uniqlabs))
    inds <- match(labs[[i]], uniqlabs)
    indmat[inds, inds] <- TRUE
    labmats[[i]] <- indmat
  }
  for(i in seq_along(matlist)) {
    mat[labmats[[i]]] <- mat[labmats[[i]]] + matlist[[i]]
    countmat[labmats[[i]]] <- countmat[labmats[[i]]] + 1
  }
  
  mat <- mat / countmat # Average overlapping values
  mat[is.nan(mat)] <- NA # Replace NaN w/ NA
  colnames(mat) <- rownames(mat) <- uniqlabs
  
  return(mat)
}