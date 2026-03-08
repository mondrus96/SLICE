#' Graphical SCAD estimator.
#'
#' Returns a sparse precision matrix estimate using the SCAD penalty.
#'
#' @param Sigma A matrix. The input covariance matrix.
#' @param lambda A numeric. The regularization parameter.
#' @param a A numeric. The SCAD concavity parameter, default = 3.7.
#' @return A matrix corresponding to the SCAD-penalized precision estimate.
#' @keywords internal
#' @noRd
gscad <- function(Sigma, lambda, a = 3.7){
  # Sigma = input covariance matrix
  # lambda = regularization parameter
  # a = concavity parameter, 3.7 by default

  # Initial L1 estimate
  L1est <- glasso::glasso(Sigma,lambda)$wi

  # Generate rhomat
  p <- nrow(Sigma)
  rhomat <- lambda*matrix(1,p,p)/(pmax(abs(L1est)^0.5,1e-5))

  # Get SCAD estimate
  rhomat <- scadrightderv(abs(L1est), a, lambda)
  SCADest <- glasso::glasso(Sigma, rhomat)$wi
  return(SCADest)
}

# Helper function
scadrightderv <- function(lamhat, a, lambda){
  return(pmax(lambda*((lamhat<=lambda)+pmax(a*lambda-lamhat,0)*(lamhat>lambda)/(a-1)/lambda),1e-4))
}
