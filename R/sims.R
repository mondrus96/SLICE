#' Generate an exponentially decaying sparse matrix.
#'
#' Returns a matrix where entries decay exponentially with distance from the
#' main diagonal, followed by a random permutation of rows and columns.
#' @param p A numeric. The dimension of the square matrix.
#' @param decay A numeric. The exponential decay rate controlling how quickly
#' values decrease away from the diagonal.
#' @param init_val A numeric. The value placed on the main diagonal.
#' @return A \code{p x p} matrix with exponentially decaying off-diagonal values.
#' @keywords internal
#' @noRd
Smat <- function(p, decay, init_val){
  # Fill the matrix with exponential decay values
  S <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) {
        # Main diagonal element
        S[i, j] <- init_val
      } else {
        # Off-diagonal elements decay exponentially
        distance <- abs(i - j)
        S[i, j] <- init_val * exp(-decay * distance)
      }
    }
  }
  perm <- sample(p) # Permute the matrix
  S <- S[perm, perm]
  
  return(S)
}

#' Generate a random low-rank community matrix.
#'
#' Returns a low-rank matrix constructed from random community assignments
#' together with the corresponding community labels.
#' @param p A numeric. The number of variables (matrix dimension).
#' @param r A numeric. The number of communities (rank).
#' @param init_val A numeric. The scaling value applied to the low-rank matrix.
#' @return A list containing:
#' \item{L}{A \code{p x p} low-rank matrix generated from community memberships.}
#' \item{z}{A vector of community assignments for each variable.}
#' @keywords internal
#' @noRd
Lrand <- function(p, r, init_val){
  probs <- stats::runif(r) # Get probabilities
  probs <- probs/sum(probs)
  
  Z <- matrix(0, p, r)
  for(i in 1:p){
    Z[i, sample(1:r, 1, prob = probs)] <- 1
  }
  z <- apply(Z, 1, function(row) which(row == 1))
  
  L <- Z %*% t(Z)
  L <- L * init_val
  
  return(list(L = L, z = z))
}

#' Simulate multilayer covariance data.
#'
#' Generates sparse and low-rank components along with corresponding covariance
#' matrices and sampled observations for a multilayer graphical model.
#'
#' @export
#'
#' @param r A numeric. The rank of the latent component.
#' @param p A numeric. The number of variables per layer.
#' @param l A numeric. The number of layers.
#' @param n A numeric. The number of observations per layer.
#' @param seed A numeric. Optional random seed for reproducibility.
#' @return A list containing:
#' \item{L_star}{The true latent component matrix.}
#' \item{Sigma_stars}{A list of true covariance matrices for each layer.}
#' \item{Xs}{A list of generated observation matrices for each layer.}
#' \item{Sigmas}{A list of sample covariance matrices computed from \code{Xs}.}
#' \item{S_stars}{A list of true sparse precision components for each layer.}

sim_multislice_data <- function(r = 2, p = 50, l = 2, n = 500, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  pobs <- makeseq(p, l)
  
  # Generate L matrix
  labs <- unique(unlist(pobs))
  Lout <- Lrand(length(labs), r, 1.5)
  ns <- rep(n, l)
  L_star <- Lout$L + 0.1
  
  colnames(L_star) <- rownames(L_star) <- labs
  L_stars <- mat2list(L_star, pobs)
  
  # Generate Sigmas for each layer
  Sigma_stars <- Xs <- Sigmas <- S_stars <- vector("list", l)
  
  for (j in 1:l) {
    S_stars[[j]] <- Smat(length(pobs[[j]]), 2, 1.5)
    S_stars[[j]][S_stars[[j]] < 0.01] <- 0
    
    Sigma_stars[[j]] <- solve(S_stars[[j]] + L_stars[[j]])
    Xs[[j]] <- MASS::mvrnorm(
      n = ns[[j]],
      mu = rep(0, length(pobs[[j]])),
      Sigma = Sigma_stars[[j]]
    )
    Sigmas[[j]] <- stats::cov(Xs[[j]])
  }
  
  return(list(
    L_star = L_star,
    Sigma_stars = Sigma_stars,
    Xs = Xs,
    Sigmas = Sigmas,
    S_stars = S_stars
  ))
}

#' Simulate single-layer covariance data.
#' 
#' Generates sparse and low-rank components along with the corresponding
#' covariance matrix and sampled observations for a single-layer graphical
#' model.
#'
#' @export
#'
#' @param r A numeric. The rank of the latent component.
#' @param p A numeric. The number of variables.
#' @param n A numeric. The number of observations.
#' @param seed A numeric. Optional random seed for reproducibility.
#' @return A list containing:
#' \item{L_star}{The true latent component matrix.}
#' \item{Sigma_star}{The true covariance matrix.}
#' \item{X}{The generated observation matrix.}
#' \item{Sigma}{The sample covariance matrix computed from \code{X}.}
#' \item{S_star}{The true sparse precision component matrix.}

sim_slice_data <- function(r = 2, p = 50, n = 500, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate latent component
  Lout <- Lrand(p, r, 1.5)
  L_star <- Lout$L + 0.1
  
  # Generate sparse component
  S_star <- Smat(p, 2, 1.5)
  S_star[S_star < 0.01] <- 0
  
  # Define true and sample covariance
  Sigma_star <- solve(S_star + L_star)
  X <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = Sigma_star
  )
  Sigma <- stats::cov(X)
  
  return(list(
    L_star = L_star,
    Sigma_star = Sigma_star,
    X = X,
    Sigma = Sigma,
    S_star = S_star
  ))
}