#' @title Forward algorithm
#' @description Computes the forward probabilities for \eqn{K} hidden states
#'          using the initial distribution, transition matrix and distribution
#'          parameters.
#' @export
#' @param probabilities A matrix containing the marginal probability
#' distribution of the \eqn{K} states within the Markov model.
#' @param trained_params A list containing the initialised shape parameters for
#' the betaHMM model (initial distribution, transition matrix, shape parameters
#' for the observed data).
#' @return A list containing:
#' \itemize{
#'    \item log_alpha - A \eqn{C*K} matrix (where \eqn{C} is the number of CpG
#'    sites and \eqn{K} is the number of hidden states) containing the
#'    logarithmized forward probabilities.
#'    \item scaled_logL - The log-likelihood calculated using the forward
#'                        probabilities.
#'    }

forward <- function(probabilities, trained_params) {
  C  <-  nrow(probabilities)
  K  <-  ncol(probabilities)
  log_alpha <- matrix(0, nrow = C, ncol = K)
  alpha <- trained_params$tau * probabilities[1,]
  sum_of_alpha <- sum(alpha)
  scaled_logL <- log(sum_of_alpha)
  alpha <- alpha / sum_of_alpha
  log_alpha[1,] <- scaled_logL + log(alpha)
  for (i in 2:C)
  {
    alpha <- alpha %*% trained_params$A * probabilities[i,]
    sum_of_alpha <- sum(alpha)
    scaled_logL <- scaled_logL + log(sum_of_alpha)
    alpha <- alpha / sum_of_alpha
    log_alpha[i,] <- scaled_logL + log(alpha)
  }

  #forward_probabilities<-list(log_alpha=log_alpha,scaled_logL=scaled_logL)
  return(list(log_alpha = log_alpha, scaled_logL = scaled_logL))
}
