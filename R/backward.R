#' @title Backward algorithm
#' @description Compute the backward probabilities for the BHMM model.
#' @details Computes the backward probabilities for the number of hidden states
#'          using the initial distribution, transition matrix and distribution
#'          parameters.
#' @export
#' @param probabilities Probability calculated based on the initialised
#'                      model parameters.
#' @param trained_params A list containing the initialised parameters for the
#'                       BHMM (initial distribution, transition matrix,
#'                      shape parameters for the observation data distribution).
#' @return A list containing:
#' \itemize{
#'    \item log_beta - The backward probabilities.
#'    \item scaled_logL - The log-likelihood calculated using the backward
#'                        probabilities.
#'    }

backward <- function(probabilities,trained_params){
  C = nrow(probabilities)
  K = ncol(probabilities)
  log_beta <- matrix(0, nrow = C, ncol = K)
  log_beta[C,] <- rep(0,K)
  beta <- rep(1 / K, K)
  scaled_logL <- log(K)
  for (i in (C-1):1)
  {
    beta <- trained_params$A %*% (probabilities[i+1,] * beta)
    log_beta[i,] <- log(beta) + scaled_logL
    sum_of_beta <- sum(beta)
    beta <- beta / sum_of_beta
    scaled_logL <- scaled_logL + log(sum_of_beta)
  }

  #backward_probabilities<-list(log_beta=log_beta,scaled_logL=scaled_logL)
  return(list(log_beta=log_beta,scaled_logL=scaled_logL))
}
