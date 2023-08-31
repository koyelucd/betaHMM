#' @title Backward algorithm
#' @description Computes the backward probabilities for \eqn{K} hidden states
#' using the initial distribution, transition matrix and distribution
#' parameters.
#' @keywords internal
#' @param probabilities A matrix containing the marginal probability
#' distribution of the \eqn{K} states within the Markov model.
#' @param trained_params A list containing the initialised shape parameters for
#' the betaHMM model (initial distribution, transition matrix, shape parameters
#' for the observed data).
#' @return A list containing:
#' \itemize{
#' \item log_beta - A \eqn{C*K} matrix (where \eqn{C} is the number of CpG
#' sites and \eqn{K} is the number of hidden states) containing the
#' logarithmized backward probabilities.
#' \item scaled_logL - The log-likelihood calculated using the backward
#' probabilities.
#' }

backward <- function(probabilities, trained_params) {
    C <- nrow(probabilities)
    K <- ncol(probabilities)
    log_beta <- matrix(0, nrow = C, ncol = K)
    log_beta[C, ] <- rep(0, K)
    beta <- rep(1/K, K)
    scaled_logL <- log(K)
    for (i in (C - 1):1) {
        beta <- trained_params$A %*% (probabilities[i + 1, ] * beta)
        log_beta[i, ] <- log(beta) + scaled_logL
        sum_of_beta <- sum(beta)
        beta <- beta/sum_of_beta
        scaled_logL <- scaled_logL + log(sum_of_beta)
    }
    return(list(log_beta = log_beta, scaled_logL = scaled_logL))
}
