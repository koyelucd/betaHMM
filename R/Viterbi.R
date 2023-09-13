#' @title Estimating the hidden states using Viterbi algorithm
#' @description The Viterbi algorithm is used to estimate the most probable
#' sequence of the hidden states utilizing the betaHMM model
#' parameters estimated by the Baum-Welch algorithm.
#' @keywords internal
#' @param data A dataframe of dimension \eqn{C \times (N \times R)} containing
#' methylation values for \eqn{C} CpG sites from \eqn{R}
#' treatment groups each having \eqn{N} replicates or each
#' collected from \eqn{N} patients.
#' @param M Number of methylation states to be identified
#' in a single DNA sample.
#' @param N Number of DNA samples (patients/replicates) collected for each
#' treatment group.
#' @param R Number of treatment groups (For. eg: Benign and Tumour).
#' @param A The transition matrix for the betaHMM model.
#' @param tau The initial distribution for the betaHMM model.
#' @param phi The shape parameters for the observation sequence data
#' in the betaHMM model.
#' @return A vector containing the most probable sequence of the hidden states
#' of the betaHMM model.
#' @importFrom stats dbeta
#'

Viterbi <- function(data, M, N, R, tau, A, phi) {
    K <- M^R
    C <- nrow(data)
    probabilities <- list()
    n1 <- 1
    n2 <- N[1]
    for (r in seq(1, R)) {
        for (n in n1:n2) {
            probabilities[[n]] <- matrix(data[, n], ncol = K, nrow=nrow(data))
            probabilities[[n]] <- t(apply(X = probabilities[[n]], MARGIN = 1,
                                            FUN = dbeta,
                                            shape1 = phi$sp_1[r, ],
                                            shape2 = phi$sp_2[r, ]))
        }
        if ((r + 1) <= R) {
            n1 <- n1 + N[r]
            n2 <- n2 + N[r + 1]
        }
    }
    prob <- matrix(probabilities[[1]], nrow = C, ncol = K)
    for (k in 2:ncol(data)) {
        temp <- as.matrix(probabilities[[k]])
        prob <- prob * temp
    }
    probabilities <- matrix(0, ncol = K, nrow = C)
    probabilities <- prob
    omega <- matrix(0, ncol = K, nrow = C)
    prob <- tau * probabilities[1, ]
    omega[1, ] <- prob/sum(prob)
    for (c in 2:C) {
        prob <- apply(omega[c - 1, ] * A, 2, max) * probabilities[c, ]
        omega[c, ] <- prob/sum(prob)
    }
    decoding <- numeric(C)
    decoding[C] <- which.max(omega[C, ])
    for (c in (C - 1):1) {
        decoding[c] <- which.max(A[, decoding[c + 1]] * omega[c, ])
    }
    return(decoding)
}
