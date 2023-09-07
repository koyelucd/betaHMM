#' @title Initialising the betaHMM model parameters
#' @description Initialise the betaHMM model parameters.
#' @details Computes the shape parameters using the threshold_function to
#' initialise the shape parameters of the model.
#' @keywords internal
#' @param data A dataframe of dimension \eqn{C \times (N \times R)} containing
#' methylation values for \eqn{C} CpG sites from \eqn{R}
#' treatment groups each having \eqn{N} replicates or each
#' collected from \eqn{N} patients.
#' @param M Number of methylation states to be identified in a DNA sample.
#' @param N Number of patients or DNA sample replicates collected for each
#' treatment group.
#' @param R Number of treatment groups (For. eg: Benign and Tumour).
#' @param seed Seed to allow for reproducibility (default = NULL).
#' @return A list containing:
#' \itemize{
#' \item A - The transition matrix for the betaHMM model.
#' \item tau - The initial distribution for the betaHMM model.
#' \item phi - The shape parameters for the observation sequence data
#' in the betaHMM model.}
#' @importFrom stats kmeans
#' @importFrom stats sd

initialise_parameters <- function(data, M, N, R, seed = NULL) {
    K <- M^R
    A <- 0.6 * diag(K) + rep(0.05/K, K)  ## transmission probabilities
    tau <- rep(1/K, K)  ## initial mixing proportions
    sh_p1 <- list()
    sh_p2 <- list()
    n1 <- 1
    n2 <- N[1]
    data_ilmnid <- data[, "IlmnID"]
    data <- subset(data, select = -c(IlmnID))
    for (r in seq(1, R)) {
        data_th <- cbind(data_ilmnid, data[, n1:n2])
        data_th <- as.data.frame(data_th)
        colnames(data_th)[1] <- "IlmnID"
        threshold_out <- threshold_identification(data_th,
                                                    package_workflow = TRUE,
                                                    M = 3, N[r],
                                                    parameter_estimation_only=
                                                    TRUE, seed = 321)
        model_params <- model_parameters(threshold_out)
        alpha <- model_params$phi$sp_1
        beta <- model_params$phi$sp_2
        alpha <- sort(alpha)
        delta <- sort(beta, decreasing = TRUE)
        alpha[2] <- 1
        delta[2] <- 1
        sh_p1[[r]] <- alpha
        sh_p2[[r]] <- delta
        if ((r + 1) <= R) {
            n1 <- n1 + N[r]
            n2 <- n2 + N[r + 1]
        }
    }
    x1 <- do.call(expand.grid, sh_p1)
    x2 <- do.call(expand.grid, sh_p2)
    x1 <- t(x1)
    x2 <- t(x2)

    phi <- list(sp_1 = x1, sp_2 = x2)
    rownames(phi$sp_1) <- NULL
    rownames(phi$sp_2) <- NULL
    return(list(A = A, tau = tau, phi = phi))
}
