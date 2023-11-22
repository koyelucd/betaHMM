#' @title Area under the curve method for calculating dissimilarities between
#' estimated distributions
#' @description The function is used to calculate the dissimilarity between the
#' cumulative distributions estimated in each hidden state.
#' @keywords internal
#' @param K Number of hidden states estimated.
#' @param M Number of methylation states to be identified in a single DNA
#' sample.
#' @param N Number of patients or DNA sample replicates collected for each
#' treatment group.
#' @param R Number of treatment groups (For. eg: Benign and Tumour).
#' @param A The transition matrix for hidden states in the betaHMM model.
#' @param tau The initial distribution for hidden states in the betaHMM model.
#' @param phi The shape parameters estimated for the observed data and the
#' estimated hidden states in the betaHMM model.
#' @return A dataframe returning the hidden state, the AUC value and the
#' distributions which resulted in highest AUC value calculated
#' for the corresponding hidden state.
#' @importFrom pROC auc
#' @importFrom stats rbeta
#'
AUC_DM_analysis <- function(M, N, R, K, tau, A, phi) {
    alpha <- t(phi$sp_1)
    delta <- t(phi$sp_2)
    comb <- factorial(R)/(factorial(2) * factorial(R - 2))
    auc_vec <- vector()
    auc_final <- vector()
    comp_final <- vector()
    n <- 1000
    auc_mat <- matrix(NA, nrow = K, ncol = comb)
    text <- vector()
    process_j <- function(j, R, alpha, delta) {
        shape_1 <- alpha[j, seq_len(R), drop = FALSE]
        shape_2 <- delta[j, seq_len(R), drop = FALSE]
        i <- rep(seq_len(R-1), each = R - 1)
        k <- rep((i + 1):R, each = R - 1)
        group_1 <- rbeta(n * (R - 1) * (R - 1), shape1 = shape_1[i],
                         shape2 = shape_2[i])
        group_2 <- rbeta(n * (R - 1) * (R - 1), shape1 = shape_1[k],
                         shape2 = shape_2[k])
        auc_dat <- data.frame(
            predictor = c(group_1, group_2),
            response = factor(rep(0:1, each = n * (R - 1) * (R - 1))))
        auc_values <- auc(predictor = auc_dat$predictor,
                          response = auc_dat$response, quiet = TRUE)
        auc_vec <- unlist(auc_values)
        max_index <- which.max(auc_vec)
        auc_final <- auc_vec[max_index]
        comp_final <- paste(i[max_index], "-->", k[max_index])
        return(c(auc_final,comp_final))
    }
    indices <- seq(1, K)
    result<-lapply(indices, process_j, R = R, alpha = alpha, delta = delta)
    result<-do.call(rbind,result)
    auc_final<-as.numeric(as.vector(result[,1]))
    comp_final<-result[,2]
    auc_df <- as.data.frame(cbind(as.character(seq(1, K)),
                                  as.numeric(round(auc_final, 3)),
                                  comp_final))
    auc_df <- auc_df[order(auc_df$V2), , drop=FALSE]
    auc_df <- t(auc_df)
    rownames(auc_df) <- c("State", "AUC","DM_Conditions")
    return(auc_df)
}
