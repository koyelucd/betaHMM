#' @title Area under the curve method for calculating dissimilarities between
#' estimated distributions
#' @description The function is used to calculate the dissimilarity between the
#' cumulative distributions estimated in each hidden state.
#' @export
#' @param K Number of hidden states estimated
#' @param M Number of methylation states to be identified in a single DNA sample
#' @param N Number of patients or DNA sample replicates collected for each
#'          treatment group.
#' @param R Number of treatment groups (For. eg: Benign and Tumour).
#' @param A The transition matrix for hidden states in the betaHMM model.
#' @param tau The initial distribution for hidden states in the betaHMM model.
#' @param phi The shape parameters estimated for the observed data and the
#'            estimated hidden states in the betaHMM model.
#' @return A dataframe returning the hidden state and the AUC value calculated
#'  for the corresponding hidden states
#' @importFrom pROC auc
#' @importFrom stats rbeta
#'
AUC_DM_analysis<-function(M, N, R, K, tau, A, phi)
{

  alpha <- t(phi$sp_1)
  delta <- t(phi$sp_2)
  comb <- factorial(R) / (factorial(2) * factorial(R - 2))
  auc_vec <- vector()
  auc_final <- vector()
  n <- 1000
  auc_mat <- matrix(NA, nrow = K, ncol = comb)
  text <- vector()
  for(j in 1:K)
  {
    shape_1 <- alpha[j, 1:R]
    shape_2 <- delta[j, 1:R]
    vec_count=1
    for(i in 1: (R-1))
    {
      for(k in (i+1):R)
      {
        group_1 <- stats::rbeta(n, shape1 = shape_1[i], shape2 = shape_2[i])
        group_2 <- stats::rbeta(n, shape1 = shape_1[k], shape2 = shape_2[k])
        auc_dat <- data.frame(predictor=c(group_1, group_2),
                              response=factor(c(rep(0,n), rep(1,n))))
        auc_value <- pROC::auc(predictor = auc_dat$predictor,
                               response  = auc_dat$response,quiet=TRUE)
        auc_vec[vec_count] <- unlist(auc_value)
        vec_count <- vec_count+1
      }
    }
    auc_final[j] <- max(auc_vec)
  }

  auc_df <- as.data.frame(cbind(as.character(seq(1:K)),
                              as.numeric(round(auc_final,3))))
  auc_df <- auc_df[order(auc_df$V2),]
  auc_df <- t(auc_df)
  rownames(auc_df) <- c("State","AUC")
  return(auc_df)
}
