#' @title Threshold indentification function
#' @description HMM for beta valued DNA data for a single treatment condition
#'
#' @export
#' @param data Methylation data and IlmnID. Maybe provided as
#'  a matrix or dataframe.
#' @param package_workflow Flag set to TRUE if method called from package
#'   workflow.
#'   If set to FALSE then the parameter annotation_file needs to be supplied.
#' @param annotation_file A dataframe containing the EPIC methylation
#' annotation file.
#'@param M Number of methylation states to be identified in a single DNA sample.
#' @param N Number of DNA samples (patients/replicates) collected for each
#'          treatment group.
#' @param parameter_estimation_only If only model parameters are to be estimated
#'                                  then value is TRUE else FALSE.
#' @param seed Seed to allow for reproducibility (default = NULL).
#' @param ... Extra parameters.
#' @return  An S4 object of class \code{threshold_Results}, where conditional
#' probabilities
#'of each CpG site belonging to a one of the \eqn{M} methylation states
#'is stored as a SimpleList of
#'assay data, and the corresponding estimated model parameters, the thresholds
#'and most probable hidden state sequence for each chromosome are
#'stored as metadata.
#'
#' @importFrom stats kmeans
#' @importFrom stats sd
#' @examples
#' ## Use simulated data for the betaHMM workflow example
#' set.seed(12345)
#' ## read methylation data
#' methylation_data=sample_methylation_file
#'
#' ## read sample EPIC annotation file
#' legacy_data=sample_annotation_file
#'
#' ##merge data
#' df=merge(legacy_data[,c("IlmnID","CHR","MAPINFO")],
#' methylation_data,by="IlmnID")
#'
#' ## sort data
#' df=df[order(df$CHR,df$MAPINFO),]
#' thr_out=threshold_identification(df[,c(1,4:7)],package_workflow=TRUE,M=3,4,
#' parameter_estimation_only=TRUE,seed=12345)
#'

threshold_identification_run <- function(data, package_workflow = TRUE,
                                         annotation_file = NULL,
                                         M, N,
                                         parameter_estimation_only = FALSE,
                                         seed = NULL, ...)
{


  if(package_workflow == FALSE)
  {
    if(is.null(annotation_file))
      stop("Annotation file cannot be empty if data is not sorted.")

    final_subset <- subset(data, select = - IlmnID)
    #
    if(abs(max(final_subset) - 1) < 1e-6)
    {
      max_f <- max(final_subset[final_subset != max(final_subset)])
    }else{
      max_f <- max(final_subset)
    }
    if(abs(min(final_subset) - 0) < 1e-6)
    {
      min_f <- min(final_subset[final_subset != min(final_subset)])
    }else{
      min_f <- min(final_subset)
    }

    final_subset[final_subset > max_f] <- max_f
    final_subset[final_subset < min_f] <- min_f
    cols <- which(colnames(data) != "IlmnID")
    data[ ,cols] <- final_subset
    data_merged <- merge(data, annotation_file[ ,c("IlmnID", "CHR", "MAPINFO")],
                         by = "IlmnID")
    col_order <- which(colnames(data_merged) == "IlmnID" |
                       colnames(data_merged) == "CHR" |
                       colnames(data_merged) == "MAPINFO")
    col_order2 <- which(colnames(data_merged) != "IlmnID" &
                        colnames(data_merged) != "CHR" &
                        colnames(data_merged) != "MAPINFO")
    data_merged <- data_merged[ ,c(col_order, col_order2)]
    complete_data <- data_merged
    complete_data <- complete_data[stats::complete.cases(complete_data), ]
    sorted_data <- complete_data[order(as.numeric(complete_data$CHR),
                                     complete_data$MAPINFO), ]
    sorted_data <- as.data.frame(sorted_data)
    data_final <- subset(sorted_data, select = -c(CHR, MAPINFO, IlmnID))
    data_return <- subset(sorted_data, select = -c(CHR, MAPINFO))

  }else{
    data_final <- subset(data, select = -c( IlmnID))
    data_return <- data

  }
initialize_parameters_th <- function(data, M, N, R, seed = NULL){
  K <- M
  # K=M
  A <- 0.6 * diag(K) + rep(0.05 / K, K)  ## transmission probabilities
  tau <- rep(1 / K, K)                  ## initial mixing proportions
  data <- as.data.frame(data)
  for(i in 1:ncol(data))
  {
    data[ ,i] <- as.numeric(data[ ,i])
  }
  data <- as.matrix(data)
  k_cluster <- stats::kmeans(data, M)
  mem <- k_cluster$cluster
  data_clust <- cbind(data, mem)
  x <- as.matrix(data_clust)
  mem <- x[ ,ncol(x)]
  data_full <- x
  x <- x[ ,-ncol(x)]
  C <- nrow(x)
  N <- ncol(x)

  ## initial model parameters
  mu <- vector("numeric", K)
  sigma <- vector("numeric", K)
  sigma_sq <- vector("numeric", K)
  alpha <- vector("numeric", K)
  beta <- vector("numeric", K)
  term <- vector("numeric", K)


  for(k in 1:K)
  {
    mu[k] <- mean(x[mem == k, ])
    sigma[k] <- stats::sd(x[mem == k, ])
    term[k] <- (mu[k] * (1 - mu[k]) / (sigma[k] ^ 2)) - 1
    alpha[k] <- mu[k] * term[k]
    beta[k] <- (1 - mu[k]) * term[k]
  }
  alpha <- sort(alpha)
  delta <- sort(beta, decreasing = T)
  alpha[2] <- 1
  delta[2] <- 1


  phi <- list(sp_1 = alpha, sp_2 = delta)

  return(list(A = A, tau = tau, phi = phi )) #,thresholds=threshold_mat))
}

BaumWelch_th  <-  function(data, M = 3, N, R, n.iter = 100, seed = NULL){
  K <- M
  C  <-  nrow(data)


  ## Getting initialized parameters
  trained_params <- initialize_parameters_th(data, M, N, R, seed)

  #v=xmat
  oldlogL <- 1
  BW_limit_accuracy  <-  1e-5
  logL_vec <- vector()
  logL_vec <- oldlogL
  R <- 1
  for(i in 1:n.iter){

    probabilities <-  list()



    for(n in 1:N)
    {
      probabilities[[n]] <- matrix(data[,n], ncol = K, nrow = nrow(data))
      probabilities[[n]] <- t(apply(X = probabilities[[n]], MARGIN = 1,
                                    FUN = dbeta,
                                    shape1 = trained_params$phi$sp_1,
                                    shape2 = trained_params$phi$sp_2))

    }

    prob <- matrix(probabilities[[1]], nrow = C, ncol = K)
    for(k in 2:(N))
    {
      temp <- as.matrix(probabilities[[k]])
      prob <- prob * temp
    }
    probabilities <- matrix(0, ncol = K, nrow = C)
    probabilities <- prob


    ## Generating alpha and beta matrices
    forward_alpha  <-  forward(probabilities, trained_params)
    backward_beta  <-  backward(probabilities, trained_params)

    ## Calculating log-likelihood
    log_alpha <- forward_alpha$log_alpha
    logL_calculated_with_alpha_T <- forward_alpha$scaled_logL

    log_beta <- backward_beta$log_beta
    logL_calculated_with_beta_1 <- backward_beta$scaled_logL

    logL_calculated_with_alpha_t_and_beta_t <- 0
    middle_t <- round(C / 2)
    t <- max(log_alpha[middle_t, ])
    for (i in 1:K)
    {
      logL_calculated_with_alpha_t_and_beta_t <-
        logL_calculated_with_alpha_t_and_beta_t +
        exp( log_alpha[middle_t,i] + log_beta[middle_t,i]) #- t)
    }
    logL_calculated_with_alpha_t_and_beta_t <-
      log(logL_calculated_with_alpha_t_and_beta_t) #+ t

    logL <- logL_calculated_with_alpha_t_and_beta_t
    logL_calculation <- "logL calculated with alpha_middle_t and beta_middle_t"
    if (logL == -Inf | logL == Inf | is.na(logL))
    {
      logL_calculation <- "logL calculated with alpha_T"
      logL <- logL_calculated_with_alpha_T
    }
    if (logL == -Inf | logL == Inf | is.na(logL))
    {
      logL <- logL_calculated_with_beta_1
      logL_calculation <- "logL calculated with beta_1"
    }

    ## Calculating the joint and conditional probabilities

    eta <- matrix(c(0), ncol = K, nrow = C)
    eta <- exp(log_alpha + log_beta - logL)


    xi <- array(NA, dim = c( (C-1), K, K))

    alpha_c <- log_alpha[1:C-1, ]
    beta_c <- log_beta[2:C, ]
    proba <- log(probabilities[2:C, ])
    for (j in 1:K)
    {
      for (i in 1:K)
      {

        xi[ ,i , j] <- exp( alpha_c[ ,i] + log(trained_params$A[i,j]) +
                             proba[ ,j]+ beta_c[,j] - logL)
      }
    }

    ## Initial distribution proportions
    tau <- eta[1, ]



    ## Estimating transmission probabilities
    A <- matrix(0, K, K)

    A1 <- apply(xi, c(2,3), sum)
    A2 <- apply(A1, FUN = sum, MARGIN = 1)
    A <- A1/A2


    # registerDoSEQ()
    #
    y1 <- vector("logical", K)
    y2 <- vector("logical", K)
    al_new2 <- vector("logical", K)
    be_new2 <- vector("logical", K)
    term <- vector("logical", K)

    if(N == 1)
    {
      y1 <- apply((eta * (log(data[ ,1:N]))), FUN = sum, MARGIN = 2)/
        (N * apply(eta, FUN = sum, MARGIN = 2))
      y2 <- apply((eta * (log(1 - data[ ,1:N]))), FUN = sum, MARGIN = 2)/
        (N * apply(eta, FUN = sum, MARGIN = 2))
    }else{
      y1 <- apply((eta * rowSums(log(data[ ,1:N]))), FUN = sum, MARGIN = 2)/
        (N * apply(eta, FUN = sum, MARGIN = 2))
      y2 <- apply((eta * rowSums(log(1 - data[ ,1:N]))), FUN = sum, MARGIN = 2)/
        (N * apply(eta, FUN = sum, MARGIN = 2))
    }
    term <- ((exp(-y1) - 1) * (exp(-y2) - 1)) - 1
    al_new2 <- 0.5 + (0.5 * exp(-y2) / term)
    be_new2 <- (0.5 * exp(-y2) * (exp(-y1) - 1)) / term




    phi <- list(sp_1 = al_new2, sp_2 = be_new2)


    trained_params$A <- A
    trained_params$tau <- tau
    trained_params$phi <- phi


    ## Comparing log-likelihood calculated with the previous log-likhelihood
    ## Exiting recurssion if the difference is below the selected epsilon
    logL_vec <- c(logL_vec, logL)
    difference_old_logL_and_new_logL  <-  (logL - oldlogL) / oldlogL

    if (difference_old_logL_and_new_logL > 0 &
        difference_old_logL_and_new_logL < BW_limit_accuracy)
    {
      reached_limit_of_accuracy  <-  TRUE
      break
    }

    oldlogL <- logL
  }
  log_vec <- logL_vec[-1]
  return(list(A = A,  tau = tau, phi = phi, log_vec = logL_vec, z = eta))
}
Viterbi_th <- function(data, M = 3, N = 4, R = 1, tau, A, phi) {
  K <- M
  C  <-  nrow(data)

  ## Generating the probabilities Pr(Y_c | X_c=k)
  probabilities <-  list()



  for(n in 1:N)
  {
    probabilities[[n]] <- matrix(data[,n], ncol = K, nrow = nrow(data))
    probabilities[[n]] <- t(apply(X = probabilities[[n]], MARGIN = 1,
                                  FUN = dbeta, shape1 = phi$sp_1,
                                  shape2 = phi$sp_2))

  }

  prob <- matrix(probabilities[[1]], nrow = C, ncol = K)
  for(k in 2:(N))
  {
    temp <- as.matrix(probabilities[[k]])
    prob <- prob * temp
  }
  probabilities <- matrix(0, ncol = K, nrow = C)
  probabilities <- prob

  omega <- matrix(0, ncol = K, nrow = C)
  prob <- tau * probabilities[1, ]
  omega[1,] <- prob / sum(prob)
  for(c in 2:C)
  {
    prob <- apply(omega[c-1,] * A, 2, max) * probabilities[c,]
    omega[c,] <- prob / sum(prob)
  }

  decoding <- numeric(C)
  decoding[C] <- which.max(omega[C,])
  for(c in (C - 1):1)
  {
    decoding[c] <- which.max(A[,decoding[c+1]] * omega[c,])
  }



  return(decoding)

}

threshold <- function(data, tau, phi){

  threshold_func_low <- function(data_x, alpha, delta, tau, i)
  {
    num <- tau[i] * stats::dbeta(data_x, alpha[i], delta[i])
    deno <- 0
    for(j in 1:length(alpha))
    {
      if(j != i)
      {

        deno <- deno + (tau[j] * stats::dbeta(data_x, alpha[j], delta[j]))
      }
    }

    r <- num / deno
    index <- which(r >= 1)
    l <- data_x[min(index)]
    return(l)
  }
  threshold_func_upper <- function(data_x, alpha, delta, tau, i)
  {
    num <- tau[i] * stats::dbeta(data_x, alpha[i], delta[i])
    deno <- 0
    for(j in 1:length(alpha))
    {
      if(j != i)
      {
        deno <- deno + (tau[j] * stats::dbeta(data_x, alpha[j], delta[j]))
      }
    }

    r <- num / deno
    index <- which(r >= 1)
    u <- data_x[max(index)]
    return(u)
  }
  data_x <- sort(data[ ,1])
  mode <- (phi$sp_1 - 1) / (phi$sp_1 + phi$sp_2 - 2)
  cluster <- c(1, 2, 3)
  mode_vec <- cbind(mode, cluster)
  mode_vec <- as.data.frame(mode_vec)
  mode_vec <- mode_vec[order(mode_vec$mode), ]
  hypo <- as.numeric(mode_vec[1,2])
  hyper <- as.numeric(mode_vec[3,2])
  th_1 <- threshold_func_upper(data_x, phi$sp_1, phi$sp_2, tau, hypo)
  th_2 <- threshold_func_low(data_x, phi$sp_1, phi$sp_2, tau, hyper)
  th_vec <- c(th_1,th_2)
  th_new1 <- unique(round(th_vec, 3))

  return(list(thresholds = th_new1))
}


 th_bw_out <- BaumWelch_th(data_final, M = M, N = N, R = 1, seed = seed)
 if(parameter_estimation_only == FALSE)
 {
   th_vit_out <- Viterbi_th(data_final, M = M, N = N, R = 1, th_bw_out$tau,
                         th_bw_out$A, th_bw_out$phi)
   tau <- as.numeric(table(th_vit_out) / (nrow(data_final)))
   threshold_bhmm <- threshold(data_final, tau, th_bw_out$phi)
 }else{
   th_vit_out <- NA
   threshold_bhmm <- NA
 }

 # threshold_out<-list(model_params=th_bw_out,states=th_vit_out,
 #                     threshold=threshold_bhmm)
 # class(threshold_out)<-"threshold_identification"
 # return(threshold_out)

 ####################################
 ## RETURN RESULTS
 ####################################
 z <- as.data.frame(th_bw_out$z)
 rownames(z) <- data_return[ ,"IlmnID"]
 data_return <- as(data_return, "DataFrame")
 select.results <- SummarizedExperiment(assays = z,
                                        metadata = list(model_parameters =
                                                          th_bw_out, K = M,
                                                        hidden_states =
                                                          th_vit_out,
                                                        threshold =
                                                          threshold_bhmm
                                        )
 )
 RES <- threshold_Results(as(select.results, "RangedSummarizedExperiment"),
                       annotatedData = data_return)
 return(RES)


}
