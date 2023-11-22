#' @title HMM parameter estimation using Baum-Welch algorithm
#' @description  The function determines the parameters of a homogeneous
#' beta hidden Markov model (betaHMM), wherein the Baum-Welch algorithm
#' constitutes a variant of the EM (Estimation-Maximization) procedure.
#' @keywords internal
#' @param data A dataframe of dimension \eqn{C \times (N \times R)} containing
#' methylation values for \eqn{C} CpG sites from \eqn{R}
#' treatment groups each having \eqn{N} DNA samples.
#' @param trained_params A list containing the initialized model parameters:
#' \itemize{
#' \item A - The transition matrix for the betaHMM model.
#' \item tau - The initial distribution for the betaHMM model.
#' \item phi - The shape parameters for the observation sequence data
#' in the betaHMM model.}
#' @param K The number of hidden states identified using the betaHMM model.
#' @param M Number of methylation states to be identified in a single
#' DNA sample.
#' @param N Number of DNA samples (patients/replicates) collected for each
#' treatment group.
#' @param R Number of treatment groups (For. eg: Benign and Tumour).
#' @param seed Seed to allow for reproducibility (default = NULL).
#' @param iterations Number of iterations for algorithm convergence.
#' @return A list containing:
#' \itemize{
#' \item A - The transition matrix estimated for the betaHMM model.
#' \item tau - The initial distribution estimated for the betaHMM model.
#' \item phi - The shape parameters estimated for the observed data
#' in the betaHMM model.
#' \item log_vec - A vector containing the log-likelihood values
#' calculated for each iteration of the algorithm.
#' \item z - A matrix of dimension \eqn{C \times K} containing the
#' conditional posterior probability of each CpG site belonging to each of
#' the \eqn{K} hidden states.
#' }
#' @importFrom stats dbeta
#'
#'

BaumWelch <- function(data, trained_params=NULL,K, M, N, R,
                      seed = NULL,iterations= 100){
    C <- nrow(data)
    if(K==M)
    {
        trained_params <- initialize_parameters_th(data, M, N, R, seed)
    }
    oldlogL <- 1
    BW_limit_accuracy <- 1e-5
    logL_vec <- vector()
    logL_vec <- oldlogL
    iter<-1
    while(iter<=iterations){
        iter<-iter+1
        probabilities <- list()
        n1 <- 1
        n2 <- N[1]
        for(r in seq(1,R))
        {
            if(K==M)
            {
                sp1<-trained_params$phi$sp_1
                sp2<-trained_params$phi$sp_2
            }else{
                sp1<-trained_params$phi$sp_1[r, ]
                sp2<-trained_params$phi$sp_2[r, ]
            }
            for(n in n1:n2)
            {
                probabilities[[n]] <- matrix(data[,n],ncol=K,nrow=nrow(data))
                probabilities[[n]] <- t(apply(X=probabilities[[n]],MARGIN =1,
                                              FUN = dbeta,
                                              shape1 = sp1,
                                              shape2=sp2))}
            if((r+1) <= R){
                n1 <- n1 + N[r]
                n2 <- n2 + N[r+1] } }
        prob <- matrix(probabilities[[1]], nrow = C, ncol = K)
        for(k in 2:ncol(data)){
            temp <- as.matrix(probabilities[[k]])
            prob <- prob*temp}
        probabilities<-matrix(0, ncol = K, nrow = C)
        probabilities<-prob
        forward_alpha <- forward(probabilities, trained_params)
        backward_beta <- backward(probabilities, trained_params)
        log_alpha <- forward_alpha$log_alpha
        logL_alpha_T <- forward_alpha$scaled_logL
        log_beta <- backward_beta$log_beta
        logL_beta_1 <- backward_beta$scaled_logL
        logL_alpha_t_beta_t <- 0
        middle_t <- round(C / 2)
        t <- max(log_alpha[middle_t,])
        logL_alpha_t_beta_t <- sum(exp(log_alpha[middle_t, seq_len(K)] +
                                           log_beta[middle_t, seq_len(K)]))
        logL_alpha_t_beta_t <- log(logL_alpha_t_beta_t)
        logL <- logL_alpha_t_beta_t
        if (logL == -Inf | logL == Inf | is.na(logL)){logL <- logL_alpha_T}
        if (logL == -Inf | logL == Inf | is.na(logL)) {logL <- logL_beta_1}
        eta <- matrix(c(0), ncol = K, nrow = C)
        eta<-exp(log_alpha+log_beta-logL)
        alpha_c <- log_alpha[seq(1,C-1),, drop=FALSE]
        beta_c  <- log_beta[2:C,, drop=FALSE]
        prob_log <- log(probabilities[2:C,, drop=FALSE])
        combinations <- expand.grid(i = seq(1, K), j = seq(1, K))
        xi <- apply(combinations, 1, function(comb) {
            i <- comb[["i"]]
            j <- comb[["j"]]
            exp(alpha_c[,i] + log(trained_params$A[i, j]) +
                    prob_log[, j] + beta_c[, j] - logL)
        })
        xi <- array(unlist(xi), dim = c((C-1), K, K))
        tau <- eta[1,]
        A <- matrix(0,K,K)
        A1 <- apply(xi, c(2,3), sum)
        A2 <- apply(A1, FUN = sum, MARGIN = 1)
        A <- A1/A2
        y1 <- matrix(0, nrow = R, ncol = K)
        y2 <- matrix(0, nrow = R, ncol = K)
        sh_p1 <- matrix(0, nrow = R, ncol = K)
        sh_p2 <- matrix(0, nrow = R, ncol = K)
        term <- matrix(0, nrow = R, ncol = K)
        process_r<-function(r,data,N,n1,n2){
            start<-n1[r]
            end<-n2[r]
            if(N[r] == 1)
            {
                y1 <- apply((eta*(log(data[,start:end, drop=FALSE]))),
                                 FUN=sum,MARGIN=2)/
                    (N[r]*apply(eta,FUN=sum,MARGIN=2))
                y2 <- apply((eta*(log(1-data[,start:end, drop=FALSE]))),
                                 FUN=sum,MARGIN=2)/
                    (N[r]*apply(eta,FUN=sum,MARGIN=2))
            }else{
                y1 <- apply((eta*rowSums(log(data[,start:end, drop=FALSE]))),
                                 FUN=sum,MARGIN=2)/
                    (N[r]*apply(eta,FUN=sum,MARGIN=2))
                y2 <- apply((eta*rowSums(log(1-data[,start:end, drop=FALSE]))),
                                 FUN=sum,MARGIN=2)/
                    (N[r]*apply(eta,FUN=sum,MARGIN=2))}
            term <- ((exp(-y1)-1) *(exp(-y2)-1)) - 1
            shp1 <- 0.5 + (0.5 * exp(-y2) /term)
            shp2 <- (0.5 * exp(-y2) *(exp(-y1) - 1)) /term
            return(list(shp1=shp1,shp2=shp2))

        }
        indices <- seq(1, R)
        n1<-vector()
        n2<-vector()
        n1[1]<-1
        n2[1]<-N[1]
        n1 <- c(1, cumsum(N[-R]+1))
        n2 <- cumsum(N)
        result<-lapply(indices, process_r, data = data, N = N,n1,n2)
        shp1_matrix <- matrix(NA, nrow = R, ncol = K)
        shp2_matrix <- matrix(NA, nrow = R, ncol = K)
        shp1_matrix <- do.call(rbind, lapply(result, function(res) res$shp1))
        shp2_matrix <- do.call(rbind, lapply(result, function(res) res$shp2))
        phi<-list(sp_1=shp1_matrix,sp_2=shp2_matrix)
        trained_params$A <- A
        trained_params$tau <- tau
        trained_params$phi <- phi
        logL_vec<-c(logL_vec,logL)
        diff_logL  <-  (logL - oldlogL) / oldlogL
        if (diff_logL > 0 &diff_logL < BW_limit_accuracy)
        {
            reached_limit_of_accuracy  <-  TRUE
            break}
         oldlogL <- logL
         }
    log_vec <- logL_vec[-1]
    return(list(A = A,  tau = tau, phi = phi,log_vec = log_vec, z = eta))
}
