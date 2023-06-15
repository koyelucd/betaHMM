#' @title Estimation Using the Baum-Welch Algorithm
#' @description The Baum-Welch algorithm is used to estimate the parameters of
#'              the beta distributed hidden Markov model (BHMM).
#' @details The Baum-Welch algorithm is used to estimate the parameters of
#'              the beta distributed hidden Markov model (BHMM).
#' @export
#' @param data A dataframe of dimension \eqn{C \times (N \times R)} containing
#'             methylation values for \eqn{C} CpG sites from \eqn{R}
#'             treatment groups each having \eqn{N} replicates or each collected
#'             from \eqn{N} patients.
#' @param trained_params A list containing the initialised model parameters:
#' \itemize{
#'    \item A - The transition matrix for the BHMM model.
#'    \item tau - The initial distribution for the BHMM model.
#'    \item phi - The shape parameters for the observation sequence data
#'                in the BHMM model.
#'    }
#' @param M Number of methylation states to be identified in a DNA sample.
#' @param N Number of patients or DNA sample replicates collected for each
#'          treatment group.
#' @param R Number of treatment groups (For. eg: Benign and Tumour).
#' @param n.iter Number of iterations for algorithm convergence.
#' @param seed Seed to allow for reproducibility (default = NULL).
#' @return A list containing:
#' \itemize{
#'    \item A - The transition matrix for the BHMM model.
#'    \item tau - The initial distribution for the BHMM model.
#'    \item phi - The shape parameters for the observation sequence data
#'                in the BHMM model.
#'    \item log_vec - A vector containing the log-likelihood values calculated
#'                    for each iteration of the algorithm.
#'    \item z - A matrix of dimension \eqn{C \times K} containing the posterior
#'              probability of each CpG site belonging to each of the
#'              \eqn{K} clusters.
#'    }
#' @importFrom stats dbeta
#'
#'

BaumWelch = function(data,trained_params,M,N,R,n.iter=100,seed=NULL){

  K=as.numeric(M)^as.numeric(R)
  C = nrow(data)
  #print(str(data))
  oldlogL <- 1
  BW_limit_accuracy = 1e-5
  logL_vec=vector()
  logL_vec=oldlogL

  for(i in 1:n.iter){

    probabilities <-  list()
    n1=1
    n2=N

    for(r in 1:R)
    {
      for(n in n1:n2)
      {
        probabilities[[n]]=matrix(data[,n], ncol = K, nrow = nrow(data))
        probabilities[[n]]=t(apply(X = probabilities[[n]], MARGIN = 1,
                                   FUN = stats::dbeta,
                                   shape1 = trained_params$phi$sp_1[r,],
                                   shape2 = trained_params$phi$sp_2[r,]))

      }
      n1=n1+N
      n2=n2+N
    }
    prob=matrix(probabilities[[1]],nrow=C,ncol=9)
    for(k in 2:(N*R))
    {
      temp=as.matrix(probabilities[[k]])
      prob=prob*temp
    }
    probabilities<-matrix(0, ncol = K, nrow = C)
    probabilities<-prob

    #print("1. probabilities calculated")

    ## Generating alpha and beta matrices
    forward_alpha = forward(probabilities, trained_params)
    backward_beta = backward(probabilities, trained_params)

    ## Calculating log-likelihood
    log_alpha=forward_alpha$log_alpha
    logL_alpha_T=forward_alpha$scaled_logL

    log_beta=backward_beta$log_beta
    logL_beta_1=backward_beta$scaled_logL

    logL_alpha_t_beta_t <- 0
    middle_t <- round(C / 2)
    t <- max(log_alpha[middle_t,])
    for (i in 1:K)
    {
      logL_alpha_t_beta_t <- logL_alpha_t_beta_t +
                             exp( log_alpha[middle_t,i] +
                             log_beta[middle_t,i])
    }
    logL_alpha_t_beta_t <- log(logL_alpha_t_beta_t)

    logL <- logL_alpha_t_beta_t
    if (logL == -Inf | logL == Inf | is.na(logL))
    {
      logL <- logL_alpha_T
    }
    if (logL == -Inf | logL == Inf | is.na(logL))
    {
      logL <- logL_beta_1
    }

    #print("LogL calculated by forward backward")
    ## Calculating the joint and conditional probabilities

    eta=matrix(c(0),ncol=K,nrow=C)
    eta=exp(log_alpha+log_beta-logL)

    #print("eta calculated")

    xi = array(NA, dim=c( (C-1),K, K))

    alpha_c=log_alpha[1:C-1,]
    beta_c=log_beta[2:C,]
    prob_log=log(probabilities[2:C,])
    for (j in 1:K)
    {
      for (i in 1:K)
      {

        xi[,i,j] <- exp(alpha_c[,i] +
                        log(trained_params$A[i,j]) +
                        prob_log[,j]+
                        beta_c[,j]-
                        logL)
      }
    }

    #print("xi calcultaed")
    ## Initial distribution proportions
    tau=eta[1,]



    ## Estimating transmission probabilities
    A=matrix(0,K,K)

    A1=apply(xi,c(2,3),sum)
    A2=apply(A1,FUN = sum,MARGIN = 1)
    A=A1/A2

    #print("A calculated")

    ## Estimating the shape parameters using the digamma approximation approach

    y1=matrix(0,nrow=R,ncol=K)
    y2=matrix(0,nrow=R,ncol=K)
    sh_p1=matrix(0,nrow=R,ncol=K)
    sh_p2=matrix(0,nrow=R,ncol=K)
    term=matrix(0,nrow=R,ncol=K)
    for(r in 1:R) ## treatment groups
    {
      if(N==1) ## if no replicates or analysing a single patient
      {
        y1[r,]=apply((eta*(log(data[,(((r-1)*N)+1):(((r-1)*N)+N)]))),
                     FUN=sum,MARGIN = 2)/(N*apply(eta,FUN=sum,MARGIN = 2))
        y2[r,]=apply((eta*(log(1-data[,(((r-1)*N)+1):(((r-1)*N)+N)]))),
                     FUN=sum,MARGIN = 2)/(N*apply(eta,FUN=sum,MARGIN = 2))
      }else{
        y1[r,]=apply((eta*rowSums(log(data[,(((r-1)*N)+1):(((r-1)*N)+N)]))),
                     FUN=sum,MARGIN = 2)/(N*apply(eta,FUN=sum,MARGIN = 2))
        y2[r,]=apply((eta*rowSums(log(1-data[,(((r-1)*N)+1):(((r-1)*N)+N)]))),
                     FUN=sum,MARGIN = 2)/(N*apply(eta,FUN=sum,MARGIN = 2))
      }
      term[r,]=((exp(-y1[r,])-1)*(exp(-y2[r,])-1))-1
      sh_p1[r,]=0.5+(0.5*exp(-y2[r,])/term[r,])
      sh_p2[r,]= (0.5*exp(-y2[r,])*(exp(-y1[r,])-1))/term[r,]
    }



    phi <- list(sp_1 = sh_p1, sp_2 = sh_p2)

    #print("phi calculated")


    trained_params$A=A
    trained_params$tau=tau
    trained_params$phi=phi


    ## Comparing log-likelihood calculated with the previous log-likelihood
    ## Exiting recursion if the difference is below the selected epsilon
    logL_vec<-c(logL_vec,logL)
    difference_old_logL_and_new_logL = (logL - oldlogL)/oldlogL

    if (difference_old_logL_and_new_logL > 0 &
        difference_old_logL_and_new_logL < BW_limit_accuracy)
    {
      reached_limit_of_accuracy = TRUE
      break
    }

    oldlogL <- logL
  }
  log_vec=logL_vec[-1]
  return(list(A = A,  tau = tau,phi=phi,log_vec=logL_vec,z=eta))
}
