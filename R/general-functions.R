globalVariables(c("betaHMMresults","CHR","IlmnID","MAPINFO"))
#' The betaHMM wrapper function
#' @description A homogeneous hidden Markov model for the beta valued DNA
#'              methylation data.
#' @details A novel approach utilizing a homogeneous hidden Markov model and
#'          effectively model untransformed beta values to identify DMCs while
#'          considering the spatial correlation of the adjacent CpG sites
#' @export
#' @param data A dataframe of dimension \eqn{C \times (N \times R)} containing
#'             methylation values for \eqn{C} CpG sites from \eqn{R}
#'             treatment groups each having \eqn{N} replicates or each collected
#'             from \eqn{N} patients.
#' @param M Number of methylation states to be identified in a DNA sample.
#' @param N Number of patients or DNA sample replicates collected for each
#'          treatment group.
#' @param R Number of treatment groups (For. eg: Benign and Tumour).
#' @param seed Seed to allow for reproducibility (default = NULL).
#' @param ... Extra arguments
#' @return The function returns an object of the
#'         \code{\link[betaHMM:betaHMM]{betaHMM}} class
#'         which contains the following values:
#' \itemize{
#' \item function_call - The parameters passed as arguments to the function
#'                       \code{\link[betaclust:betaclust]{betaclust}}.
#' \item K - The number of hidden states identified using the BHMMs.
#' \item C - The number of CpG sites analysed using the BHMMs.
#' \item N - The number of patients or DNA replicates corresponding to each
#'           treatment group analysed using the BHMMs.
#' \item R - The number of treatment groups analysed using the BHMMs.
#' \item A - The transition matrix for the BHMM model.
#' \item tau - The initial distribution for the BHMM model.
#' \item phi - The shape parameters for the observation sequence data
#'                in the BHMM model.
#' \item llk - A vector containing the log-likelihood values calculated
#'                    for each iteration of the algorithm.
#' \item z - A matrix of dimension \eqn{C \times K} containing the posterior
#'              probability of each CpG site belonging to each of the
#'              \eqn{K} clusters.
#' \item hidden_states - The vector containing the estimated hidden states
#'                          for each CpG sites.
#'    }
#' @importClassesFrom S4Vectors DataFrame
#' @importMethodsFrom S4Vectors metadata
#'
betaHMMrun<-function(data,M,N,R,seed=NULL,...)
{
  #call_function<-match.call()
  K=M^R
  C=nrow(data)
  complete_data=data
  sorted_data<-complete_data[order(complete_data$CHR,complete_data$MAPINFO),]
  sorted_data=as.data.frame(sorted_data)

  # Remove using subset
  data <- subset(sorted_data, select = -c(CHR, MAPINFO, IlmnID))

  data=as.data.frame(data)
  ## Getting initialized parameters
  trained_params=initialise_parameters(data,M,N,R,seed)
  #print("initialised")
  ## Baum-Welch algorithm for estimating BHMM parameters
  out_baumwelch = BaumWelch(data,trained_params,M,N,R,seed=seed)
  #print("BW done")
  out_viterbi=Viterbi(data,M,N,R,out_baumwelch$tau,out_baumwelch$A,
                      out_baumwelch$phi)

  # betaHMM_out<-list(function_call=call_function,K=K,C=C,N=N,R=R,
  #                   A = out_baumwelch$A,
  #                   tau = out_baumwelch$tau,
  #                   phi=out_baumwelch$phi,
  #                   log_vec=out_baumwelch$log_vec,
  #                   z=out_baumwelch$z,
  #                   hidden_states=out_viterbi)


  ####################################
  ## RETURN RESULTS
  ####################################
  z=as.data.frame(out_baumwelch$z)
  rownames(z)<-sorted_data$IlmnID
  sorted_data=as(sorted_data,"DataFrame")
  select.results <- SummarizedExperiment(assays=z,
                                         metadata=list(K=K,N=N,R=R,
                                                       A = out_baumwelch$A,
                                                       tau = out_baumwelch$tau,
                                                       phi=out_baumwelch$phi,
                                                       llk=out_baumwelch$log_vec,
                                                       hidden_states=out_viterbi))
  RES <- betaHMMResults(as(select.results, "RangedSummarizedExperiment"),annotatedData=sorted_data)
  return(RES)

}

#### DMR identification
#' @title DMR identification from DMCs
#' @description Function to identify the DMRs for the DMCs identified using BHMM
#' @details Function to identify the DMRs for the DMCs identified using BHMM
#' @export
#' @param object A \code{\link[betaHMM:betaHMM]{betaHMM}} object.
#' @param diff_meth_cluster The clusters identified as differentially methylated
#' @param ... extra arguments
#'
#' @return A dataframe containing the following columns:
#' \itemize{
#' \item start_CpG - The starting CpG site in the particular DMR
#' \item end_CpG -  The ending CpG site in the particular DMR
#' \item DMR_size - Number of CPG sites identified in the DMR
#' \item chr_dmr - The chromosome corresponding to the CpG sites in the DMR.
#' \item map_start - MAPINFO of starting CpG site in the particular DMR
#' \item map_end -  MAPINFO of the ending CpG site in the particular DMR}
#' @importFrom stats complete.cases
dmr_identification_run<-function(object,diff_meth_cluster,...)
{

  x=object
  data=as.data.frame(annotatedData(object))
  C=nrow(data)
  df_dmr<-cbind(data,as.vector(hidden_states(object)))
  colnames(df_dmr)[ncol(df_dmr)]="hidden_states"
  df_dmr$true_dmc<-ifelse(df_dmr$hidden_states %in% diff_meth_cluster , 1,0)
  block_counter <- 0
  block_start <- 0
  block_end<-0
  mat<-matrix(NA,C,3)
  block_length <- 0

  # Loop through the sequence of numbers
  for (i in 1:C) {

    # Check if the current number is equal to 3
    if (df_dmr[i,"true_dmc"] == 1) {

      # If this is the start of a new block, save the starting index and
      #set the block length to 1
      if (block_length == 0) {
        block_start <- i
        block_length <- 1

        # If this is part of an existing block, increment the block length
      } else {
        if(df_dmr[i,"CHR"]==df_dmr[i-1,"CHR"])
        {block_length <- block_length + 1}
        else{

          if(block_length>=2)
          {
            mat[block_counter,]<-c(block_start,i-1,block_length)
          }
          block_start <- 0
          block_length <- 0
        }


      }

      # If the block length is greater than 2,
      #increment the block counter and print the
      #starting and ending index of the block
      if (block_length >= 2) {
        block_counter <- block_counter + 1
        #cat("Block", block_counter, "starts at index", block_start,
        #"and ends at index", i, "\n")
        #mat[block_counter,]<-c(block_counter,block_start,i)
      }

      # If this is not part of a block, reset the block_start and
      #block_length variables
    } else {
      if(block_length>=2)
      {
        mat[block_counter,]<-c(block_start,i-1,block_length)
      }
      block_start <- 0
      block_length <- 0
    }
  }
  mat<-mat[stats::complete.cases(mat),]
  mat<-as.data.frame(mat)


  df_unique<-mat
  colnames(df_unique)<-c("start_CpG","end_CpG","DMR_size")
  df_unique2=df_unique
  cpg_names<-as.vector(data[,"IlmnID"])
  start_cpg<-sapply(df_unique2$start_CpG,function(x){cpg_names[x] })
  end_cpg<-sapply(df_unique2$end_CpG,function(x){cpg_names[x] })
  df_unique$start_CpG<-start_cpg
  df_unique$end_CpG<-end_cpg

  CHR<-sapply(df_unique$start_CpG,function(x)
  {data[data$IlmnID==x,"CHR"]})
  map_start<-sapply(df_unique$start_CpG,function(x)
  {data[data$IlmnID==x,"MAPINFO"]})
  map_end<-sapply(df_unique$end_CpG,function(x)
  {data[data$IlmnID==x,"MAPINFO"]})
  df_unique$CHR<-CHR
  df_unique$map_start<-map_start
  df_unique$map_end<-map_end
  #
  # return(df_unique)
  ####################################
  ## RETURN RESULTS
  ####################################
  dmr_df=df_unique
  rownames(dmr_df)=dmr_df$start_CpG
  select.results <- SummarizedExperiment(assays=dmr_df)
  RES <- dmrResults(as(select.results, "RangedSummarizedExperiment"))
  return(RES)
}

