#' @title Threshold indentification function
#' @description HMM for beta valued DNA data for a single treatment condition
#'
#' @export
#' @param data Methylation data and IlmnID. Maybe provided as
#' a matrix or dataframe.
#' @param package_workflow Flag set to TRUE if method called from package
#' workflow.
#' If set to FALSE then the parameter annotation_file needs to be supplied to
#' the function.
#' @param annotation_file A dataframe containing the EPIC methylation
#' annotation file.
#' @param M Number of methylation states to be identified in a
#' single DNA sample.
#' @param N Number of DNA samples (patients/replicates) collected for each
#' treatment group.
#' @param parameter_estimation_only If only model parameters are to be
#' estimated then value is TRUE else FALSE.
#' @param seed Seed to allow for reproducibility (default = NULL).
#' @param ... Extra parameters.
#' @return  An S4 object of class \code{threshold_Results}, where conditional
#' probabilities
#' of each CpG site belonging to a one of the \eqn{M} methylation states
#' is stored as a SimpleList of
#' assay data, and the corresponding estimated model parameters, the thresholds
#' and most probable hidden state sequence for each chromosome are
#' stored as metadata.
#'
#' @importFrom stats kmeans
#' @importFrom stats sd
#' @examples
#' ## Use simulated data for the betaHMM workflow example
#' set.seed(12345)
#' library(betaHMM)
#'
#' ## read files
#' data(sample_methylation_file)
#' head(sample_methylation_file)
#' data(sample_annotation_file)
#' head(sample_annotation_file)
#' ##merge data
#' df=merge(sample_annotation_file[,c('IlmnID','CHR','MAPINFO')],
#' sample_methylation_file,by='IlmnID')
#'
#' ## sort data
#' df=df[order(df$CHR,df$MAPINFO),]
#' thr_out=threshold_identification(df[,c(1,4:7)],package_workflow=TRUE,M=3,4,
#' parameter_estimation_only=TRUE,seed=12345)
#'



threshold_identification_run<-function(data,package_workflow=TRUE,
                                        annotation_file=NULL,M,N,
                                        parameter_estimation_only=FALSE,
                                        seed = NULL,...)
    {
    if (package_workflow == FALSE) {if (is.null(annotation_file))
            stop("Annotation file cannot be empty if data is not sorted.")
        final_subset <- subset(data, select = -IlmnID)
        if (abs(max(final_subset) - 1) < 1e-06) {
            max_f <- max(final_subset[final_subset != max(final_subset)])
        } else {max_f <- max(final_subset)}
        if (abs(min(final_subset) - 0) < 1e-06) {
            min_f <- min(final_subset[final_subset != min(final_subset)])
        } else {min_f <- min(final_subset)}
        final_subset[final_subset > max_f] <- max_f
        final_subset[final_subset < min_f] <- min_f
        cols<-which(colnames(data)!="IlmnID"); data[, cols] <- final_subset
        data_merged<-merge(data,annotation_file[,c("IlmnID","CHR",
                                                    "MAPINFO")],by="IlmnID")
        col_order <- which(colnames(data_merged) == "IlmnID" |
                                colnames(data_merged) == "CHR" |
                                colnames(data_merged) == "MAPINFO")
        col_order2 <- which(colnames(data_merged) != "IlmnID" &
                                colnames(data_merged) != "CHR" &
                                colnames(data_merged) != "MAPINFO")
        data_merged<-data_merged[,c(col_order,col_order2)]
        complete_data<-data_merged
        complete_data <- complete_data[stats::complete.cases(complete_data), ]
        sorted_data<-complete_data[order(as.numeric(complete_data$CHR),
                                        as.numeric(complete_data$MAPINFO)), ]
        sorted_data <- as.data.frame(sorted_data)
        data_final <- subset(sorted_data, select = -c(CHR, MAPINFO, IlmnID))
        data_return <- subset(sorted_data, select = -c(CHR, MAPINFO))
    } else {data_final<-subset(data,select=-c(IlmnID)); data_return<-data}
    th_bw_out <- BaumWelch_th(data_final, M = M, N = N, R = 1, seed = seed)
    if (parameter_estimation_only == FALSE) {
        th_vit_out<-Viterbi_th(data_final,M=M,N=N,R=1,th_bw_out$tau,
                                th_bw_out$A, th_bw_out$phi)
        tau <- as.numeric(table(th_vit_out)/(nrow(data_final)))
        threshold_bhmm <- threshold_values(data_final, tau, th_bw_out$phi)
    } else {th_vit_out <- NA; threshold_bhmm <- NA }
    z <- as.data.frame(th_bw_out$z); rownames(z) <- data_return[, "IlmnID"]
    data_return <- as(data_return, "DataFrame")
    select.results <- SummarizedExperiment(assays = z, metadata =
                                            list(model_parameters = th_bw_out,
                                            K=M,hidden_states=th_vit_out,
                                            threshold = threshold_bhmm))
    RES <- threshold_Results(as(select.results, "RangedSummarizedExperiment"),
                                annotatedData = data_return)
    return(RES)}
