globalVariables(c("betaHMMresults", "CHR", "IlmnID", "MAPINFO", "chr",
                    "cond_prob_threshold"))
#' The betaHMM model parameter estimation function
#' @description A homogeneous hidden Markov model for the beta valued DNA
#' methylation data.
#' @details The betaHMM function employs initially set parameters (utilizing a
#' basic 3-state beta hidden Markov model) to estimate the parameters of the
#' homogeneous hidden Markov model, adapted for beta-valued DNA methylation
#' data, through implementation of the Baum-Welch algorithm. Subsequently,
#' the derived parameters are utilized to ascertain the most probable sequence
#' of hidden states using the Viterbi algorithm.
#' @export
#' @param methylation_data A dataframe of dimension \eqn{(C \times (N \times R)
#' )+1} containing methylation values for \eqn{C} CpG sites from \eqn{R}
#' treatment groups each having \eqn{N} DNA samples and the IlmnID
#' for each CpG site.
#' @param annotation_file A dataframe containing the EPIC methylation
#' annotation file.
#' @param M Number of methylation states to be identified in a
#' single DNA sample.
#' @param N Number of DNA samples (patients/replicates) collected for each
#' treatment group.
#' @param R Number of treatment groups (For. eg: Benign and Tumour).
#' @param treatment_group The names of each treatment groups/
#' conditions. If no value is passed then default values of sample names,
#' e.g. Sample 1, Sample 2, etc are used as legend text (default = NULL).
#' @param parallel_process The 'TRUE' option results in parallel processing of
#' the models for each chromosome for increased computational efficiency.
#' The default option has been set as 'FALSE' due to package testing
#' limitations.
#' @param seed Seed to allow for reproducibility (default = NULL).
#' @param iterations Number of iterations for algorithm convergence
#' (default=100).
#' @param ... Extra arguments
#' @return The function returns an object of the
#' \code{\link[betaHMM:betaHMMResults]{betaHMMResults}} class
#' which contains a SimpleList of assay data containing the posterior
#' probability of each CpG site belonging to each of the
#' \eqn{K} hidden states and the following values as metadata:
#' \itemize{
#' \item K - The number of hidden states identified using the betaHMM model.
#' \item C - The number of CpG sites analysed using the betaHMM model.
#' \item N - The number of DNA samples corresponding to each
#' treatment group analysed using the betaHMM model.
#' \item R - The number of treatment groups analysed using the betaHMM model.
#' \item A - The transition matrix estimated for the betaHMM model.
#' \item tau - The initial distribution estimated for the betaHMM model.
#' \item treatment_group - The names of the treatment
#' groups/conditions analysed.
#' \item phi - The shape parameters estimated for the observed data
#' in the betaHMM model.
#' \item llk - A vector containing the log-likelihood values calculated
#' for each iteration of the algorithm.
#' \item hidden_states - The vector containing the estimated hidden states
#' for each CpG sites. }
#' @importClassesFrom S4Vectors DataFrame
#' @importMethodsFrom S4Vectors metadata
#' @importFrom stats complete.cases
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @example inst/examples/betaHMM_package.R
#'
#'
betaHMMrun <- function(methylation_data, annotation_file, M, N, R,
                        treatment_group = NULL, parallel_process = FALSE,
                        seed = NULL,iterations=100, ...) {
    if (is.null(methylation_data)) {
        stop("Please provide a dataframe containing methylation values.")
    } else if (is.null(annotation_file)) {
        stop("Please provide a dataframe containing
            annotation for EPIC array.")} else if (M > 3 | M < 2) {
        stop("M cannot be more than 3 or less than 2 as a CpG can be any of
            the 3 methylation states. ")
    } else if (R < 2) { stop("R cannot be < 2 to identify DMCs between
    2 or more conditions.") } else if (any(N < 1)) { stop("N cannot be
    less than 1 as one or more DNA replicates need to be analysed.")
    } else if (M!=round(M)|N!=round(N)|R!=round(R)){
            stop("M, N and R has to be whole numbers.")}
    if (is.null(treatment_group)) {treatment_group <-
        vapply(seq(1, R), function(x) { paste0("Sample ", x)}, character(1))
    } else if (length(treatment_group) > R) {
        stop("Treatment groups cannot be more than the number of conditions
            entered.")}
    meth_data <- methylation_data
    anno_file <- annotation_file
    final_subset <- subset(meth_data, select = -IlmnID)
    row <- nrow(final_subset);col <- ncol(final_subset)
    numeric_vec <- as.numeric(as.matrix(final_subset))
    final_subset <- matrix(numeric_vec, ncol = col)
    if (abs(max(final_subset) - 1) < 1e-06) {
        max_f <- max(final_subset[final_subset != max(final_subset)])
    } else { max_f <- max(final_subset) }
    if (abs(min(final_subset) - 0) < 1e-06) {
        min_f <- min(final_subset[final_subset != min(final_subset)])
    } else { min_f <- min(final_subset)}
    final_subset[final_subset > max_f] <- max_f
    final_subset[final_subset < min_f] <- min_f
    cols <- which(colnames(meth_data) != "IlmnID")
    meth_data[, cols] <- final_subset
    meth_data <- as.data.frame(meth_data)
    data<-merge(meth_data,anno_file[,c("IlmnID","CHR","MAPINFO")],by="IlmnID")
    col_order <- which(colnames(data) == "IlmnID" | colnames(data) == "CHR" |
                        colnames(data) == "MAPINFO")
    col_order2 <- which(colnames(data) != "IlmnID" & colnames(data) != "CHR" &
                        colnames(data) != "MAPINFO")
    data <- data[, c(col_order, col_order2)]; K <- M^R ; C <- nrow(data)
    is.scalar <- function(x) is.atomic(x) && length(x) == 1L && Im(x) == 0
    if (is.scalar(N)) {
        N <- rep(N, R)} else if (length(N) < R) {
        n <- R - length(N)
        N <- c(N, rep(N[1], n))}
    complete_data <- data
    complete_data <- complete_data[complete.cases(complete_data), ]
    sorted_data <- complete_data[order(as.numeric(complete_data$CHR),
                                        complete_data$MAPINFO), ]
    sorted_data <- as.data.frame(sorted_data)
    chr_unique <- unique(sorted_data$CHR)
    chr_unique <- as.character(sort(as.numeric(chr_unique)))
    ncores <- ifelse(parallel_process == FALSE, 2L, detectCores())
    my.cluster <- makeCluster(ncores - 1)
    registerDoParallel(cl = my.cluster)
    `%dopar%` <- foreach::`%dopar%`; `%do%` <- foreach::`%do%`
    betaHMM_workflow <- function(data_chr, K, M, N, R, chr, seed = NULL,
                                    iterations=100) {
        data <- subset(data_chr, select = -c(CHR, MAPINFO, IlmnID))
        data_w_ilmnid <- subset(data_chr, select = -c(CHR, MAPINFO))
        data <- as.data.frame(data)
        data_w_ilmnid <- as.data.frame(data_w_ilmnid)
        trained_params <- initialise_parameters(data_w_ilmnid, M, N, R, seed)
        out_baumwelch<-BaumWelch(data,trained_params,M,N,R,seed,iterations)
        out_viterbi <- Viterbi(data, M, N, R, out_baumwelch$tau,
                                out_baumwelch$A, out_baumwelch$phi)
        return(list(A = out_baumwelch$A, tau = out_baumwelch$tau,
                    sp_1 = out_baumwelch$phi$sp_1,
                    sp_2 = out_baumwelch$phi$sp_2, z = out_baumwelch$z,
                llk=out_baumwelch$log_vec,hiddenStates=out_viterbi,chr=chr))}
    comb <- function(x, ...) {lapply(seq_along(x),
                function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))}
    beta_chr_out <- foreach(chr = seq(1, length(chr_unique)),
                                    .combine = "comb", .multicombine = TRUE,
                                    .init = list(list(), list(), list(),
                                                list(), list(), list(),
                                                list(), list())) %dopar%
    betaHMM_workflow(sorted_data[sorted_data$CHR == chr_unique[chr], ],
                    K, M, N, R, chr_unique[chr],seed,iterations)
    stopCluster(cl = my.cluster)
    z_final_out <- as.data.frame(do.call(rbind, beta_chr_out[[5]]))
    chromosome_number <- unlist(beta_chr_out[[8]])
    chromosome_number_list <- paste("chr", chromosome_number)
    for (i in seq(1, 7)){ names(beta_chr_out[[i]]) <- chromosome_number_list }
    beta_out_phi <- list(); beta_out_phi[["sp_1"]] <- beta_chr_out[[3]]
    beta_out_phi[["sp_2"]]<-beta_chr_out[[4]];z<-as.data.frame(z_final_out);
    rownames(z)<-sorted_data$IlmnID; sorted_data<-as(sorted_data,"DataFrame")
    select.results<-SummarizedExperiment(assays = z,
                                        metadata = list(K = K, N = N, R = R,
                                                    A=beta_chr_out[[1]],
                                                    tau=beta_chr_out[[2]],
                                                    phi = beta_out_phi,
                                                    treatment_group =
                                                        treatment_group,
                                                    llk=beta_chr_out[[6]],
                                                    hidden_states =
                                                        beta_chr_out[[7]],
                                                    chromosome_number =
                                                        unlist(
                                                            beta_chr_out[[8]]
                                                            )))
    RES <- betaHMMResults(as(select.results, "RangedSummarizedExperiment"),
                            annotatedData = sorted_data)
    return(RES)}


#### DMC identification
#' @title DMC identification from estimated betaHMM model parameters
#' @description The dissimilarities between the cumulative distributions
#' calculated for each hidden state are determined through employment of the
#' area-under-curve (AUC) technique. By incorporating user-defined threshold
#' values for AUC alongside the associated uncertainties in membership within
#' that hidden state, the aim is to pinpoint the most distinctively methylated
#' states. This process facilitates the identification of CpGs that exhibit the
#' most notable differential methylation, guided by the predefined threshold
#' criteria.
#' @export
#' @param betaHMM_object A \code{\link[betaHMM:betaHMMResults]{betaHMMResults}}
#' object.
#' @param AUC_threshold The threshold for AUC metric for each chromosome.
#' @param uncertainty_threshold The threshold for uncertainty of belonging to
#' a particular hidden state, for each chromosome.
#' @param ... extra arguments
#'
#' @return The function returns an object of the
#' \code{\link[betaHMM:dmcResults]{dmcResults}} class
#' which contains a SimpleList of assay data which contains the
#' following values:
#' \itemize{
#' \item CHR Chromosome number
#' \item MAPINFO Mapinfo
#' \item IlmnID IlmnID
#' \item N*R columns containing methylation states
#' \item hidden_state The assigned hidden_state
#' \item DMC The value is 1 if the CpG is a DMC else 0.}
#' The object contains the following values as the metadata:
#' \itemize{
#' \item A list containing the AUC values for K hidden states for each
#' chromosome.
#' \item A list containing the conditional probability values for K hidden
#' states for each chromosome.
#' \item The treatment group labels.
#' \item K The number of hidden states estimated.
#' \item N The number of DNA replicates/patients for each treatment group.
#' \item R The number of treatment groups to be compared.
#' }
#' @importFrom stats complete.cases
#' @example inst/examples/betaHMM_package.R

dmc_identification_run <- function(betaHMM_object, AUC_threshold = 0.8,
                                    uncertainty_threshold = 0.2, ...) {
    object <- betaHMM_object; chr_unique <- chromosome_number(object)
    M <- 3; N <- N(object); R <- R(object); K <- K(object); A <- A(object)
    tau <- tau(object); phi <- phi(object); z <- assay(object)
    treatment_group <- treatment_group(object)
    is.scalar <- function(x) is.atomic(x) && length(x) == 1L && Im(x) == 0
    if (is.scalar(AUC_threshold)) {
        AUC_threshold <- rep(AUC_threshold, length(chr_unique))
    } else if (length(AUC_threshold) != length(chr_unique)) {
        stop("AUC threshold cannot be a vector with length less/greater
            than the no. of chromosomes.") }
    if (is.scalar(uncertainty_threshold)) {
        uncertainty_threshold <- rep(uncertainty_threshold, length(chr_unique))
    } else if (length(uncertainty_threshold) != length(chr_unique)) {
        stop("Uncertainty threshold cannot be a vector with length less
            than the no. of chromosomes.")}
    cond_prob_threshold <- 1 - uncertainty_threshold
    df <- as.data.frame(annotatedData(object)); col <- ncol(df)
    df_complete <- matrix(NA, nrow = 1, ncol = (col + 2))
    colnames(df_complete) <- c(colnames(df), "hidden_state", "DMC")
    auc_list <- list(); uncertainty_list <- list()
    for (i in seq(1, length(chr_unique))) {
        phi_chr<-list(); phi_chr[["sp_1"]]<-phi[[1]][[i]]
        phi_chr[["sp_2"]]<-phi[[2]][[i]]
        auc_df <- AUC_DM_analysis(M, N, R, K, tau[[i]], A[[i]], phi_chr)
        x <- as.data.frame(t(auc_df)); colnames(x) <- c("State", "AUC")
        x$AUC <- as.numeric(x$AUC); auc_list[[i]] <- x
        x_co <- x[x$AUC >= AUC_threshold[i], ]; diff_meth_x <- x_co$State
        df_chr <- df[df$CHR == chr_unique[i], ]
        if (nrow(x_co) != 0) {
            hs <- hidden_states(object)[[i]]
            df_chr$hidden_state <- hs
            z_mat <- z[row.names(z) %in% df_chr$IlmnID, ]
            if (length(diff_meth_x) > 1) {
                z_sum2 <- as.vector(rowSums(z_mat[, as.numeric(diff_meth_x)]))
            } else {  z_sum2 <- as.vector(z_mat[, as.numeric(diff_meth_x)]) }
            df_chr$DMC <- ifelse(z_sum2 >= cond_prob_threshold[i], 1, 0)
            uncertainty_list[[i]] <- 1 - z_sum2
        } else { df_chr$hidden_state <- 0; df_chr$DMC <- 0;
        uncertainty_list[[i]] <- 1}; df_complete<-rbind(df_complete, df_chr)}
    df_complete<-df_complete[-1,];rownames(df_complete)<-df_complete$IlmnID
    chrom_list <- paste("chr", chr_unique)
    names(auc_list) <- chrom_list; names(uncertainty_list) <- chrom_list
    select.results <-
        SummarizedExperiment(assays=df_complete,
        metadata =list(K=K,N=N,R=R,treatment_group=treatment_group,
                        AUC = auc_list, uncertainty =uncertainty_list))
    RES <- dmcResults(as(select.results, "RangedSummarizedExperiment"))
    return(RES)}

#### DMR identification
#' @title DMR identification from DMCs identified
#' @description Function to identify the DMRs from the DMCs identified in each
#' chromosome.
#' @export
#' @param dmc_identification_object a \code{dmcResults} object or
#' the assay data from the dmcResults.
#' @param DMC_count The minimal number of consecutive CpGs in a DMR.
#' @param parallel_process The 'TRUE' option results in parallel processing of
#' the DMCs from each chromosome for increased computational efficiency.
#' The default option has been set as 'FALSE' due to package testing
#' limitations.
#' @param ... extra arguments
#'
#' @return A \code{\link[betaHMM:dmrResults]{dmrResults}}
#' object containing a SimpleList of assay data containing the following data:
#' \itemize{
#' \item start_CpG - The starting CpG site IlmnID in the particular DMR
#' \item end_CpG -  The ending CpG site IlmnID in the particular DMR
#' \item DMR_size - Number of CPG sites identified in the DMR
#' \item chr_dmr - The chromosome corresponding to the CpG sites in the DMR.
#' \item map_start - MAPINFO of starting CpG site in the particular DMR
#' \item map_end -  MAPINFO of the ending CpG site in the particular DMR}
#' The object also returns the chromosomes analysed by the betaHMM model as
#' the metadata.
#' @importFrom stats complete.cases
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @example inst/examples/betaHMM_package.R

dmr_identification_run <- function(dmc_identification_object, DMC_count = 2,
                                    parallel_process = FALSE, ...) {
    dmc_df <- dmc_identification_object; chr_unique <- unique(dmc_df$CHR)
    ncores <- ifelse(parallel_process == FALSE, 2L, detectCores())
    my.cluster <- makeCluster(ncores - 1)
    registerDoParallel(cl = my.cluster)
    `%dopar%` <- foreach::`%dopar%`;`%do%` <- foreach::`%do%`
    dmr_parallel <- function(dmc_df_chr, DMC_count) {
        block_counter<-0; block_start <- 0; block_end <- 0; block_length <- 0
        C <- nrow(dmc_df_chr); mat <- matrix(NA, C, 3)
        for (i in seq(1, C)) { if (dmc_df_chr[i, "DMC"] == 1) {
                if (block_length == 0) {
                    block_start <- i ;block_length <- 1 } else {
                    block_length <- block_length + 1}
                if (block_length >= DMC_count) {
                    block_counter <- block_counter + 1 }
            } else { if (block_length >= DMC_count) {
                mat[block_counter, ] <- c(block_start, i - 1, block_length)}
                block_start <- 0; block_length <- 0}}
        mat<-mat[complete.cases(mat),];mat<-as.data.frame(mat);df_unique<-mat
        colnames(df_unique) <- c("start_CpG", "end_CpG", "DMR_size")
        cpg_names <- as.vector(dmc_df_chr[, "IlmnID"])
        start_cpg <- vapply(df_unique$start_CpG, function(x) {
            cpg_names[x]}, character(1))
        end_cpg<-vapply(df_unique$end_CpG,
                        function(x) {cpg_names[x]},character(1))
        df_unique$start_CpG <- start_cpg; df_unique$end_CpG <- end_cpg
        map_start <- vapply(df_unique$start_CpG, function(x) {
        dmc_df_chr[dmc_df_chr$IlmnID == x, "MAPINFO"]}, numeric(1))
        map_end <- vapply(df_unique$end_CpG, function(x) {
        dmc_df_chr[dmc_df_chr$IlmnID == x, "MAPINFO"]}, numeric(1))
        df_unique$CHR <- unique(dmc_df_chr$CHR)
        df_unique$map_start <- map_start; df_unique$map_end <- map_end
        for (j in seq(1, nrow(df_unique))) {
            start <- which(dmc_df_chr$MAPINFO == df_unique[j, 5] &
                            dmc_df_chr$CHR == df_unique[j, 4])
            end <- which(dmc_df_chr$MAPINFO == df_unique[j, 6] &
                            dmc_df_chr$CHR == df_unique[j, 4])
            x <- dmc_df_chr[start:end, "IlmnID"]
            df_unique[j,"DMCs"]<-paste(x,collapse=",")}; return(df_unique)}
    dmr_mat <- foreach(chr = seq(1, length(chr_unique)),
                                .combine = rbind) %dopar%
        dmr_parallel(dmc_df[dmc_df$CHR == chr_unique[chr],], DMC_count)
    stopCluster(cl = my.cluster)
    dmr_df <- dmr_mat[complete.cases(dmr_mat), ]
    rownames(dmr_df) <- dmr_df$start_CpG
    select.results <- SummarizedExperiment(assays = dmr_df,
                        metadata = list(chromosome_number =chr_unique))
    RES <- dmrResults(as(select.results, "RangedSummarizedExperiment"))
    return(RES)}

