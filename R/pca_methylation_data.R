#' @title Simulated DNA methylation data
#'
#' @description A subset of the dataset containing beta methylation values from
#' \eqn{R=2} sample types (Benign and Tumour), collected from \eqn{N=4}
#' patients from the a prostate cancer study. The dataset contains methylation
#' values corresponding to chromosome 7.
#'
#' @details The array data were then normalized and and probes located outside
#' of CpG sites and on the sex chromosome were filtered out. The CpG sites with
#' missing values were removed from the resulting dataset. A subset of the
#' complete dataset has been uploaded in the package for testing purposes.
#' The complete dataset is available on
#' \href{https://github.com/koyelucd/betaclust}{GitHub}.
#' @seealso \code{\link{annotation_data}}
#' @format A data frame with 38672 rows and 9 columns. The data contain no
#' missing values.
#' \itemize{
#' \item{IlmnID: The unique identifier from the Illumina CG database,
#' i.e. the probe ID.}
#' \item{Benign_Patient_1: Methylation values from benign tissue from
#' patient 1.}
#' \item{Benign_Patient_2: Methylation values from benign tissue from
#' patient 2.}
#' \item{Benign_Patient_3: Methylation values from benign tissue from
#' patient 3.}
#' \item{Benign_Patient_4: Methylation values from benign tissue from
#' patient 4.}
#' \item{Tumour_Patient_1: Methylation values from tumor tissue from
#' patient 1.}
#' \item{Tumour_Patient_2: Methylation values from tumor tissue from
#'  patient 2.}
#'  \item{Tumour_Patient_3: Methylation values from tumor tissue from
#'  patient 3.}
#'  \item{Tumour_Patient_4: Methylation values from tumor tissue from
#'  patient 4.}}
#' @usage data(pca_methylation_data)
#' @return A data frame containing a subset of methylation data
#' from real study.
"pca_methylation_data"
