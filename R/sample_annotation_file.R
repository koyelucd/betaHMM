#' MethylationEPIC manifest data.
#'
#' A dataset containing a subset of the manifest data from the
#' Illumina MethylationEPIC beadchip array. A subset of the complete dataset
#' has been uploaded in the package for testing purpose. The complete dataset
#' is available on \href{https://github.com/koyelucd/betaHMM}{GitHub}.
#'
#' @seealso \code{\link{sample_methylation_file}}
#' @format A data frame with 100 rows and 9 columns.
#' \itemize{
#' \item{IlmnID: The unique identifier from the Illumina CG database,
#' i.e. the probe ID.}
#' \item{Genome_Build: The genome build referenced by the Infinium
#' MethylationEPIC manifest.}
#' \item{CHR: The chromosome containing the CpG (Genome_Build = 37).}
#' \item{MAPINFO: The chromosomal coordinates of the CpG sites.}
#' \item{UCSC_RefGene_Name: The target gene name(s), from the UCSC database.
#' Note: multiple listings of the same gene name indicate splice variants.}
#' \item{UCSC_CpG_Islands_Name: The chromosomal coordinates of the CpG
#' Island from UCSC.}}
#' @usage data(sample_annotation_file)
#' @return A data frame containing annotations for simulated data.
"sample_annotation_file"