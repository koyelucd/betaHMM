#' MethylationEPIC manifest data.
#'
#' A dataset containing a subset of the manifest data from the Illumina
#' MethylationEPIC beadchip array. A subset of the complete dataset has been
#' uploaded in the package for testing purpose. The complete dataset is
#' available on \href{https://github.com/koyelucd/betaclust}{GitHub}.
#'
#' @seealso \code{\link{pca_methylation_data}}
#' @format A data frame with 100 rows and 9 columns.
#' \itemize{
#'         \item{IlmnID: The unique identifier from the Illumina CG database,
#'         i.e. the probe ID.}
#'         \item{Genome_Build: The genome build referenced by the Infinium
#'         MethylationEPIC manifest.}
#'         \item{CHR: The chromosome containing the CpG (Genome_Build = 37).}
#'         \item{MAPINFO: The chromosomal coordinates of the CpG sites.}
#'         \item{UCSC_RefGene_Name: The target gene name(s), from the UCSC
#'         database. Note: multiple listings of the same gene name indicate
#'         splice variants.}
#'         \item{UCSC_RefGene_Accession: The UCSC accession numbers of the
#'         target transcripts. Accession numbers are in the same order as the
#'         target gene transcripts.}
#'         \item{UCSC_RefGene_Group: Gene region feature category describing
#'         the CpG position, from UCSC. Features are listed in the same order
#'         as the target gene transcripts.}
#'         \item{UCSC_CpG_Islands_Name: The chromosomal coordinates of the CpG
#'         Island from UCSC.}
#'         \item{Relation_to_UCSC_CpG_Island: The location of the CpG relative
#'         to the CpG island.}
#'         }
#' @usage data(annotation_data)
#' @return A data frame containing the array design for Illumina’s
#' Human Methylation EPIC microarray for the chromosome 7. Based on the v1.0b2
#' version of the manifest file.
"annotation_data"
