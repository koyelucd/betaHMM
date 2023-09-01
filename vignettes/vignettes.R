## ----setup, include = FALSE---------------------------------------------------
knitr::include_graphics("betaHMM_hex.png")
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(eval = TRUE)

#knitr::opts_chunk$set(dev = 'png')
#knitr::opts_chunk$set(dpi=100)

## ----package, include=TRUE, echo=TRUE, message=FALSE,warning=FALSE------------
library(betaHMM)

## ----data,include=TRUE, echo=TRUE---------------------------------------------
data(pca_methylation_data)
head(pca_methylation_data)
data(annotation_data)
head(annotation_data)



## ----betaHMM,include=TRUE, echo=TRUE------------------------------------------
M <- 3 ## No. of methylation states in a DNA sample type
N <- 4 ## No. of patients
R <- 2 ## No. of treatment conditions
my.seed <- 321 ## set seed for reproducibility

betaHMM_out <- betaHMM::betaHMM(pca_methylation_data,
                                annotation_data,
                                M = 3,
                                N = 4,
                                R = 2,
                                parallel_process = FALSE,
                                seed = my.seed,
                                treatment_group = c("Benign","Tumour"))

## ----betaHMMclass,include=TRUE, echo=TRUE-------------------------------------
class(betaHMM_out)

## ----betaHMMaccessor,include=TRUE, echo=TRUE----------------------------------
## transition matrix estimated for all chromosomes
A(betaHMM_out)

## Shape parameters estimated for a certain chromosome
phi(betaHMM_out)

## Hidden states assigned to all CpG sites for a certain chromosome
head(hidden_states(betaHMM_out)[["chr 7"]])

## ----betaHMMsummary,include=TRUE, echo=TRUE-----------------------------------
summary(betaHMM_out)

## ----dmc,include=TRUE, echo=TRUE----------------------------------------------
dmc_out <- dmc_identification(betaHMM_out)
dmc_df <- assay(dmc_out)
head(dmc_df)

## ----dmcsummary,include=TRUE, echo=TRUE---------------------------------------
summary(dmc_out)

## ----betaHMMplot,include=TRUE,echo=TRUE,fig.width=6.5,fig.height=5,dev='png'----
AUC_chr <- AUC(dmc_out)
plot(betaHMM_out, chromosome = "7", what = "fitted density", AUC = AUC_chr)

## ----betaHMMplot2,include=TRUE,echo=TRUE,fig.width=6,fig.height=5,dev='png'----
plot(betaHMM_out, chromosome = "7", what = "uncertainty",
        uncertainty_threshold = 0.2)

## ----dmr,include=TRUE, echo=TRUE----------------------------------------------
dmr_out <- dmr_identification(dmc_out, parallel_process = FALSE)
dmr_df <- assay(dmr_out)
head(dmr_df)

## ----dmrsummary,include=TRUE, echo=TRUE---------------------------------------
summary(dmr_out)

## ----dmrplot,include=TRUE,echo=TRUE,fig.width=7,fig.height = 5, dev = 'png'----
p<-plot(dmc_out, start_CpG = "cg17750844", end_CpG = 15)

## ----threshold,include=TRUE, echo=TRUE----------------------------------------
threshold_out <- threshold_identification(pca_methylation_data[,1:5],
                                            package_workflow = FALSE,
                                            annotation_file = annotation_data,
                                            M = 3,
                                            N = 4,
                                            parameter_estimation_only = FALSE,
                                            seed = my.seed)
threshold(threshold_out)

## ----thresholdplot,include=TRUE, echo=TRUE,fig.width = 5, fig.height = 4,dev = 'png'----
plot(threshold_out, what = "fitted density")

## -----------------------------------------------------------------------------
sessionInfo()

