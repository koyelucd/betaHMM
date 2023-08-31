## Use simulated data for the betaHMM workflow example
set.seed(12345)
library(betaHMM)

## read files
data(sample_methylation_file)
head(sample_methylation_file)
data(sample_annotation_file)
head(sample_annotation_file)
## Run betaHMM function
beta_out <- betaHMM(sample_methylation_file, sample_annotation_file,
                    M = 3, N = 4, R = 2,
                    parallel_process = FALSE, seed = 12345,
                    treatment_group = c("Benign","Tumour"))

## Run dmc_identification function
dmc_out <- dmc_identification(beta_out)
dmc_df <- assay(dmc_out)

# Run dmr_identification function
dmr_out <- dmr_identification(dmc_out, parallel_process = FALSE)
#
# Plot functions
# Get the AUC values calculated for each hidden state
AUC_chr <- AUC(dmc_out)

## plot the fitted density estimates
plot(beta_out, chromosome = "1", what = "fitted density", AUC = AUC_chr)

## Plot the DMR graph
plot(dmc_out, start_CpG = "cg14817997", end_CpG = 8)
