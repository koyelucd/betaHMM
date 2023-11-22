## Use simulated data for the betaHMM workflow example
set.seed(12345)

## read files
data(sample_methylation_file)
data(sample_annotation_file)
# Run betaHMM function
beta_out <- betaHMM(sample_methylation_file[1:50,],
                    sample_annotation_file[1:50,],
                    M = 3, N = 4, R = 2,iterations=2,
                    parallel_process = FALSE, seed = 12345,
                    treatment_group = c("Benign","Tumour"))

## Run dmc_identification function
dmc_out <- dmc_identification(beta_out)

# Run dmr_identification function
dmr_out <- dmr_identification(dmc_out, parallel_process = FALSE)

# Plot functions
# Get the AUC values calculated for each hidden state
AUC_chr <- AUC(dmc_out)

## plot the uncertainty for each hidden state
plot(beta_out, chromosome = "1", what = "uncertainty")

