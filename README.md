# betaHMM: a hidden Markov model to identify differential methylation sites and
 regions from beta-valued methylation data 

Authors: KOYEL MAJUMDAR, ISOBEL CLAIRE GORMLEY, THOMAS BRENDAN MURPHY, ROMINA
SILVA, ANTOINETTE SABRINA PERRY, RONALD WILLIAM WATSON,
FLORENCE JAFFR ́EZIC, ANDREA RAU

A hidden Markov model approach for identifying differentially methylated 
sites and regions for beta-valued DNA methylation data. The workflow consists
of 3 main functions:
\itemize{
\item betaHMM: to estimate the hidden Markov model parameters and the 
hidden methylation states of each CpG sites.
\item dmc_identification: uses the output from the above function to identify
CpG sites that are mostly differentially methylated between DNA samples 
collected from \eqn{R} conditions.
\item dmr_identification: to identify adjacent DMCs forming DMRs.
}

A typical call to betaHMM to apply the Baum-Welch algorithm and Viterbi 
algorithm takes the following form:
```
library(betaHMM)
### read the methylation file and annotation file
data(sample_methylation_file)
head(sample_methylation_file)
data(sample_annotation_file)
head(sample_annotation_file)
betaHMM_out <- betaHMM(sample_methylation_file,
                                sample_annotation_file,
                                M = 3,
                                N = 4,
                                R = 2,
                                parallel_process = FALSE,
                                seed = 321,
                                treatment_group = c("Benign","Tumour"))
```
where `M` represents the number of methylation state in a DNA sample, `N`
represents the number of subjects or DNA replicates, `R` represents the number
of conditions from where DNA is extracted.

A typical call to dmc_identification takes the following form:
```
dmc_out <- dmc_identification(betaHMM_out)
dmc_df <- assay(dmc_out)
```

A typical call to dmr_identification takes the following form:
```
dmr_out <- dmr_identification(dmc_out)
dmr_df <- assay(dmr_out)
```


The output of the `betaHMM` function is an
S4 object of class `betaHMMResults`, the output of the `dmc_identification`
function is an S4 object of class `dmcResults` and the output of the 
`dmr_identification` function is an S4 object of class `dmrResults` 
on which standard `plot` and `summary` functions can be directly applied; 
the former uses functionalities from the 
[ggplot2]( https://cran.r-project.org/package=ggplot2) package.

### Reference


### License

The betaHMM package is free software; you can redistribute it and/or modify 
it under the terms of the GNU General Public License, version 3, 
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without 
any warranty; without even the implied warranty of merchantability or fitness 
for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, 
is available at http://www.r-project.org/Licenses/GPL-3.
