# Simulation study

Author: Martin Papenberg
Year: 2025

This directory contains code and data to reproduce the simulation reported in the manuscript

"Anticlustering for Sample Allocation To Minimize Batch Effects"

The script `Simulation.R` implements the simulation, loading some functions defined in `functions.R`. It sets a seed via `set.seed()` for exact reproducibility. Note that the aggregated results reported in the paper do not depend on a seed, so feel free to turn it off when playing with the code. 

The script `analysis_simulation.R` implements the analysis reported in the paper (the Rmd file in subdirectory Supplementary_Materials also has this code, slightly modified so that the values computed in R can directly be inserted into the PDF file).

## Codebook

The file `results.csv` contain the results of the 10000 simulation runs. Each row corresponds to one data set that was generated during the simulation and which was processed via anticlust, OSAT, and propensity score batch assignment. It has the following columns:

- "N" (the number of samples/rows in the data set that was generated), varies between 50 and 500; valid values depend on "K" because only "N"s were generated that were divisible by "K".
- "M" (the number of variables/columns in the data set that was generated), varies between 2 and 5.
- "K" (the number of batches), is 2, 4, or 10
- "P" (the number of categories per categorical variable), varies between 2 and 5.
- "p_osat1" - "p_osat5" (*p* values for the OSAT method for each of the 2-5 variables; "p_osat3" - "p_osat5" can be NA)
- "p_anticlust1" - "p_anticlust5" (*p* values for the anticlustering method for each of the 2-5 variables; "p_anticlust3" - "p_anticlust5" can be NA)
- "p_anticlust_c1" - "p_anticlust_c5" (*p* values for the must-link anticlustering method for each of the 2-5 variables; "p_anticlust_c3" - "p_anticlust_c5" can be NA)
- "p_ps1" - "p_ps5" (*p* values for the propensity score method for each of the 2-5 variables; "p_ps3" - "p_ps5" can be NA)
- "OSAT_t" (time in seconds for OSAT)
- "ANTICLUST_t" (time in seconds for unconstrained anticlustering)
- "ANTICLUST_t_c" (time in seconds for must-link anticlustering) 
- "PS_t" (time in seconds for propensity score assignment)
- "diversity_unconstrained" (diversity of the unconstrained anticlustering assignment)
- "diversity_constrained"  (diversity of the must-link anticlustering assignment)
