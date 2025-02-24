# README

Author: Martin Papenberg
Year: 2025

This repository contains additional documentation, code and data for the manuscript

"Anticlustering for Sample Allocation To Minimize Batch Effects"

It also has the code for the [anticlust Shiny App](https://comphealth.ucsf.edu/app/anticlustering/) that was first described in that paper. 

- Directory `Supplementary_Materials` has a document that has supplementary information and analysis on the project. It is a reproducible (via R markdown) document that generates all analyses, Figures and Tables reported in the paper, as well as some additional analyses not reported in the paper. 
- Directory `Application` contains code and data to reproduce the example application.
- Directory `Simulation` implements a simulation that compares the quality of the anticlustering method and alternative approaches. It contains code to do the simulation, data from the simulation, and code to analyze the simulation data .
- Directory `Running_Time_Optimal_Algorithm` has code that investigates the feasibility of an optimal algorithm for must-link anticlustering
- Directory `2PML` contains an R script that solves the (NP-hard) subset sum problem for up to 27 elements via complete enumeration, to generate combinations of numbers that are used in Phase 2 of the 2PML algorithm
- Directory `Web app` has the code for the Shiny App.
- Directory `Shiny_Interface` includes the function that is used for anticlustering in the Shiny App. It "smartly" selects an anticlustering method depending on the input (primarily, depending on the size of the data set). 

## Dependencies

The analyses in the paper were generated using R version 4.4.2 on an Intel i7-10700 personal computer (4.800GHz x 8) with 16 GB RAM, running Ubuntu 20.04.6 LTS.

Additionally, the following R packages / scripts were used in the analyses:

**R packages**

- `anticlust` (version 0.8.9-1)
  * `install.packages("anticlust")`
- `OSAT` (version 1.54.0; for Simulation Study)
  * [Available from Bioconductor](https://bioconductor.org/packages/release/bioc/html/OSAT.html)
  * `if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")`
  * `BiocManager::install("OSAT")`
- Propensity Score Batch Assignment
  * Functions available from https://github.com/carryp/PS-Batch-Effect
  * Last commit from when I accessed the files: d9ba5cd0d9a240e761d045a99557b715355ddfb8
- `here` (version 1.0.1; for handling relative paths in this project; the working directory is always assumed to be the top-level directory)
  * `install.packages("here")`
- `dplyr` (version 1.1.4)
  * `install.packages("dplyr")` 
- `rmarkdown` (version 2.28)
  * `install.packages("rmarkdown")` 
- `knitr` (version 1.48)
  * `install.packages("knitr")` 
- `tableone` (version 0.13.2)
  * `install.packages("tableone")` 
- `kableExtra` (version 1.4.0)
  * `install.packages("kableExtra")` 
- `santoku` (version 1.0.0)
  * `install.packages("santoku")` 
- `tidyr` (version 1.3.1)
  * `install.packages("tidyr")` 
- `ggplot2` (version 3.5.1)
  * `install.packages("ggplot2")` 

**Additional dependencies for `Running_Time_Optimal_Algorithm`**

In the directory `Running_Time_Optimal_Algorithm`, we evaluate an optimal algorithm for anticlustering using must-link constraints. To use the optimal method, we need a "solver" for integer linear programming. We use [gurobi](https://www.gurobi.com/), which generally outperforms open source solvers (i.e., it can be used to process larger data sets). However, gurobi is generally not free to use and requires a license for usage. For academics, [free to use licenses are available](https://www.gurobi.com/academia/academic-program-and-licenses/). The gurobi software also ships an R package gurobi, which has to be installed to be used as solver for optimal anticlustering algorithms in `anticlust`. To use gurobi for anticlustering, however, we need to install `anticlust` from a separate branch on Github, and I do not expect that the anticlust that supports gurobi will be distributed via CRAN:

```
if (!require("remotes", quietly = TRUE)) 
  install.packages("remotes")
remotes::install_github("m-Py/anticlust", ref = "gurobi")
```

I made the analysis using version 0.8.9-1.604 (latest commit 9ccde8bc1ec2c925cb9340ab2710b14c382b96ca), which is version 0.8.9-1 with the only difference that it can use the gurobi solver as backend for ILP algorithms.

