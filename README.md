# must-link-anticlustering

- Directory `Application` contains code for the example application (synthetic data provided by Tomiko Oskotsky)
- Directory `Simulation_vs_OSAT` implements a simulation that compares the quality of the OSAT method, the propensity score method, and anticlustering
- Directory `Supplementary_Materials` has a document that has supplementary information and analysis on the project. 
- Directory `Running_Time_Optimal_Algorithm` has code that investigates the feasibility of an optimal algorithm for must-link anticlustering

## Dependencies

The following R packages are needed to reproduce the analyses / documents:

**R packages**

- `anticlust`
  * I used version 0.8.7, which is available from CRAN.
- `OSAT` (For Simulation Study)
  * [Available from Bioconductor](https://bioconductor.org/packages/release/bioc/html/OSAT.html)
  * `if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")`
  * `BiocManager::install("OSAT")`
- `here` (for handling relative paths in this project; the working directory is always assumed to be the top-level directory)
  * `install.packages("here")`
- tidyverse packages (for data analysis)
  * `install.packages("tidyverse")` 
- `rmarkdown` (for rendering the supplementary materials as PDF)
  * `install.packages("rmarkdown")` 

**Additional dependencies for `Running_Time_Optimal_Algorithm`**

In the directory `Running_Time_Optimal_Algorithm`, we evaluate an optimal algorithm for anticlustering using must-link constraints. To use the optimal method, we need a "solver" for integer linear programming. We use [gurobi](https://www.gurobi.com/), which generally outperforms open source solvers (i.e., it can be used to process larger data sets). However, gurobi is generally not free to use and requires a license for usage. For academics, [free to use licenses are available](https://www.gurobi.com/academia/academic-program-and-licenses/). The gurobi software also ships an R package gurobi, which has to be installed to be used as solver for optimal anticlustering algorithms in `anticlust`. To use gurobi for anticlustering, however, we need to install `anticlust` from a separate branch on Github, and I do not expect that the anticlust that supports gurobi will be distributed via CRAN:

```
if (!require("remotes", quietly = TRUE)) 
  install.packages("remotes")
remotes::install_github("m-Py/anticlust", ref = "gurobi")
```

