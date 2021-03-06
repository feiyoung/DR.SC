# DR.SC
DR.SC: Joint dimension reduction and spatial clustering for single-cell/spatial transcriptomics data 

# Installation

To install the the packages "DR.SC", firstly, install the 'remotes' package. Besides, "DR.SC" depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.
```{Rmd}
# Method 1: Install it from CRAN
install.packages("DR.SC")

# Method 2: Install it from github
install.packages("remotes")
remotes::install_github("feiyoung/DR.SC")
```
## Setup on Linux or MacOS system
For running big data, users can use the following system command to set the C_stack unlimited in case of `R Error: C stack usage is too close to the limit`.
```{Linux}
ulimit -s unlimited
```

# Demonstration

For an example of typical DR.SC usage, please see our [Package vignette](https://feiyoung.github.io/DR.SC/index.html) for a demonstration and overview of the functions included in DR.SC.

