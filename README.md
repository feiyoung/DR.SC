# DR.SC
DR-SC: Joint dimension reduction and spatial clustering for single-cell/spatial transcriptomics data 

DR.SC (Method name is DR-SC) is a package for analyzing  spatially resolved transcriptomics (SRT) datasets, developed by the Jin Liu's lab. It is not only computationally efficient and scalable to the sample size increment, but also is capable of choosing the smoothness parameter and the number of clusters as well.

Check out our [NAR paper](https://doi.org/10.1093/nar/gkac219) for a more complete description of the methods and analyses. 

DR.SC can be used to analyze experimental dataset from different technologies with different resolutions, for instance:

* ST plaform
* 10X Visium platform
* SeqFISH, MerFISH,  etc
* Slide-seq, Slide-seqV2, etc.
* Other platforms...

Once DR-SC model is fitted, the package provides functionality for further data exploration, 
analysis, and visualization. Users can:

* Identify clusters
* Extract low-dimensional embeddings
* Find significant  gene markers
* Visualize clusters and gene expression 2-dim using spatial coordinates or tSNE and UMAP

To further investigate transcriptomic properties, combining the results from DR.SC and other packages, users can:

* Infer the cell/domain lineages
* Infer RNA velocity if splicing and unsplicing matrix (can obtained from raw fastq data) are available
* Detect conditional spatially variational genes
* Conduct cell-deconvolution


## Installation

To install the the packages "DR.SC", firstly, install the 'remotes' package. Besides, "DR.SC" depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.
```{Rmd}
# Method 1: Install it from CRAN
install.packages("DR.SC")

# Method 2: Install it from github
install.packages("remotes")
remotes::install_github("feiyoung/DR.SC")
```


## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Toy example](https://feiyoung.github.io/DR.SC/articles/DR.SC.Simu.html)
* [DLPFC data analysis](https://feiyoung.github.io/DR.SC/articles/DR.SC.DLPFC.html)

## Setup on Linux or MacOS system
For running big data, users can use the following system command to set the C_stack unlimited in case of `R Error: C stack usage is too close to the limit`.
```{Linux}
ulimit -s unlimited
```

# Demonstration

For an example of typical DR.SC usage, please see our [Package vignette](https://feiyoung.github.io/DR.SC/index.html) for a demonstration and overview of the functions included in DR.SC.

