\name{DR.SC}
\alias{DR.SC}
\alias{DR.SC.Seurat}
\title{Joint dimension reduction and spatial clustering}
\description{
  Joint dimension reduction and spatial clustering for scRNA-seq and spatial transcriptomics data
}
\usage{

\method{DR.SC}{Seurat}(seu, K, q=15,  platform= "Visium",...)
}
\arguments{
  \item{seu}{an object of class "Seurat". The details of this object are given under 'Details'.}
  \item{q}{a positive integer, specify the number of latent features to be extracted, default as 15.}
  \item{K}{a positive  integer or integer vector, specify the number of clusters. When \code{K} is a vector, it is automatically selected by MBIC criteria. Users can also use BIC and AIC  to select K or refine model with MBIC by argument \code{pen.const} in the function \link{selectModel}.}
  \item{platform}{a string, specify the platform of the provided data, default as "Visium". There are more platforms to be chosen, including ("Visuim", "ST", "seqfish", 'merfish', 'slide-seqv2', 'seqscope', 'hdst') and "scRNAseq", where the first group means there are spatial coordinates information in the metadata of seu, named "row" and "col" and a Hidden Markov random field is used to model the spatial coordinates, and the last one "scRNAseq" means there is no spatial information in object seu and a multinomial model is used to model the unobserved class labels. The platform helps to calculate the adjacency matrix.}
  \item{...}{Other arguments to pass into \link{DR.SC_fit} function.}
}
\details{
  seu is an object named \code{Seurat}, thich can easily created by R package \link{Seurat}.
If the data is collected by the spatial transcriptomics technologies such as 10X Visium, ST, seqFISH, MERFISH and  Slide-seq, then there are spatial coordinates information in the metadata of seu, named "row" and "col". DR-SC model uses a Hidden Markov random field to model the spatial coordinates.  If the data is collected by the single cell RNA sequencing techonologies which means there is no spatial information in object seu then a multinomial model is used to model the unobserved class labels. 

}
\value{
  DR.SC returns a revised  \code{Seurat} object. There are two revisions in the seu. 1. the metadata is added a new column named \code{spatial.drsc.cluster} that represents the clustering results from DR-SC model, and the \code{Idents(seu)} is assigned with \code{spatial.drsc.cluster}. 2. a DimReduc object named \code{dr-sc} is added in the slot \code{reductions}, which represents the features extracted by DR-SC model.

}
\references{
 \href{https://www.biorxiv.org/content/10.1101/2021.12.25.474153v4}{Wei Liu, Xu Liao, Yi Yang, Huazhen Lin, Joe Yeong, Xiang Zhou, Xingjie Shi and Jin Liu. (2022) Joint dimension reduction and clustering analysis for single-cell RNA-seq and spatial transcriptomics data}
}
\author{
Wei Liu
}
\note{
nothing
}
\seealso{
 None
}
\examples{

## we generate the spatial transcriptomics data with lattice neighborhood, i.e. ST platform.
seu <- gendata_RNAExp(height=10, width=10,p=50, K=4,platform="ST")
library(Seurat)
seu <- NormalizeData(seu, verbose=FALSE)
# choose 100 highly variable features
# seu <- FindVariableFeatures(seu, nfeatures = 100)
# maxIter = 2 is only used for illustration, and user can use default.
# seu1 <- DR.SC(seu, K=4, platform = 'ST', maxIter=2,verbose=FALSE)

# choose spatially variable features (SVGs)
seu <- FindSVGs(seu, nfeatures = 40, verbose=FALSE)
# use SVGs to fit DR.SC model
# maxIter = 2 is only used for illustration, and user can use default.
seu1 <- DR.SC(seu, K=4,platform = 'ST', maxIter=2, verbose=TRUE)

}
\keyword{Dimension reduction}
\keyword{Spatial clustering}
