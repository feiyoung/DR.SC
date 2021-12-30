\name{DR.SC}
\alias{DR.SC}
\alias{DR.SC.Seurat}
\title{Joint dimension reduction and spatial clustering}
\description{
  Joint dimension reduction and spatial clustering for scRNA-seq and spatial transcriptomics data
}
\usage{
### S3 method for class "Seurat"
DR.SC(seu, q=15, K=NULL, platform= "Visium",
       nfeatures=2000,K_set = seq(2, 10), variable.type="HVGs", ...)
}
\arguments{
  \item{seu}{an object of class "Seurat". The details of this object are given under 'Details'.}
  \item{q}{a positive integer, specify the number of latent features to be extracted, default as 15.}
  \item{K}{a positive  integer, specify the number of clusters, default as \code{NULL}. When \code{K=NULL}, it is automatically selected by MBIC criteria.}
  \item{platform}{a string, specify the platform of the provided data, default as "Visium". There are three platforms to be chosen, including ("Visuim", "ST", "seqfish", 'merfish', 'slide-seqv2', 'seqscope') and "scRNAseq", where the first five means there are spatial coordinates information in the metadata of seu, named "row" and "col" and a Hidden Markov random field is used to model the spatial coordinates, and the last one "scRNAseq" means there is no spatial information in object seu and a multinomial model is used to model the unobserved class labels. The platform helps to calculate the adjoint matrix.}
  \item{nfeatures}{a positive integer, means how many highly variable or spatially variable genes used for DR-SC model, default as 2000. If there are less than 2000 features in seu, then all features are used for DR-SC model.}
  \item{K_set}{an optional vector of positive integer, means the candidates of number of clusters used for MBIC.}
  \item{variable.type}{an optional string, specify whether use highly variable genes (HVGs) or spatially variable genes (SVGs) to fit DR-SC model, default as "HVGs".}
  \item{...}{Other arguments to pass into \link{DR.SC_fit} function.}
}
\details{
  seu is an object named \code{Seurat}, thich can easily created by R package \link{Seurat}.
  If the data is collected by the spatial transcriptomics technologies such as 10X Visium, ST, seqFISH, MERFISH and  Slide-seq, then there are spatial coordinates information in the metadata of seu, named "row" and "col". DR-SC model uses a Hidden Markov random field to model the spatial coordinates.  If the data is collected by the single cell RNA sequencing techonologies which means there is no spatial information in object seu then a multinomial model is used to model the unobserved class labels.

}
\value{
  DR.SC returns a revised  \code{Seurat} object. There are two revisions in the seu. 1. the metadata is added a new column named \code{spatial.drsc.cluster} that represents the clustering results from DR-SC model. 2. a DimReduc object named \code{dr-sc} is added in the slot \code{reductions}, which represents the features extracted by DR-SC model.

}
\references{
  None
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
seu <- NormalizeData(seu)
# choose 100 highly variable features
# seu <- FindVariableFeatures(seu, nfeatures = 100)
# maxIter = 2 is only used for illustration, and user can use default.
# seu1 <- DR.SC(seu, K=4, platform = 'ST', maxIter=2,verbose=FALSE)

# choose spatially variable features
seu <- FindSVGs(seu, nfeatures = 40)
# use SVGs to fit DR.SC model
# maxIter = 2 is only used for illustration, and user can use default.
seu1 <- DR.SC(seu, K=4,platform = 'ST', maxIter=2,variable.type = 'SVGs',verbose=FALSE)

}
\keyword{Dimension reduction}
\keyword{Spatial clustering}
