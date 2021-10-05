\name{DR.SC}
\alias{DR.SC}
\title{Joint dimension reduction and spatial clustering}
\description{
  Joint dimension reduction and spatial clustering for scRNA-seq and spatial transcriptomics data
}
\usage{
  DR.SC(object, ...)

### S3 method for class "Seurat"
  DR.SC(seu, q=15, K=NULL, platform= c("Visium", "ST", 'scRNAseq'),
                         nfeatures=2000,K_set = seq(2, 10), ...)
}
\arguments{
  \item{seu}{an object of class "Seurat". The details of this object are given under 'Details'.}
  \item{q}{a positive integer, specify the number of latent features to be extracted, default as 15.}
  \item{K}{a positive  integer, specify the number of clusters, default as \code{NULL}. When \code{K=NULL}, it is automatically selected by MBIC criteria.}
  \item{platform}{a string, specify the platform of the provided data, default as "Visium". There are three platforms to be chosen, including "Visuim", "ST" and "scRNAseq", where the first two means there are spatial coordinates information in the metadata of seu, named "row" and "col" and a Hidden Markov random field is used to model the spatial coordinates, and the last one "scRNAseq" means there is no spatial information in object seu and a multinomial model is used to model the unobserved class labels. The platform helps to calculate the adjoint matrix.}
  \item{nfeatures}{a positive integer, means how many highly variable or spatially variable genes used for DR-SC model, default as 2000. If there are less than 2000 features in seu, then all features are used for DR-SC model.}
  \item{K_set}{a vector of positive integer, means the candidates of number of clusters used for MBIC.}
  \item{...}{Other arguments to pass into \link{DR.SC.fit} function.}
}
\details{
  seu is an object named \code{Seurat}, thich can easily created by R package \link{Seurat}.
  If the data is collected by the spatial transcriptomics technologies such as 10X Visium and ST, then there are spatial coordinates information in the metadata of seu, named "row" and "col". DR-SC model uses a Hidden Markov random field to model the spatial coordinates.  If the data is collected by the single cell RNA sequencing techonologies which means there is no spatial information in object seu then a multinomial model is used to model the unobserved class labels.

}
\value{
  DR.SC returns a revised  \code{Seurat} object. There are two revisions in the seu. 1. the metadata is added a new column named \code{spatial.drsc.cluster} that represents the clustering results from DR-SC model. 2. a DimReduc object named \code{dr-sc} is added in the slot \code{reductions}, which represents the features extracted by DR-SC model.

}
\references{
  
}
\author{
Wei Liu
}
\note{
nothing
}
\seealso{

}
\examples{
## we generate the spatial transcriptomics data with lattice neighborhood, i.e. ST platform.
seu <- gendata_spatial(height=30, width=30,p=500, K=4)
library(Seurat)
seu <- NormalizeData(seu)
# choose 2000 variable features using Seurat
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu1 <- DR.SC(seu, K=4, platform = 'ST', verbose=F)
}
\keyword{Dimension reduction, spatial clustering and cell type identification.}
