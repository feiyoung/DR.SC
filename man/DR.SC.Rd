\name{DR.SC}
\alias{DR.SC}
\alias{DR.SC.Seurat}
\title{Joint dimension reduction and spatial clustering}
\description{
  Joint dimension reduction and spatial clustering for scRNA-seq and spatial transcriptomics data
}
\usage{

\method{DR.SC}{Seurat}(seu, K, q=15,  platform= c('Visium', "ST", "Other_SRT", "scRNAseq"),...)
}
\arguments{
  \item{seu}{an object of class "Seurat". The details of this object are given under 'Details'.}
  \item{q}{a positive integer, specify the number of latent features to be extracted, default as 15.}
  \item{K}{a positive  integer or integer vector, specify the number of clusters. When \code{K} is a vector, it is automatically selected by MBIC criteria. Users can also use BIC and AIC  to select K or refine model with MBIC by argument \code{pen.const} in the function \link{selectModel}.}
  \item{platform}{a string, specify the platform of the provided data, default as "Visium". There are more platforms to be chosen, including ("Visuim", "ST", "Other_SRT") and "scRNAseq", where the first group means there are spatial coordinates information in the metadata of seu, named "row" and "col" and a Hidden Markov random field is used to model the unobserved class label using spatial coordinates ("Other_SRT" represents the other SRT platforms except for 'Visium' and 'ST'), and the other group "scRNAseq" means there is no spatial information in object seu and a multinomial model is used to model the unobserved class labels. The platform helps to calculate the adjacency matrix.}
  \item{...}{Other arguments to pass into \link{DR.SC_fit} function.}
}
\details{
seu is an object named \code{Seurat}, thich can easily created by R package \link{Seurat}.
DR-SC  model can be applied to analyze both single cell RNA sequencing data and spatially resoved transcriptomics (SRT) data. 

If the data are collected by the single cell RNA sequencing techonologies which means there is no spatial information in object seu then a multinomial model is used to model the unobserved class labels. 


If the data is collected by the spatial transcriptomics technologies, then there are spatial coordinates information in the metadata of seu, named "row" and "col". DR-SC model uses a Hidden Markov random field to model the spatial coordinates.  \code{DR.SC} supports different platforms of SRT data, such as `Visium`, `ST` and any other platforms `Other_SRT`.

For lattice grids in ST platform (ST), the interior spot has four neighbors (left, right, up and down),the boundary spot has three neighbors, and the spot in the corner has two neighbors. For hexagon grids, such as spatial coordinate in 10X Visium platform (Visium), the interior spot has six neighbors. For the irregular coordinates in other platforms (Other_SRT), Euclidean distance is adopted to decide whether a spot is an neighbor of another spot. For example, if the Euclidean distance between spot A and B is less than a radius, then A is taken as the neighbourhood of B. See function \link{getAdj} for more details.

}
\value{
  DR.SC returns a revised  \code{Seurat} object. There are two revisions in the seu: 1. the metadata is added a new column named \code{spatial.drsc.cluster} that represents the clustering results from DR-SC model, and the \code{Idents(seu)} is assigned with \code{spatial.drsc.cluster}. 2. a DimReduc object named \code{dr-sc} is added in the slot \code{reductions}, which represents the features extracted by DR-SC model.

}
\references{
 \href{https://academic.oup.com/nar/article/50/12/e72/6555431}{Wei Liu, Xu Liao, Yi Yang, Huazhen Lin, Joe Yeong, Xiang Zhou, Xingjie Shi & Jin Liu (2022). Joint dimension reduction and clustering analysis of single-cell RNA-seq and spatial transcriptomics data, Nucleic Acids Research, gkac219.}
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
