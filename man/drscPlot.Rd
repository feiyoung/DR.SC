\name{drscPlot}
\alias{drscPlot}
\title{tNSE or UMAP plot visualization}
\description{
  Intuitive way of visualizing how cell types changes across the spatial locations.
}
\usage{
  drscPlot(seu, dims=1:5, visu.method='tSNE',...)
}
\arguments{
 \item{seu}{an object of class "Seurat" obtained by \link{DR.SC}.}
 \item{dims}{a positive integer to specify the number of latent features for visualiztion.}
 \item{visu.method}{a string including 'tSNE' or "UMAP".}
 \item{...}{Other arguments passing to \link{DimPlot} function.}
}
\details{
  Nothing
}
\value{
 return a ggplot2 object.
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
seu <- gendata_RNAExp(height=10, width=10,p=50, K=4)
library(Seurat)
seu <- NormalizeData(seu)
# choose spatially variable features
seu <- FindSVGs(seu)

# use SVGs to fit DR.SC model
# maxIter = 2 is only used for illustration, and user can use default.
seu1 <- DR.SC(seu, K=4,platform = 'ST', variable.type = 'SVGs',maxIter = 2,verbose=FALSE)
drscPlot(seu1)
}
\keyword{tSNE}
\keyword{UMAP}
\concept{DR-SC}
