\name{drscPlot}
\alias{drscPlot}
\title{tNSE or UMAP plot visualization}
\description{
  Intuitive way of visualizing how cell types changes across the spatial locations.
}
\usage{
  drscPlot(seu)
}
\arguments{
 \item{seu}{an object of class "Seurat" obtained by \link{DR.SC}.}
}
\details{
  Nothing
}
\value{
 return a ggplot2 object.
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
seu <- gendata_RNAExp(height=30, width=30,p=500, K=4)
library(Seurat)
seu <- NormalizeData(seu)
# choose spatially variable features
seu <- FindSVGs(seu)
# use SVGs to fit DR.SC model
seu1 <- DR.SC(seu, K=4,platform = 'ST', variable.type = 'SVGs',verbose=F)
drscPlot(seu1)
}
\keyword{DR-SC, tSNE, UMAP.}
