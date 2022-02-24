\name{spatialPlotClusters}
\alias{spatialPlotClusters}
\title{Spatial coordinates plot visualization}
\description{
  Intuitive way of visualizing how cell types changes across the spatial locations.
}
\usage{
  spatialPlotClusters(seu)
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
    # choose spatially variable features using Seurat
    seu <- FindSVGs(seu)
    # use SVGs to fit DR.SC model
    # maxIter = 2 is only used for illustration, and user can use default.
    seu1 <- DR.SC(seu, K=4,platform = 'ST', maxIter=2,verbose=FALSE)
    spatialPlotClusters(seu1)
}
\keyword{Spatial plot}
