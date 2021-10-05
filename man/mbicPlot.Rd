\name{mbicPlot}
\alias{mbicPlot}
\title{MBIC plot visualization}
\description{
  Intuitive way of visualizing how modified BIC values changes across different number of clusters
}
\usage{
  mbicPlot(seu)
}
\arguments{
 \item{seu}{an object of class "Seurat" revised by \link{DR.SC} with argument \code{K=NULL}.}
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
seu <- gendata_spatial(height=30, width=30,p=500, K=4)
library(Seurat)
seu <- NormalizeData(seu)
# choose 2000 variable features using Seurat
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu1 <- DR.SC(seu, K=NULL, K_set=2:6,platform = 'ST', verbose=F)
mbicPlot(seu1)
}
\keyword{Dimension reduction, spatial clustering and cell type identification.}
