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
seu <- gendata_RNAExp(height=30, width=30,p=500, K=4)
library(Seurat)
seu <- NormalizeData(seu)
# choose spatially variable features
seu <- FindSVGs(seu)
# use SVGs to fit DR.SC model
seu1 <- DR.SC(seu, K_set=2:6,platform = 'ST',variable.type = 'SVGs', verbose=F,pen.const=0.5)
mbicPlot(seu1)
}
\keyword{Dimension reduction, spatial clustering and cell type identification.}
