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
  seu <- gendata_RNAExp(height=20, width=20,p=100, K=4)
  library(Seurat)
  seu <- NormalizeData(seu)
  # choose spatially variable features
  seu <- FindSVGs(seu)
  ## Just for illustrating the usage of mbicPlot
  seu@tools$icMat <- data.frame(K=2:5, mbic=c(105, 101, 99, 108))
  mbicPlot(seu)

}
\concept{Modified BIC}
