\name{mbicPlot}
\alias{mbicPlot}
\title{MBIC plot visualization}
\description{
  Intuitive way of visualizing how modified BIC values changes across different number of clusters
}
\usage{
  mbicPlot(seu, criteria="MBIC")
}
\arguments{
 \item{seu}{an object of class \code{Seurat} revised by \link{DR.SC} with argument \code{K=NULL}.}
 \item{criteria}{a string specifying the information criteria such as AIC, BIC and MBIC to be plotted, default as MBIC.}
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
  seu[["RNA"]]@misc[['icMat']] <- data.frame(K=2:5, MBIC=c(105, 101, 99, 108))
  mbicPlot(seu)

}
\concept{Modified BIC}
