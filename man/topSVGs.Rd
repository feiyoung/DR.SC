\name{topSVGs}
\alias{topSVGs}
\title{Return the top n SVGs}
\description{
  Return top n spatially variable genes given a Seurat object performed by \link{FindSVGs}.
}
\usage{
  topSVGs(seu, ntop=5)
}
\arguments{
  \item{seu}{an object of class "Seurat".}
  \item{ntop}{an optional positive integer, means how many spatially variable genes to access. }
}
\details{
  Nothing
}
\value{
  return a \link{character} vector including the names of SVGs.
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
  \link{topSVGs}
}
\examples{
  seu <- gendata_RNAExp(height=20, width=20,p=200, K=4)
  seu <- FindSVGs(seu, nfeatures=100)
  (genes <- topSVGs(seu, ntop=10))
  
}
\keyword{SVGs}
