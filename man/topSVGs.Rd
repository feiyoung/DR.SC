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
  
}
\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  \link{topSVGs}, \link{SPARK}, \link{sparkx}
}
\examples{
  seu <- gendata_RNAExp(height=40, width=40,p=1000, K=4)
  seu <- FindSVGs(seu, nfeatures=800)
  (genes <- topSVGs(seu, ntop=10))
}
\keyword{SVGs}
