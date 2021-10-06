\name{FindSVGs}
\alias{FindSVGs}
\title{Find spatially variable genes}
\description{
  Identifies features that have spatially variation along spots.
}
\usage{
  FindSVGs(seu, nfeatures=2000,num_core=1, verbose=T)
}
\arguments{
  \item{seu}{an object of class "Seurat".}
  \item{nfeatures}{a positive integer, means how many highly variable or spatially variable genes to be chosen. If there are less than 2000 features in seu, then all features are identified.}
   \item{num_core}{an optional positive integer, specify the cores used for identifying the SVGs in parallel.}
   \item{verbose}{an optional logical value, whether output the mediate information.}
}
\details{
  Nothing
}
\value{
  return a revised Seurat object by adding a column named "is.SVGs" in the meta.features  of default Assay.
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
  topSVGs(seu)
}
\keyword{Input of spatial data}
