\name{FindSVGs}
\alias{FindSVGs}
\title{Find spatially variable genes}
\description{
  Identifies features that have spatially variation along spots.
}
\usage{
  FindSVGs(seu, nfeatures=2000,num_core=1, verbose=TRUE)
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
  return a revised Seurat object by adding three columns named "is.SVGs", "order.SVGs" and "adjusted.pval.SVGs" in the meta.features  of default Assay.
}
\references{
   Zhu, J., Sun, S., Zhou, X.: Spark-x: non-parametric modeling enables scalable and robust detection of spatialexpression patterns for large spatial transcriptomic studies. Genome Biology 22(1), 1–25 (2021)
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
  topSVGs(seu)

}
\keyword{Input of spatial data}
