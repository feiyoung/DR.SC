\name{FindSVGs}
\alias{FindSVGs}
\title{Find spatially variable genes}
\description{
  Identifies features that have spatially variation along spots using SPARK-X.
}
\usage{
  FindSVGs(seu, nfeatures=2000, covariates=NULL, num_core=1, verbose=TRUE)
}
\arguments{
  \item{seu}{an object of class "Seurat".}
  \item{nfeatures}{a positive integer, means how many spatially variable genes to be chosen. If there are less than 2000 features in seu, then all features are identified.}
  \item{covariates}{a covariate matrix named control variable matrix whose number of rows is equal to the number of columns of seu.}
   \item{num_core}{an optional positive integer, specify the cores used for identifying the SVGs in parallel.}
   \item{verbose}{an optional logical value, whether output the related information.}
}
\details{
  Nothing
}
\value{
  return a revised Seurat object by adding three columns named "is.SVGs", "order.SVGs" and "adjusted.pval.SVGs" in the meta.features  of default Assay.
}
\references{
   Zhu, J., Sun, S., Zhou, X.: Spark-x: non-parametric modeling enables scalable and robust detection of spatialexpression patterns for large spatial transcriptomic studies. Genome Biology 22(1), 1-25 (2021)
}
\note{
  nothing
}
\seealso{
  \link{topSVGs}
}
\examples{
  seu<-gendata_RNAExp(height=20, width=20,p=200, K=4)
  seu<-FindSVGs(seu, nfeatures=100)
  topSVGs(seu)

}