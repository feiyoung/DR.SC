\name{selectModel}
\alias{selectModel}
\alias{selectModel.drscObject}
\alias{selectModel.Seurat}
\title{Select the number of clusters}
\description{
  Select the number of clusters by specified criteria.
}
\usage{
  
  selectModel(obj, criteria = 'MBIC', pen.const=1)
  \method{selectModel}{drscObject}(obj, criteria = 'MBIC', pen.const=1)
  \method{selectModel}{Seurat}(obj, criteria = 'MBIC', pen.const=1)
  
  
}
\arguments{S
  \item{obj}{an object with class \code{Seurat} by \link{DR.SC} or class \code{drscObject} by \link{DR.SC_fit}.}
  \item{criteria}{a string, specify the criteria used for selecting the number of clusters, supporting "MBIC", "BIC" and "AIC".}
  \item{pen.const}{an optional positive value, the adjusted constant used in the MBIC criteria. It usually takes value between 0.1 to 1.}
}

\value{
   For S3 method of \code{Seurat}, it return a  revised "Seurat" object with updated \code{Idents(seu)}, \code{spatial.drsc.cluster} in the metadata and DimReduc object named \code{dr-sc} in the slot \code{reductions}. For S3 method of \code{drscObject}, it returns a list with the following components:
   \item{bestK}{the selected number of clusters.}
   \item{cluster}{inferred class labels}
   \item{hZ}{extracted latent features.}
   \item{icMat}{a numeric matrix including the criteria value for each number of clusters K.}
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
  \link{DR.SC}, \link{DR.SC_fit}.
}
\examples{
seu <- gendata_RNAExp(height=10, width=10,p=50, K=4)
library(Seurat)
seu <- NormalizeData(seu, verbose=FALSE)
# or choose 40 spatailly variable features using FindSVGs in DR.SC
seu <- FindSVGs(seu, nfeatures = 40, verbose=FALSE)
# users define the adjacency matrix
Adj_sp <- getAdj(seu, platform = 'ST')
var.features <- seu@assays$RNA@var.features
X <- Matrix::t(seu[["RNA"]]@data[var.features,])
# maxIter = 2 is only used for illustration, and user can use default.
drscList <- DR.SC_fit(X,Adj_sp=Adj_sp ,K=4, maxIter=2, verbose=TRUE)
drsc1 <- selectModel(drscList)
str(drsc1)
}

