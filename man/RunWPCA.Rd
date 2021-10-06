\name{RunWPCA}
\alias{RunWPCA}
\title{Run Weighted Principal Component Analysis}
\description{
  Run a weighted PCA dimensionality reduction
}
\usage{
  RunWPCA(object, ...)
  ### S3 method for class "Seurat"
  #RunWPCA(seu, q=15)
  
  ### S3 method for class "matrix"
  #RunWPCA(X, q=15)
  
  ### S3 method for class "dgCMatrix"
  #RunWPCA(X, q=15)
}
\arguments{
  \item{seu}{an object of class "Seurat" with slot "scale.data".}
  \item{q}{an optional positive integer, specify the number of features to be extracted.}
  \item{X}{a object of class "matrix" or "dgCMatrix" whose row is gene and column is spot/cell.}
}
\details{
  Nothing
}
\value{
  For Seurat object, return a Seurat object. For objcet "matrix" and "dgCMatrix", return a object "matrix" with rownames same as the colnames of \code{X}, and colnames "WPCA1" to "WPCAq".
}
\references{
  Bai, J. and Liao, Y. (2017). Inferences in panel data with interactive effects using large
covariance matrices. Journal of Econometrics, 200(1):59–78.
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
  library(Seurat)
  seu <- gendata_RNAExp(height=40, width=40,p=1000, K=4)
  ## log-normalization
  seu <- NormalizeData(seu)
  ##
  seu <- FindVariableFeatures(seu, nfeatures=800)
  ## Scale
  seu <- ScaleData(seu)
  ## Run WPCA
  seu <- RunWPCA(seu)
  seu
  ## Run tSNE based on wpca
  seu <- RunTSNE(seu, reduction='wpca')
  seu
  ## Find SVGs
  seu <- FindSVGs(seu, nfeatures=800)
  (genes <- topSVGs(seu, ntop=10))
  Idents(seu) <- factor(paste0("cluster", seu$true_clusters), levels=paste0("cluster",1:4))
  RidgePlot(seu, features = genes[1:2], ncol = 2)
  FeaturePlot(seu, features = genes[1:2], reduction = 'tsne' ,ncol=2)
}
\keyword{Dimension reduction, WPCA.}
