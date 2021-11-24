\name{RunWPCA}
\alias{RunWPCA}
\title{Run Weighted Principal Component Analysis}
\description{
  Run a weighted PCA dimensionality reduction
}
\usage{
  RunWPCA(object, q=15)
  ### S3 method for class "Seurat"
  ## RunWPCA(object, q=15)
  
  ### S3 method for class "matrix"
  ## RunWPCA(object, q=15)
  
  ### S3 method for class "dgCMatrix"
  ## RunWPCA(object, q=15)
}
\arguments{
  \item{object}{an object named "Seurat", "maxtrix" or "dgCMatrix". The object of class "Seurat" must include slot "scale.data".}
  \item{q}{an optional positive integer, specify the number of features to be extracted.}
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
  None
}
\examples{
\dontrun{
  library(Seurat)
  seu <- gendata_RNAExp(height=20, width=20,p=100, K=4)
  ## log-normalization
  seu <- NormalizeData(seu)
  ##
  seu <- FindVariableFeatures(seu, nfeatures=80)
  ## Scale
  seu <- ScaleData(seu)
  ## Run WPCA
  seu <- RunWPCA(seu)
  seu
  ## Run tSNE based on wpca
  seu <- RunTSNE(seu, reduction='wpca')
  seu
  ## Find SVGs
  seu <- FindSVGs(seu, nfeatures=80)
  (genes <- topSVGs(seu, ntop=10))
  Idents(seu) <- factor(paste0("cluster", seu$true_clusters), levels=paste0("cluster",1:4))
  RidgePlot(seu, features = genes[1:2], ncol = 2)
  FeaturePlot(seu, features = genes[1:2], reduction = 'tsne' ,ncol=2)
  }
}
\keyword{Dimension reduction}
\keyword{Weighted PCA}
