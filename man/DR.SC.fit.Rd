\name{DR.SC_fit}
\alias{DR.SC_fit}
\title{Joint dimension reduction and spatial clustering}
\description{
  Joint dimension reduction and spatial clustering for scRNA-seq and spatial transcriptomics data
}
\usage{
  DR.SC_fit(X, K, Adj_sp=NULL, q=15,
             error.heter= TRUE, beta_grid=seq(0.5, 5, by=0.5),
             maxIter=25, epsLogLik=1e-5, verbose=FALSE, maxIter_ICM=6,
             wpca.int=FALSE, coreNum = length(K))
}
\arguments{
 \item{X}{a  sparse matrix with class \code{dgCMatrix} or \code{matrix}, specify the log-normalization gene expression matrix used for DR-SC model.}
  \item{K}{a positive  integer allowing scalar or vector, specify the number of clusters in model fitting.}
 \item{Adj_sp}{an optional sparse matrix with class \code{dgCMatrix}, specify the adjoint matrix used for DR-SC model. We provide this interface for those users  who would like to define the adjoint matrix by their own.}
  \item{q}{a positive integer, specify the number of latent features to be extracted, default as 15.}
  \item{error.heter}{an optional logical value, whether use the heterogenous error for DR-SC model, default as \code{TRUE}. If \code{error.heter=FALSE}, then the homogenuous error is used for probabilistic PCA model in DR-SC.}
  \item{beta_grid}{an optional vector of positive value, the candidate set of the smoothing parameter to be searched by the grid-search optimization approach.}
   \item{maxIter}{an optional positive value, represents the maximum iterations of EM.}
  \item{epsLogLik}{an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.}
  \item{verbose}{an optional logical value, whether output the information of the ICM-EM algorithm.}
  \item{maxIter_ICM}{an optional positive value, represents the maximum iterations of ICM.}
  \item{wpca.int}{an optional logical value, means whether use the weighted PCA to obtain the initial values of loadings and other paramters, default as \code{FALSE} which means the ordinary PCA  is used.}
  \item{coreNum}{an optional positive integer, means the number of thread used in parallel computating.}
}
\details{
  Nothing
}
\value{
  DR.SC_fit returns a \link{list} with class "drscObject" with the following three components:
  \item{Objdrsc}{a list including the model fitting results, in which the number of elements is same as the length of K.}
  \item{out_param}{a numeric matrix used for model selection in MBIC.}
  \item{K_set}{a scalar or vector equal to input argument K.}
  
In addition, each element of "Objdrsc" is a list with the following comoponents:
\item{cluster}{inferred class labels}
\item{hZ}{extracted latent features.}
\item{beta}{estimated smoothing parameter}
\item{Mu}{mean vectors of mixtures components.}
\item{Sigma}{covariance matrix of mixtures components.}
\item{W}{estimated loading matrix}
\item{Lam_vec}{estimated variance of errors in probabilistic PCA model}
\item{loglik}{pseudo observed log-likelihood.}
}
\references{
\href{https://www.biorxiv.org/content/10.1101/2021.12.25.474153v4}{Wei Liu, Xu Liao, Yi Yang, Huazhen Lin, Joe Yeong, Xiang Zhou, Xingjie Shi and Jin Liu. (2022) Joint dimension reduction and clustering analysis for single-cell RNA-seq and spatial transcriptomics data}
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
seu <- gendata_RNAExp(height=10, width=10,p=50, K=4)
library(Seurat)
seu <- NormalizeData(seu, verbose=FALSE)
# choose 40 highly variable features using FindVariableFeatures in Seurat
# seu <- FindVariableFeatures(seu, nfeatures = 40)
# or choose 40 spatailly variable features using FindSVGs in DR.SC
seu <- FindSVGs(seu, nfeatures = 40, verbose=FALSE)
# users define the adjacency matrix
Adj_sp <- getAdj(seu, platform = 'ST')
var.features <- seu@assays$RNA@var.features
X <- Matrix::t(seu[["RNA"]]@data[var.features,])
# maxIter = 2 is only used for illustration, and user can use default.
drscList <- DR.SC_fit(X,Adj_sp=Adj_sp, K=4, maxIter=2, verbose=TRUE)

}
\keyword{Dimension reduction}
\concept{spatial clustering and cell type identification}
