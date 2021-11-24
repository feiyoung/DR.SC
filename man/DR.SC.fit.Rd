\name{DR.SC_fit}
\alias{DR.SC_fit}
\title{Joint dimension reduction and spatial clustering}
\description{
  Joint dimension reduction and spatial clustering for scRNA-seq and spatial transcriptomics data
}
\usage{
  DR.SC_fit(X,Adj_sp=NULL, q=15, K= NULL,error.heter= TRUE, K_set = seq(2, 10),
    beta_grid=seq(0.5, 5, by=0.5),maxIter=25, epsLogLik=1e-5, verbose=FALSE, 
    maxIter_ICM=6,pen.const=1,wpca.int=FALSE, parallel='parallel', num_core=5)
}
\arguments{
 \item{X}{an optional sparse matrix with class \code{dgCMatrix}, specify the log-normalization gene expression matrix used for DR-SC model.}
 \item{Adj_sp}{an optional sparse matrix with class \code{dgCMatrix}, specify the adjoint matrix used for DR-SC model. We provide this interface for those users  who would like to define the adjoint matrix by their own.}
  \item{q}{a positive integer, specify the number of latent features to be extracted, default as 15.}
  \item{K}{a positive  integer, specify the number of clusters, default as \code{NULL}. When \code{K=NULL}, it is automatically selected by MBIC criteria.}
  \item{K_set}{a vector of positive integer, means the candidates of number of clusters used for MBIC.}
  \item{error.heter}{an optional logical value, whether use the heterogenous error for DR-SC model, default as \code{TRUE}. If \code{error.heter=FALSE}, then the homogenuous error is used for model.}
  \item{beta_grid}{an optional vector of positive value, the candidate set of the smoothing parameter to be searched by the grid-search optimization approach.}
   \item{maxIter}{an optional positive value, represents the maximum iterations of EM.}
  \item{epsLogLik}{an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.}
  \item{verbose}{an optional logical value, whether output the information of the ICM-EM algorithm.}
  \item{maxIter_ICM}{an optional positive value, represents the maximum iterations of ICM.}
  \item{pen.const}{an optional positive value, the adjusted constant used in the MBIC criteria. It usually takes value between 0.1 to 1.}
  \item{wpca.int}{an optional logical value, means whether use the weighted PCA to obtain the initial values of loadings and other paramters, default as \code{FALSE} which means the conventional PCA  is used.}
  \item{parallel}{a optional string, specify the parallel way to choose the number of clusters by MBIC. We provide two methods: 1. \code{parallel="parrallel"} uses \code{parrallel} R package to conduct the parallel schema; 2.\code{parallel=NULL} doesn't use parallel computation.}
  \item{num_core}{an optional positive integer, means the cores used in parallel computating.}
}
\details{
  Nothing
}
\value{
  DR.SC_fit returns a \link{list} with the following components:
\item{cluster}{inferred class labels}
\item{hZ}{extracted latent features.}
\item{beta}{estimated smoothing parameter}
\item{Mu}{mean vectors of mixtures components.}
\item{Sigma}{covariance matrix of mixtures components.}
\item{W}{estimated loading matrix}
\item{Lam_vec}{estimated variance of errors in probabilistic PCA model}
\item{loglik}{pseudo observed log-likelihood}
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
seu <- gendata_RNAExp(height=10, width=10,p=50, K=4)
library(Seurat)
seu <- NormalizeData(seu)
# choose 2000 variable features using Seurat
seu <- FindVariableFeatures(seu, nfeatures = 40)
# users define the adjoint matrix
Adj_sp <- getAdj(seu, platform = 'ST')
var.features <- seu@assays$RNA@var.features
X <- Matrix::t(LogNormalize(seu@assays$RNA@counts[var.features,],))
# maxIter = 2 is only used for illustration, and user can use default.
drscList <- DR.SC_fit(X,Adj_sp=Adj_sp ,K=4, maxIter=2, verbose=TRUE)

}
\keyword{Dimension reduction}
\concept{spatial clustering and cell type identification}
