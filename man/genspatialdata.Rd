\name{gendata_RNAExp}
\alias{gendata_RNAExp}
\title{Generate simulated data}
\description{
  Generate simulated spatial transcriptomics data or scRNAseq data.
}
\usage{
  gendata_RNAExp(height=30, width=30, platform="ST", p =100, q=10, K=7, 
                            G=4,sigma2=1, tau=8, seed=1, view=FALSE)
}
\arguments{
 \item{height, width}{Height and width of lattice grids for generating spatial coordinates. n=height * width cells for scRNAseq data}
 \item{platform}{Set the platform for the simulated data.}
 \item{p}{number of genes to generate.}
 \item{q}{number of true latent features to generate gene expression}
 \item{K}{number of clusters (cell types).}
 \item{seed}{random seed for generate data}
 \item{G}{the number of neighbors. The latter must be one of G = 4 or G = 8, which respectively correspond to a first order and a second order dependency structure. By default, G = 4.}
 \item{sigma2}{Variance of error term in probabilitic PCA model.}
 \item{tau}{a positive factor of mixture mean values.}
 \item{view}{Logical value indicating whether the draw should be printed. Do not display the optional borders.}
}
\details{
  Nothing
}
\value{
 return a "Seurat" object. If \code{platform="ST"}, then the metadata of this Seurat object will include two columns with names "row" and "col" which are the spatial coordinates; If \code{platform="scRNAseq"}, then the metadata of this Seurat object will not have them.
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
seu <- gendata_RNAExp(height=20, width=20,p=200, K=4)
seu

## generate scRNAseq data
seu <- gendata_RNAExp(height=20, width=20, platform="scRNAseq", p=100, K=4)
seu
}
