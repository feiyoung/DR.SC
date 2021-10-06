\name{gendata_RNAExp}
\alias{gendata_RNAExp}
\title{Generate simulated data}
\description{
  Generate simulated spatial transcriptomics data or scRNAseq data.
}
\usage{
  gendata_RNAExp(height=30, width=30, platform="ST", p =100, q=10, K=7, 
                            G=4,sigma2=1, tau=8, seed=1, view=F)
}
\arguments{
 \item{height, width}{Height and width of lattice grids for generating spatial coordinates. n=height * width cells for scRNAseq data}
 \item{platform}{Set the platform for the simulated data.}
 \item{p}{number of genes to generate.}
 \item{q}{number of true latent features to generate gene expression}
 \item{K}{number of clusters (cell types).}
 \item{seed}{random seed for generate data}
}
\details{
  Nothing
}
\value{
 return a "Seurat" object. If \code{platform="ST"}, then the metadata of this Seurat object will include two columns with names "row" and "col" which are the spatial coordinates; If \code{platform="scRNAseq"}, then the metadata of this Seurat object will not have them.
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

}
\examples{
## we generate the spatial transcriptomics data with lattice neighborhood, i.e. ST platform.
seu <- gendata_RNAExp(height=30, width=30,p=500, K=4)
seu

## generate scRNAseq data
seu <- gendata_RNAExp(height=30, width=30, platform="scRNAseq", p=500, K=4)
seu

}
\keyword{Simulated data.}
