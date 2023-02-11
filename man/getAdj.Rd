\name{getAdj}
\alias{getAdj}
\alias{getAdj.Seurat}
\title{Calculate the adjacency matrix given the spatial coordinates}
\description{
  Calculate the adjacency matrix for the spatial transcriptomics data measured on 10X Visium, ST or other platforms as a Seurat object.  
}
\usage{

  \method{getAdj}{Seurat}(obj, platform = c('Visium', "ST", "Other_SRT"), ...) 
  
}
\arguments{
  \item{obj}{an object with class "Seurat", there are spatial coordinates information in the metadata of obj, named "row" and "col", where first column is x-axis coordinate, the second column is y-axis coordinate. \link{getAdj_manual} and \link{getAdj_auto}  supports multi-dimensional spatial coordinates with a matrix as input.}
  \item{platform}{a string, specify the platform of the provided data, default as "Visium". There are more platforms to be chosen, including ("Visuim", "ST", "Other_SRT"), which means there are spatial coordinates information in the metadata of obj, named "row" and "col". The platform helps to calculate the adjacency matrix by defining the neighborhoods.}
  \item{...}{Other arguments to pass into \link{getAdj_auto} function.}
}
\details{
  For lattice grids in ST platform (ST), the interior spot has four neighbors (left, right, up and down),the boundary spot has three neighbors, and the spot in the corner has two neighbors. For hexagon grids, such as spatial coordinate in 10X Visium platform (Visium), the interior spot has six neighbors. For the irregular coordinates in other platforms (Other_SRT), Euclidean distance is adopted to decide whether a spot is an neighbor of another spot. For example, if the Euclidean distance between spot A and B is less than a radius, then A is taken as the neighbourhood of B. See functions \link{getAdj_auto} and \link{getAdj_manual} for more details.
}
\value{
   Return a \code{dgCMatrix} object recording the information of neighborhoods about each spot.
}
\references{
  \href{https://academic.oup.com/nar/article/50/12/e72/6555431}{Wei Liu, Xu Liao, Yi Yang, Huazhen Lin, Joe Yeong, Xiang Zhou, Xingjie Shi & Jin Liu (2022). Joint dimension reduction and clustering analysis of single-cell RNA-seq and spatial transcriptomics data, Nucleic Acids Research, gkac219.}
}
\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  \link{getAdj_auto}, \link{getAdj_manual}.
}
\examples{
  ## S3 method for class "Seurat"
  seu <- gendata_RNAExp(height=20, width=20,p=200, K=4)
  Adj_sp <- getAdj(seu, platform = 'ST')
}
\concept{Lattice grids}
\concept{Hexagon grids}
\concept{Adjacency matrix}
