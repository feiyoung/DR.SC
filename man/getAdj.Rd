\name{getAdj}
\alias{getAdj}
\title{Calculate the Adjoint matrix given the spatial coordinates}
\description{
  Calculate the Adjoint matrix for the spatial transcriptomics data measured on 10X Visium or ST platform as a Seurat object.  
}
\usage{

  getAdj(object, ...) 
  ### S3 method for class "Seurat"
  #getAdj(seu, platform = c('Visium', 'ST'))
  
  ### S3 method for class "matrix"
  #getAdj(pos, radius)
}
\arguments{
  \item{seu}{an object of class "Seurat".}
  \item{platform}{a string, specify the platform of the provided data, default as "Visium". There are three platforms to be chosen, including "Visuim", "ST" and "scRNAseq", where the first two means there are spatial coordinates information in the metadata of seu, named "row" and "col" , and the last one "scRNAseq" means there is no spatial information in object seu. The platform helps to calculate the adjoint matrix by defining the neighborhoods.}
  \item{pos}{a matrix with d columns, typical one is the two columns where first column is x-axis coordinate, the second column is y-axis coordinate. This function support multi-dimensional spatial coordinates.}
  \item{radius}{a positive real, speficify the radius of Euclidean distance between two spots to define the neighbors. The radius is larger, then the neighbors will be more.}
}
\details{
  For lattice grids, i.e., two-dimensional coordinates, the interior spot has four neighbors (left, right, up and down),the boundary spot has three neighbors, and the spot in the corner has two neighbors. For hexagon grids, such as spatial coordinate in 10X Visium platform, the interior spot has six neighbors. More flexible definition can be used if there are some additional information. And then use \code{getAdj.matrix} function to evalute the Adjoint matrix as a input for \link{DR.SC.fit} to run DR-SC model.
}
\value{
   Return a \code{dgCMatrix} object recording the information of neighborhoods about each spot.
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
  ## S3 method for class "Seurat"
  seu <- gendata_RNAExp(height=30, width=30,p=500, K=4)
  Adj_sp <- getAdj(seu, platform = 'ST')
  
  ## S3 method for class "matrix"
  pos <- cbind(seu$row, seu$col)
  Adj_sp <- getAdj(pos, radius=2)
  cat("The number of neighborhoods for each spot is about: ", summary(apply(Adj_sp, 1, function(x) sum(x==1)))[4], '\n')
  ## To ensure thatST platform has 4 neighborhoods, we tune the 'radius' parameters to 1.4
  Adj_sp <- getAdj(pos, radius=1.4)
  cat("The number of neighborhoods for each spot is about: ", summary(apply(Adj_sp, 1, function(x) sum(x==1)))[4], '\n')
}
\keyword{Lattice grids, Hexagon grids, Adjoint matrix.}
