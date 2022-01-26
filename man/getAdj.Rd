\name{getAdj}
\alias{getAdj}
\alias{getAdj.Seurat}
\title{Calculate the Adjoint matrix given the spatial coordinates}
\description{
  Calculate the Adjoint matrix for the spatial transcriptomics data measured on 10X Visium or other platforms as a Seurat object.  
}
\usage{

  getAdj(obj, platform ='Visium')
  ## getAdj(obj, platform = c('Visium', 'ST', "seqfish", 
  ## 'merfish', 'slide-seqv2', 'seqscope'))
  \method{getAdj}{Seurat}(obj,platform ='Visium')
  
}
\arguments{
  \item{obj}{an object with class "Seurat", there are spatial coordinates information in the metadata of obj, named "row" and "col", where first column is x-axis coordinate, the second column is y-axis coordinate. \link{getAdj_manual} and \link{getAdj_auto}  supports multi-dimensional spatial coordinates with a matrix as input.}
  \item{platform}{a string, specify the platform of the provided data, default as "Visium". There are many platforms to be supported, including ("Visuim", "ST", "seqfish", 'merfish', 'slide-seqv2', 'seqscope'), which means there are spatial coordinates information in the metadata of seu, named "row" and "col". The platform helps to calculate the adjoint matrix by defining the neighborhoods.}
}
\details{
  For lattice grids, i.e., two-dimensional coordinates, the interior spot has four neighbors (left, right, up and down),the boundary spot has three neighbors, and the spot in the corner has two neighbors. For hexagon grids, such as spatial coordinate in 10X Visium platform, the interior spot has six neighbors. More flexible definition can be used if there are some additional information. And then use \link{getAdj_manual} function to evalute the Adjoint matrix as a input for \link{DR.SC_fit} to run DR-SC model.
}
\value{
   Return a \code{dgCMatrix} object recording the information of neighborhoods about each spot.
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
  \link{getAdj_auto}, \link{getAdj_manual}, \link{getAdj_reg}.
}
\examples{
  ## S3 method for class "Seurat"
  seu <- gendata_RNAExp(height=20, width=20,p=200, K=4)
  Adj_sp <- getAdj(seu, platform = 'ST')
}
\concept{Lattice grids}
\concept{Hexagon grids}
\concept{Adjoint matrix}
