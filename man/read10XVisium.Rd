\name{read10XVisium}
\alias{read10XVisium}
\title{Read the spatial transcriptomics data measured on 10X Visium platform}
\description{
  Read the spatial transcriptomics data measured on 10X Visium platform  as a Seurat object, where the spatial coordinates are saved  in the metadata, named "row" and "col".
}
\usage{
  read10XVisium(dirname)
}
\arguments{
  \item{dirname}{A string, the dictory of Visium datasets}
}
\details{
  Nothing
}
\value{
  return a Seurat object.
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
\dontrun{
  ## set your file directory, then read it.
  data_name <- "D/HCC"
  HCC1 <- read10XVisium(data_name)
  }
}
\keyword{Input of spatial data}
