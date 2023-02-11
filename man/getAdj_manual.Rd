\name{getAdj_manual}
\alias{getAdj_manual}
\title{Calculate adjacency matrix by user-specified radius}
\usage{
getAdj_manual(pos, radius)
}
\arguments{
\item{pos}{is a n-by-d matrix of position, where n is the number of spots, and d is the dimension of coordinates.}

\item{radius}{is a threashold of Euclidean distance to decide whether a spot is an neighborhood of another spot. For example, if the Euclidean distance between spot A and B is less than radius, then A is taken as the neighbourhood of B.}
}
\value{
A sparse matrix containing the neighbourhood
}
\seealso{
  \link{getAdj_auto}, \link{getAdj}.
}
\description{
an efficient function to find the neighbors based on the matrix of position and a pre-defined radius.
}
