\name{getAdj_auto}
\alias{getAdj_auto}
\title{Calculate adjacency matrix by automatically choosing radius}
\usage{
getAdj_auto(pos)
}
\arguments{
\item{pos}{is a n-by-2 matrix of position.}
}
\value{
A sparse adjoint matrix containing the neighbourhood.
}
\seealso{
  \link{getAdj}, \link{getAdj_manual}.
}
\description{
an efficient function to find the radius by bi-section method then find neighborhoods based on the matrix of position and the chosen radius, which ensures that each spot has about 4~6 neighborhoods in the sense of median and shows it works well.
}
