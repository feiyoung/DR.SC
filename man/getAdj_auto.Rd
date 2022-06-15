\name{getAdj_auto}
\alias{getAdj_auto}
\title{Calculate adjacency matrix by automatically choosing radius}
\usage{
getAdj_auto(pos, lower.med=4, upper.med=6, radius.upper= 100)
}
\arguments{
\item{pos}{a n-by-2 matrix of position.}
\item{lower.med}{an integer, the lower bound of median number of neighbors among all spots.}
\item{upper.med}{an integer, the upper bound of median number of neighbors among all spots.}
\item{radius.upper}{a real, the upper bound of radius.}
}
\value{
A sparse adjacency matrix containing the neighbourhood.
}
\seealso{
  \link{getAdj}, \link{getAdj_manual}.
}
\description{
an efficient function to find the radius by bi-section method then find neighborhoods based on the matrix of position, which ensures that each spot has about lower.med~upper.med neighborhoods in the sense of median.
}
