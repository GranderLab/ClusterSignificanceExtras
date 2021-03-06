% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrixTools.R
\name{visualizeOverlaps}
\alias{visualizeOverlaps}
\title{Plot Overlaps}
\arguments{
\item{mat}{The matrix to be checked.}

\item{groups}{A character vector indicating the groups location within the
matrix. Can be generated with the makeGroups function.}
}
\description{
Plots the group overlaps in 2 dimensions for a given matrix and groups.
}
\examples{

x1 <- 0
x2 <- x1+1
mat <- makeMatrix(x=c(x1,x2), y=c(0,1), dim=2, points=10)
groups <- makeGroups(mat, names=c("grp1", "grp2"))

visualizeOverlaps(mat, groups)


}
\author{
Jason Serviss
}
\keyword{visualizeOverlaps}
