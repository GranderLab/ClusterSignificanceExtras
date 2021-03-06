% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{zero}
\alias{zero}
\title{Zero overlap dataset.}
\format{list of arrays
\describe{
\item{list names}{The number of points per group.}
\item{list elements}{An array containing all repititions of the
   ClusterSignificance input matrix for the specific points number indicated
   in the list element name}
\item{array dimension 1}{Rows of the input matrix to ClusterSignificance.}
\item{array dimension 2}{Columns of the input matrix to ClusterSignificance.}
\item{array dimension 3}{The repititions for each points per group amounts.}
}}
\usage{
zero
}
\value{
list of arrays
}
\description{
Output from the zeroOrHundred function with overlap argument = 0 and other
default arguments.
}
\examples{
zero

}
