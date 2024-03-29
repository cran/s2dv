% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InsertDim.R
\name{InsertDim}
\alias{InsertDim}
\title{Add a named dimension to an array}
\usage{
InsertDim(data, posdim, lendim, name = NULL)
}
\arguments{
\item{data}{An array to which the additional dimension to be added.}

\item{posdim}{An integer indicating the position of the new dimension.}

\item{lendim}{An integer indicating the length of the new dimension.}

\item{name}{A character string indicating the name for the new dimension. 
The default value is NULL.}
}
\value{
An array as parameter 'data' but with the added named dimension.
}
\description{
Insert an extra dimension into an array at position 'posdim' with length 
'lendim'. The array repeats along the new dimension.
}
\examples{
a <- array(rnorm(15), dim = c(a = 3, b = 1, c = 5, d = 1))
res <- InsertDim(InsertDim(a, posdim = 2, lendim = 1, name = 'e'), 4, c(f = 2))
dim(res)

}
