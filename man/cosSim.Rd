% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cosSim.R
\name{cosSim}
\alias{cosSim}
\alias{cosSim.default}
\alias{cosSim.matrix}
\alias{cosSim.data.frame}
\title{Cosine similarity}
\usage{
cosSim(x, ...)

\method{cosSim}{default}(x, y)

\method{cosSim}{matrix}(x, y, all.vs.all = FALSE, by.row = TRUE)

\method{cosSim}{data.frame}(...)
}
\arguments{
\item{x}{A numeric vector or matrix}

\item{y}{A numeric vector or matrix}

\item{all.vs.all}{(When x and y are matrices) If TRUE, all rows of x are compared with all rows
of y. If FALSE, rows of x are compared with rows of y 'side-by-side'}

\item{by.row}{(When x and y are matrices) If TRUE, compare rows. If FALSE, compare columns}
}
\value{
A numeric vector or matrix
}
\description{
Calculate cosine similarity between two vectors of the same length, or cosine
similarity between two matrices
}
\examples{
set.seed(1)
x <- matrix(runif(200), nrow=20)
y <- matrix(runif(200), nrow=20)

## Works with matrices
cosSim(x, y)

## ...or vectors
cosSim(x[1,], y[1,])

}
