\name{make.classes}
\alias{make.classes}
\title{make topographic index classes from a topographic index map}
\description{This function splits a dataset in n evenly distributed classes and calculates the number of elements of each class. It is very similar to \code{hist()}, but \code{hist()} does not always keep the number of breaks requested}
\usage{make.classes(array,n)}
\arguments{
  \item{array}{A numeric vector or matrix}
  \item{n}{Number of breaks}
}
%\details{}
\value{The function returns a data matrix with 2 columns. Column \code{breaks} contains the upper limits of each class, while \code{counts} contains the number of elements}
%\references{}
\author{Wouter Buytaert, Imperial College London}
%\note{}
\seealso{\code{\link{topidx}}, \code{\link{topmodel}}}
%\examples{}
\keyword{models}
