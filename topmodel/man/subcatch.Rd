\name{subcatch}
\alias{subcatch}
\title{Identify a hydrological catchment based on a single direction flow algorithm}
\description{Identify a hydrological catchment based on a single direction flow algorithm}
\usage{
  subcatch(DEM,outlet)
}
\arguments{
  \item{DEM}{A matrix representing a digital elevation model [m] with equally sized pixels and equal x and y resolution}
  \item{outlet}{A vector containing the row and column indices of the pixel
  representing the catchment outlet.}
}
%\details{}
\value{A binary matrix of the same size as \code{DEM}.}
\references{See \url{https://github.com/ICHydro/topmodel} for examples.}

\author{Wouter Buytaert, Imperial College London, based on an implementation from the Hydrology Group of Lancaster University}
%\note{}

\seealso{\code{\link{topmodel}}}
%\examples{}
%\keyword{}
