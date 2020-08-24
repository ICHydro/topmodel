\name{flowlength}
\alias{flowlength}
\title{Calculates the flow distance towards the catchment outlet}
\description{
Given a digital elevation model and the coordinates of the catchment outlet,
this function calculates the flowlength of each gridcell to the outlet, based
on a single flow direction algorithm (D8).
}
\usage{
flowlength(DEM,outlet)
}
\arguments{
  \item{DEM}{A matrix representing a digital elevation model [m] with equally sized pixels and equal NS and EW resolution}
  \item{outlet}{A vector containing the row and column indices of the pixel
  representing the catchment outlet.}
}
\details{
The function returns the flowlength in cell size units. So you have to multiply by the map resolution to get the flow in meters.
}

\value{A matrix of the same size as \code{DEM}}

\references{See \url{http://paramo.cc.ic.ac.uk/topmodel_tutorial} for examples.}
\author{Wouter Buytaert, Imperial College London, based on an implementation from the Hydrology Group of Lancaster University}
%\note{}
\seealso{sinkfill, river}
%\examples{}
