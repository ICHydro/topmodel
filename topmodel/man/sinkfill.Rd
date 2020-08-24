\name{sinkfill}
\alias{sinkfill}
\title{Fill sinks in a digital elevation model}
\description{Removes sinks in a digital elevation model by filling depressions}
\usage{
  sinkfill(DEM,res,degree) 
}
\arguments{
  \item{DEM}{A matrix representing a digital elevation model [m] with equally sized pixels and equal x and y resolution}
  \item{res}{Resolution of the digital elevation model [m]}
  \item{degree}{Minimum slope to be kept between cells when filling [degrees]}
}
\details{
For deep sinks or large maps, it may be possible that not all sinks are filled in one run. Then the function should be applied repeatedly over the same object.
}
\value{A matrix of the same size as the DEM.}

\references{See \url{http://paramo.cc.ic.ac.uk/topmodel_tutorial} for examples.}

\author{Wouter Buytaert, Imperial College London, based on an implementation from the Hydrology Group of Lancaster University}

%\note{}

\seealso{topmodel, topidx}

\examples{
data(huagrahuma.dem)
filled.dem <- sinkfill(huagrahuma.dem, 25, 0.1)}

%\keyword{}
