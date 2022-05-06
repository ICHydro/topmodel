\name{topidx}
\alias{topidx}
\title{Calculation of the topographic index from a DEM raster}
\description{Calculation of the topographic index from a DEM raster}
\usage{
topidx(DEM, resolution, river = NULL)
}
\arguments{
  \item{DEM}{A matrix representing a digital elevation model with equally sized pixels and equal NS and EW resolution}
  \item{resolution}{Resolution of the digital elevation model (m)}
  \item{river}{A matrix representing a river map (optional). All cells with values higher than 0 will be treated as river cells and are not included in the calculations.}
}
\details{The river map should not contain negative values}
\value{The function returns a list, containing 2 rasters. Atb are the topographic index values of each grid cell. Area contains the contributing area}
\references{See \url{https://github.com/ICHydro/topmodel} for examples.}
\author{Wouter Buytaert, Imperial College London, based on routines developed by the hydrology group of Lancaster University}
%\note{}

\seealso{\code{\link{topmodel}}}
\examples{
  data(huagrahuma.dem)
  topidx <- topidx(huagrahuma.dem, resolution= 25)$atb
  image(topidx)
}
\keyword{models}
