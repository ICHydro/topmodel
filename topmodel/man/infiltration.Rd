\name{infiltration}
\alias{infiltration}
\title{Infiltration routine of TOPMODEL}
\description{
  Infiltration routine used in TOPMODEL, based on the Green-Ampt model.
}
\usage{
  infiltration(rain, parameters)
}
\arguments{
  \item{rain}{A vector of rain data (m per timestep)}
  \item{parameters}{A vector containing 3 parameters (see below for the exact structure)}
}
\details{
This function gives direct access to the infiltration routine implemented in \code{topmodel()}.
The function needs three parameters: \code{c(dt, CD, K0, m)}, where:

\tabular{ll}{
dt 	\tab The timestep (hours)\cr
CD 	\tab capillary drive, see Morel-Seytoux and Khanji (1974)\cr
K0	\tab Surface hydraulic conductivity (m/h)\cr
m 	\tab Model parameter controlling the rate of decline of transmissivity in the soil profile, see Beven, 1984\cr
}

}
\value{
The function returns a vector with the same length as the input vector rain representing infiltration.
}
\references{Morel-Seytoux, H.J., Khanji, J., 1974. Derivation of an Equation of Infiltration. Water Resources Research, 10, 795-800.
Beven, K., 1984. Infiltration into a Class of Vertically Non-Uniform Soils. Hydrological Sciences Journal 29, 425-434.

See also \url{https://github.com/ICHydro/topmodel} for examples.
}
\author{Wouter Buytaert, Imperial College London, based on original FORTRAN code from Lancaster University}
%\note{}

\seealso{\code{\link{topmodel}}}
%\examples{}
\keyword{models}
