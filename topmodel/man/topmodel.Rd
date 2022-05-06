\name{topmodel}
\alias{topmodel}
\title{Implementation of the hydrological model TOPMODEL}
\description{
  Implementation of the 1995 Fortran version of TOPMODEL by Keith Beven.
}
\usage{
topmodel(parameters, topidx, delay, rain, ETp, verbose = F, Qobs=NA)
}
\arguments{
  \item{parameters}{A vector or a matrix containing the input parameters (see below for the exact structure)}
  \item{topidx}{A 2 column matrix with respectively the topographic index classes and values (see below for the exact structure)}
  \item{delay}{Delay function for overland flow (see below)}
  \item{rain}{A vector of rain data (m per timestep)}
  \item{ETp}{A vector of potential evapotranspiration data (m per timestep)}
  \item{verbose}{If set to TRUE, returns in addition to predicted discharge also overland flow, base flow and saturated zone storage}
  \item{Qobs}{If Qobs is given, normal output is suppressed and only the Nash and Sutcliffe efficiency is returned (m per timestep)}
}
\details{
topmodel() automatically implements a Monte Carlo simulation. If the parameter argument is a single vector, the model is run once. If the parameter argument is a matrix, each row should represent a parameter set. In that case, the model is run with each parameter set (see the examples below).

A single parameter set consists of: \code{c(qs0,lnTe,m,Sr0,SrMax,td,vch,vr,k0,CD,dt)}, with:

\tabular{ll}{
qs0 \tab Initial subsurface flow per unit area [m]\cr
lnTe 	\tab log of the areal average of T0 [m2/h]\cr
m	\tab Model parameter controlling the rate of decline of transmissivity in the soil profile, see Beven, 1984\cr
Sr0 	\tab Initial root zone storage deficit [m]\cr
Srmax 	\tab Maximum root zone storage deficit [m]\cr
td 	\tab Unsaturated zone time delay per unit storage deficit [h/m]\cr
vch 	\tab channel flow outside the catchment [m/h] (currently not used)\cr
vr 	\tab channel flow inside catchment [m/h]\cr
k0 	\tab Surface hydraulic conductivity [m/h]\cr
CD 	\tab capillary drive, see Morel-Seytoux and Khanji (1974)\cr
dt 	\tab The timestep [h]\cr
}

The topidx dataframe can be derived conveniently with \code{make.classes()}. It should contain 2 columns. The first column should give the lower boundary of each topographic index class, and the second column should give the respective area fraction. The second column must sum to 1.

k0 and CD are used only for the infiltration excess routine. Set k0 to a very high value to avoid infiltration excess overland flow

Flow is routed through a delay function which represents the time spent in the channel system. The parameter \code{delay} is used for this. Delay is a matrix with 2 columns. The first column gives the cumulative relative area. The second column gives the average distance towards the outlet (m).

}
\value{
The function returns an array of observed discharges. If more than one parameter set is given, a matrix is returned, with each column representing a discharge set coinciding with the parameter sets. If Qobs is given, the function returns an array of Nash-Sutcliffe efficiencies, 1 for each parameter sets.

If verbose output is requested, a list is returned, with the modelled discharge (Q), overland flow (qo), subsurface flow (qs), storage deficit (S), infiltration excess overland flow (fex), and actual evapotranspiration (Ea) for each time step.

Be aware that invoking \code{topmodel()} without Q for a large number of runs, may require a large amount of memory.
}
\references{Beven, K. J., Kirkby, M. J., 1979. A physically based variable contributing area model of basin hydrology. Hydrol. Sci. Bull. 24, 43-69.

Beven K, Lamb R, Quinn P, Romanowicz R, Freer J, 1995. TOPMODEL. In: Sing VP (Ed), Computer Models of Watershed Hydrology. Water Resources Publications, Colorado. pp. 627-668.

Morel-Seytoux, H.J., Khanji, J., 1974. Derivation of an Equation of Infiltration. Water Resources Research, 10, 795-800.

Beven, K., 1984. Infiltration into a Class of Vertically Non-Uniform Soils. Hydrological Sciences Journal 29, 425-434.

See also \url{https://github.com/ICHydro/topmodel} for a more examples on how to run topmodel in R.
}
\author{Wouter Buytaert, Imperial College London}
% \note{}

\seealso{\code{\link{topidx}}}
\examples{
data(huagrahuma)
attach(huagrahuma)

## returns the simulated runoff (Qobs not given)
Qsim <- topmodel(parameters, topidx, delay, rain, ETp)

## returns a list of simulated runoff (Q), overland flow (qo), subsurface flow (qs) and storage (S):
Qsim <- topmodel(parameters, topidx, delay, rain,ETp, verbose = TRUE)

## plot observed and simulated discharge:
plot(Qobs)
points(Qsim$Q, col="red", type="l")

## For a Monte carlo sampling from a uniform distribution, we construct a parameter matrix:

runs<-10

qs0   <- runif(runs, min=0, max=4e-5)
lnTe  <- runif(runs, min=-2, max=1)
m     <- runif(runs, min=0, max=0.2)
Sr0   <- runif(runs, min=0, max=0.02)
Srmax <- runif(runs, min=0, max=2)
td    <- runif(runs, min=0, max=3)
vch   <- 1000
vr    <- runif(runs, min=100, max=2500)
k0    <- runif(runs, min=0, max=0.01)
CD    <- runif(runs, min=0, max=5)
dt    <- 0.25

parameters<-cbind(qs0,lnTe,m,Sr0,Srmax,td,vch,vr,k0,CD,dt)

## returns an array of 10 Nash Sutcliffe efficiencies; one for each parameter set:
result<-topmodel(parameters, topidx, delay, rain, ETp, Qobs = Qobs)

}
\keyword{models}
