
## install and load the required packages:

install.packages("topmodel")
install.packages("Hmisc")
library(topmodel)
library(Hmisc)

############ PART 0: topographical analysis ##############

# Topmodel requires to types of topographical information: the topographic index distribution, and a channel flow delay function
# As topmodel is a semidistributed model, it does not use the maps directly, but instead uses histograms.
# Many GIS packages include routines to calculate these maps, but you can do it directly in R as well,
# starting from a digital elevation model (DEM).

# The DEM has to be imported as a matrix, which can then be processed by topidx().
# Take for instance this minimalistic DEM, saved in a test file called "DEM.txt"
# Values outside the catchment are given the value -9999 (this can be any other value that, obviously, does not occur in the DEM values):

-9999 -9999 828.9 835.6 -9999
818.3 826.0 830.7 834.5 836.0
817.1 824.0 825.2 833.3 836.9
816.5 820.0 824.1 330.8 -9999
810.7 815.6 822.2 -9999 -9999

# This file can be imported and processed in R with:

DEM <- read.table("c:/DEM.txt")
DEM <- as.matrix(DEM)

# Remove the values outside the catchment:

DEM[DEM==-9999] <- NA

# You may want to plot the DEM to see whether everything looks OK:

image(DEM)

# Then calculate the topographic index, the resolution should be in [m].
# Here we use the DEM from Huagrahuma as an example:

data(huagrahuma.dem)
DEM <- sinkfill(huagrahuma.dem, cellsize=25, degree=0.1)
topindex <- topidx(DEM, resolution=25)

# The values need to be split into a set of classes, since topmodel() is a semidistributed model that lumps hydrologically similar areas into the same hydrological response units.
# Here we define 16 hydrological response units:

topidx <- make.classes(topindex,16)

# the delay function is a bit more tricky because this requires cumulative fractions, but you generate it as follows:

n <- 5 # number of classes; a higher number will result in a smoother histogram
delay <- flowlength(huagrahuma.dem)*25 # TODO: add the outlet coordinates; otherwise the flowlength will be calculated to the edge of the map.
delay <- make.classes(delay, n)
delay <- delay[n:1,]
delay[,2] <- c(0, cumsum(delay[1:(n-1),2]))

############ PART 1: running the rainfall-runoff model ##############

## Load the example dataset from the Huagrahuma catchment
## and attach it to the search path

data(huagrahuma)
attach(huagrahuma)

## Initial exploration of the data:

str(huagrahuma)
topidx
parameters
rain

plot(rain, type="h")

## run the model and visualise the outcome:

Qsim <- topmodel(parameters, topidx, delay, rain, ET0)
plot(Qsim, type="l", col="red")
points(Qobs)

## Evaluate the model with a performance metric

NSeff(Qobs, Qsim)

############ PART 2: Sensitivity analysis ##############

## let's try first to vary only one parameter
## The function runif() samples randomly from a uniform distribution

parameters["m"] <- runif(1, min = 0, max = 0.1)
parameters["m"]

## Run the model and evaluate with the Nash – Sutcliffe efficiency metric:

Qsim <- topmodel(parameters, topidx, delay, rain, ET0)
NSeff(Qobs, Qsim)

## What value do you get? Do you think this is a good simulation?
## Verify by plotting:

plot(Qsim, type="l", col="red")
points(Qobs)

## Now sample all parameters at random. We take a sample size of 100

n <- 100

qs0 <- runif(n, min = 0.0001, max = 0.00025)
lnTe <- runif(n, min = -2, max = 3)
m <- runif(n, min = 0, max = 0.1)
Sr0 <- runif(n, min = 0, max = 0.2)
Srmax <- runif(n, min = 0, max = 0.1)
td <- runif(n, min = 0, max = 3)
vch <- runif(n, min = 100, max = 2500)
vr <- runif(n, min = 100, max = 2500)
k0 <- runif(n, min = 0, max = 10)
CD <- runif(n, min = 0, max = 5)
dt <- 0.25

parameters <- cbind(qs0,lnTe,m,Sr0,Srmax,td,vch,vr,k0,CD,dt)

## run the model and evaluate with the Nash – Sutcliffe efficiency metric:
## Note: the function accepts a table of parameter sets
## (one parameter set per row)

NS <- topmodel(parameters, topidx, delay, rain, ET0, Qobs = Qobs)
max(NS)

## visualisation of the sensitivity using dotty plots:

plot(lnTe, NS, ylim = c(0,1))

############ PART 3: GLUE uncertainty analysis ##############

## choose a behavioural threshold and remove the “bad” parameter sets:

parameters <- parameters[NS > 0.3,]
NS <- NS[NS > 0.3]

## generate predictions for the behavioural parameter sets:

Qsim <- topmodel(parameters,topidx,delay,rain,ET0)

## (have a look at the predictions for the first time step:)

hist(Qsim[1,])

## construct weights based on the performance measure

weights <- NS - 0.3
weights <- weights / sum(weights)

## make prediction boundaries by weighted quantile calculation
## (we need the Hmisc package for that)

limits <- apply(Qsim, 1, "wtd.quantile", weights = weights,
probs = c(0.05,0.95), normwt=T)

plot(limits[2,], type="l")
points(limits[1,], type="l")
points(Qobs, col="red")

## how many measurements fall outside of the prediction limits?

outside <- (Qobs > limits[2,]) | (Qobs < limits[1,])
summary(outside)

## width of the prediction boundaries

mean(limits[2,] - limits[1,]) / mean(Qobs, na.rm=T)

