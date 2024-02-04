## preparations #############################################################

load("Chapter15.RData")
str(maxima14days)           ## temperatures measured in 0.1 degree Celsius

## estimation of parameters  ################################################

library(SpatialExtremes)

## prepare data for input format for fitmaxstab
station.no <- nrow(stations)
max.data <- matrix(maxima14days$temp,nc=station.no) ## maxima14days already sorted by stn
max.data <- max.data[complete.cases(max.data),]

## euclidean coordinates of locations
locations <- as.matrix(stations[,c("lon","lat")])
coord <- locations
coord[,"lat"] <- lat.factor * locations[,"lat"]
colnames(coord)[2] <- "modif.lat"

## check if euclidean approximation is OK
library(geosphere)
geosphere.distances <- distm(stations[,c("lon","lat")],stations[,c("lon","lat")], fun=distVincentyEllipsoid)
euclidean.distances <- outer(1:18,1:18, function(i,j){sqrt((coord[i,"lon"]-coord[j,"lon"])^2 + (coord[i,"modif.lat"]-coord[j,"modif.lat"])^2)})
plot(geosphere.distances[all.pairs]/1e3,euclidean.distances[all.pairs]*dist.factor*1e2)
abline(0,1)

## weigh contribution from long distances down
weights <- exp(-geosphere.distances[all.pairs]/1e5) 
plot(geosphere.distances[all.pairs]/1e3,weights,ylim=c(0,1))

## fitting the max-stable process via pairwise composite likelihood
if (FALSE){
  logLikelihood <- rep(NA,6)
  for (i in 1:6){
    fit <- fitmaxstab(data=max.data, 
                      coord = coord, ## euclidean coord, see above/below
                      iso = TRUE,
                      cov.mod = "brown",
                      loc.form = loc ~ lon + lat + alt,
                      scale.form = scale~1, 
                      shape.form = shape ~ 1,
                      fit.marge = TRUE,
                      marg.cov = stations[c("lon","lat","alt")],
                      start=values.start,
                      weights=weights)
    logLikelihood[i] <- fit$logLik
    values.start <- as.list(fit$fitted.values)
  }
  plot(logLikelihood)
  ## convergence after 5 iterations  
}


fit.values <- fit$fitted.values
## marginal parameters
mu0 <- fit.values["locCoeff1"]
mu1 <- fit.values["locCoeff2"]
mu2 <- fit.values["locCoeff3"]
mu3 <- fit.values["locCoeff4"]
scale <- fit.values["scaleCoeff1"]
shape <- fit.values["shapeCoeff1"]
## variogram function
fitted.vario <- function(h){
  (sqrt(sum(h^2))/fit.values["range"])^fit.values["smooth"]
}

## evaluate estimated GEV location function on the inland grid
location.param <- with(inland.grid, mu0 + mu1*lon + mu2*lat + mu3*elevation_geonames)

## simulation from the fitted max-stable process ###################################

## make simulation algorithms available
source("Algorithms2.R") 

## inland coordinates
inland.grid.locations <- data.frame(lon = inland.grid$lon, modif.lat = lat.factor * inland.grid$lat)
coord.simu <- as.matrix(inland.grid.locations)

## the following example code for (just) 3 simulations may take 
## a few minutes / an hour (indicative time with i7 CORE, 8th Gen)
if (FALSE){
  simu <- simu_sumnormal(no.simu=3, coord=coord.simu, vario=fitted.vario, type="brownresnick",
                         rel.thresh=0.1, calc.err=F)
} 

## instead of simply setting no.simu=30000
## we ran 10 times 3000 simulations on a remote server and merged the output files
## exemplarily we take a look at 30 simulations
load("inlandBRsimu_30.RData")
simu.res <- simu$res

## transform the simulations to the estimated marginal distributions
no.simu <- dim(simu.res)[1]
simu.transformed <- simu.res
for (i in 1:no.simu){
  simu.transformed[i,] <- qgev(extRemes::pevd(simu.res[i,], loc=1, scale=1, shape=1,type="GEV"),
                               loc=location.param, scale=scale, shape=shape)
}

## question of interest ####################################

## joint distribution of the simulated maximum temperatures (over a 14 day period) 
## across the three subregions of the inland grid
## REMINDER: temperatures have been measured in 0.1 degree Celsius
if (FALSE){
  S1.max <- apply(simu.transformed[,inland.grid$region=="S1"], 1, max)/10
  S2.max <- apply(simu.transformed[,inland.grid$region=="S2"], 1, max)/10
  S3.max <- apply(simu.transformed[,inland.grid$region=="S3"], 1, max)/10
  area.maxima <- data.frame(S1=S1.max,S2=S2.max,S3=S3.max)
}
## already pre-computed for the 30000 simulations
str(area.maxima)

## threshold
thresh <- 40
joint.exceedance <- with(area.maxima, S1>thresh & S2 > thresh & S3> thresh)
p <- mean(joint.exceedance)
p

## estimated probability of a summer, during which
## there is an exceedance in all three regions S1, S2, S3
## simultaneously during at least one of the six 14 day intervals
1-(1-p)^6 
6*p

## REMINDER: CAUTION!
## This analysis is only indicative.
## See also Warning in Chapter 15.

