## preparations ######################################################################

load("Chapter15.RData")
library(rnaturalearth)  ## for spatial information in plots
source("AuxiliaryPlots.R") ## auxiliary plot functions

## get an overview over the data  ####################################

plot.all.stations(display="numbers")
## summer temperature data
plot.summer.temperature(stations$stn[1])
plot.summer.temperature(stations$stn[2])
## block maxima
plot.max14days.temperature(stations$stn[1])
plot.max14days.temperature(stations$stn[2])

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

## marginal checks
boxplot(matrix(with(maxima14days, 
                    temp - (mu0 + mu1 * lon + mu2 * lat + mu3 * alt)),
               nc=station.no))
abline(h=qgev(p=c(0.25,0.5,0.75), loc = 0, scale=scale, shape=shape),
       col="blue",lwd=2)

## QQ plots
location.scale.transformed.temp <- with(maxima14days, 
                                        (temp - (mu0 + mu1*lon + mu2*lat + mu3*alt))/scale) 
location.scale.transformed.empirical.quantiles <- sort(location.scale.transformed.temp)

QQ.plot(stations$stn[1])
QQ.plot(stations$stn[2])

## evaluate estimated GEV location function on the inland grid
location.param <- with(inland.grid, mu0 + mu1*lon + mu2*lat + mu3*elevation_geonames)
inland.grid.plot <- inland.grid
inland.grid.plot <- cbind(inland.grid,location.param)
colnames(inland.grid.plot)[c(1,2)]<-c("longitude","latitude")

library(latticeExtra)
## REMINDER: temperatures have been measured in 0.1 degree Celsius
## plot location parameter across inland grid points
levelplot(location.param/10 ~ longitude * latitude, inland.grid.plot, 
          panel=panel.levelplot.raster,  col.regions = fields::tim.colors(64), 
          par.settings=list(
            layout.widths=list(left.padding=0, right.padding=0.5),
            layout.heights=list(top.padding=0, bottom.padding=0, strip=0)))

## dependence checks

## theoretical bivariate extremal coefficient function
theoretical.ecf <- function(h){
  2*pnorm(sqrt(2*fitted.vario(h))/2)
}

## non-parametric estimates of bivariate extremal coefficients
no.pairs <- nrow(all.pairs)
ecf <-  rep(NA,no.pairs)
for (k in 1:no.pairs){
  ecf[k] <- F.madogram.estimator(max.data[,all.pairs[k,]])
}

## comparison of madogram estimates of extremal coefficients with model ecf
plot(geosphere.distances[all.pairs]/1e3,ecf, 
     xlim=c(0,250),ylim=c(1,1.3),pch=20,cex=1.5,
     ylab="extremal coefficients estimates",
     xlab="distance [km]",col=rgb(0,0,0,0.2))
lines((1:400)/100*dist.factor*1e2,
      Vectorize(theoretical.ecf)((1:400)/100),lwd=2)

## simulation from the fitted max-stable process ###################################

## load a fast random number generator
library(dqrng)
rnorm <- dqrnorm

## make simulation algorithms available
source("Algorithms2.R") 

#########################################
## first: another sanity check: simulations at the given points
stations.locations <- data.frame(lon = stations$lon, modif.lat = lat.factor * stations$lat)
coord.simu <- as.matrix(stations.locations)
control.simu  <- simu_extrfcts(type="brownresnick", no.simu=1e3, coord=coord.simu, vario=fitted.vario)
mean(1/control.simu$res)
## re-estimate extremal coefficients
ecf.simu <-  rep(NA,no.pairs)
for (k in 1:no.pairs){
  ecf.simu[k] <- F.madogram.estimator(control.simu$res[,all.pairs[k,]])
}
plot(geosphere.distances[all.pairs]/1e3,ecf.simu, 
     xlim=c(0,250),ylim=c(1,1.3),pch=21,
     ylab="extremal coefficients estimates",
     xlab="distance [km]",col=rgb(0,0,0,0.9))
lines((1:400)/100*dist.factor*1e2,Vectorize(theoretical.ecf)((1:400)/100),lwd=2)
## transform the simulations at these locations to the estimated marginal distributions
no.simu.at.stations <- dim(control.simu$res)[1]
control.simu.transformed <- control.simu$res
stations.location.parameters <- with(stations, mu0 + mu1*lon + mu2*lat + mu3*alt)
for (i in 1:no.simu.at.stations){
  control.simu.transformed[i,] <- qgev(pevd(control.simu$res[i,], loc=1, scale=1, shape=1,type="GEV"),
                               loc=stations.location.parameters, scale=scale, shape=shape)
}
## re-estimate parameters
control.fit <- fitmaxstab(data=control.simu.transformed, 
                  coord = coord, ## euclidean coord, see above/below
                  iso = TRUE,
                  cov.mod = "brown",
                  loc.form = loc ~ lon + lat + alt,
                  scale.form = scale~1, 
                  shape.form = shape ~ 1,
                  fit.marge = TRUE,
                  marg.cov = stations[c("lon","lat","alt")],
                  weights=weights)
control.fit.values <- control.fit$fitted.values
## for comparison
fit.values
control.fit.values
fit$std.err
##########################################

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

## simple marginal check on a sample of 30 locations
sample.grid.locations <- sample(1:dim(simu.res)[2],30,rep=F)
boxplot(log(simu.res[,sample.grid.locations])) 
abline(h=evd::qgumbel(c(0.25,0.5,0.75)),col="red")

## transform the simulations to the estimated marginal distributions
no.simu <- dim(simu.res)[1]
simu.transformed <- simu.res
for (i in 1:no.simu){
  simu.transformed[i,] <- qgev(extRemes::pevd(simu.res[i,], loc=1, scale=1, shape=1,type="GEV"),
                               loc=location.param, scale=scale, shape=shape)
}

## plot of three simulations
no.plot <- 3 
sample.simulations <- 1:3 ## sample(1:dim(simu.transformed)[1], no.plot, rep=F)
data <- data.frame(lon = rep(inland.grid$lon,no.plot),
                   lat = rep(inland.grid$lat,no.plot),
                   simu.tr = as.vector(t(simu.transformed[sample.simulations,]))/10,
                   no.simu = rep(1:no.plot,each=dim(inland.grid)[1]))
levelplot(simu.tr ~ lon * lat | factor(no.simu), data, 
          panel = panel.levelplot.raster, col.regions = fields::tim.colors(64))

## question of interest ####################################

## we consider three regions
par(mfrow=c(1,1))
plot.all.stations(display="numbers")
points(inland.grid[with(inland.grid,region=="S1"),], pch=".", col=rgb(1,0,0,0.5),cex=2)
points(inland.grid[with(inland.grid,region=="S2"),], pch=".", col=rgb(1,0,0,0.5),cex=2)
points(inland.grid[with(inland.grid,region=="S3"),], pch=".", col=rgb(1,0,0,0.5),cex=2)

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

## complementary information from the original data
## note: one station can typically not represent an areal maximum
## this serves as a sanity check only
S1.max.data <- subset(maxima14days,stn==stations[8,"stn"])$temp/10
S2.max.data <- subset(maxima14days,stn==stations[9,"stn"])$temp/10
S3.max.data <- subset(maxima14days,stn==stations[17,"stn"])$temp/10
station.representative.maxima.data <- data.frame(S1=S1.max.data,S2=S2.max.data,S3=S3.max.data)

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

## further illustration of the result
xlim <- ylim <- range(area.maxima)

par(mfrow=c(1,3))
plot(area.maxima$S1,area.maxima$S2,pch=".",cex=2,
     col=rgb(0,0,0,0.1),xlab="S1",ylab="S2",xlim=xlim,ylim=ylim)
abline(v=thresh,lty="dashed",col="red")
abline(h=thresh,lty="dashed",col="red")
abline(0,1,lty="dotted")
points(area.maxima$S1[joint.exceedance],area.maxima$S2[joint.exceedance],
     pch="+",cex=0.7,col="red")
points(station.representative.maxima.data$S1,
       station.representative.maxima.data$S2,pch="*",cex=1,col="magenta")

plot(area.maxima$S2,area.maxima$S3,pch=".",cex=2,
     col=rgb(0,0,0,0.1),xlab="S2",ylab="S3",xlim=xlim,ylim=ylim)
abline(v=thresh,lty="dashed",col="red")
abline(h=thresh,lty="dashed",col="red")
abline(0,1,lty="dotted")
points(area.maxima$S2[joint.exceedance],area.maxima$S3[joint.exceedance],
       pch="+",cex=0.7,col="red")
points(station.representative.maxima.data$S2,
       station.representative.maxima.data$S3,pch="*",cex=1,col="magenta")

plot(area.maxima$S1,area.maxima$S3,pch=".",cex=2,
     col=rgb(0,0,0,0.1),xlab="S1",ylab="S3",xlim=xlim,ylim=ylim)
abline(v=thresh,lty="dashed",col="red")
abline(h=thresh,lty="dashed",col="red")
abline(0,1,lty="dotted")
points(area.maxima$S1[joint.exceedance],area.maxima$S3[joint.exceedance],
       pch="+",cex=0.7,col="red")
points(station.representative.maxima.data$S1,
       station.representative.maxima.data$S3,pch="*",cex=1,col="magenta")


