## Auxiliary plot functions for Chapter 15

plot.aggregation.2groups <- function(data, VAR, FUN1, VAR1, FUN2, VAR2, aggrFUN=mean, ... ){
  require(fields)
  stopifnot(c(VAR,VAR1,VAR2) %in% colnames(data))
  new.data <- aggregate(data[,VAR], 
                        by = list(Fct1 = FUN1(data[,VAR1]), 
                                  Fct2 = FUN2(data[,VAR2])), 
                        FUN = aggrFUN)
  new.data <- new.data[with(new.data, order(Fct2,Fct1)),]
  one.values <- unique(as.vector(new.data[,"Fct1"]))
  two.values <- unique(as.vector(new.data[,"Fct2"]))
  new.matrix <- matrix(new.data[,"x"],nr=length(two.values),nc=length(one.values),byrow=T)
  two.labels <- as.character(two.values)
  one.labels <- as.character(one.values)
  ltwo <- length(two.labels)
  lone <- length(one.labels)
  image.plot(1:lone,1:ltwo,t(new.matrix),xaxt="n",yaxt="n",xlab="",ylab="",...)
  axis(LEFT<-2, at=1:ltwo, labels=two.labels, las=HORIZONTAL <- 1,cex.axis=0.7)
  axis(BELOW<-1, at=1:lone, labels=one.labels, las=VERTICAL <-2, cex.axis=0.7)
}

plot.all.stations <- function(display="",title.text="Stations",
                              display.stations=T,display.lakes=T){
  stopifnot(exists("NL_sp") & exists("lakes_sp") & exists("coast_sp"))
  stopifnot(exists("inland.grid") & exists("stations"))
  plot(NL_sp, col="whitesmoke", 
       xlim=range(inland.grid$lon)+c(-0.1,0.1), ylim=range(inland.grid$lat)+c(-0.1,0.3)) 
  if(display.lakes){plot(lakes_sp, col="lightblue", add=TRUE)}
  plot(coast_sp, add = TRUE)
  box()
  title(title.text)
  ## points(inland.grid[,c("lon","lat")],pch=".")
  if (display.stations){
    points(stations[,c("lon","lat")],pch=20)
    if (display=="names") {text(stations[,c("lon","lat")],stations[,"name"],pos=4,cex=0.7,offset = 0.1)}
    if (display=="numbers") {text(stations[,c("lon","lat")],rownames(stations),pos=4,cex=0.7,offset = 0.1)}
  }
}

plot.one.station <- function(station){
  stopifnot(exists("stations"))
  one.station <- subset(stations,stn==station)
  plot.all.stations(title.text = one.station[,"name"])
  points(one.station[,c("lon","lat")],pch=21,cex=1.5,col="red",lwd=2)
}

plot.summer.temperature <- function(station){
  require(lubridate)
  stopifnot(exists("stations") & exists("summer.temperature"))
  par(mfrow=c(1,2))
  plot.one.station(station)
  plot.aggregation.2groups(data=subset(summer.temperature, stn==station), 
                           VAR="temp", 
                           FUN1=year, VAR1="date", 
                           FUN2=function(d){yday(d)-leap_year(d)-155}, VAR2="date")
  par(mfrow=c(1,1))
}

plot.max14days.temperature <- function(station){
  require(lubridate)
  stopifnot(exists("stations") & exists("maxima14days"))
  par(mfrow=c(1,2))
  plot.one.station(station)
  plot.aggregation.2groups(data=subset(maxima14days, stn==station), 
                           VAR="temp", 
                           FUN1=year, VAR1="date", 
                           FUN2=function(d){(yday(d)-leap_year(d)-156)/14+1}, VAR2="date")
  par(mfrow=c(1,1))
}

QQ.plot <- function(station,horiz=T){
  require(extRemes)
  require(car)
  if (horiz) {par(mfrow=c(1,2))} else {par(mfrow=c(2,1))}
  ## plot station
  if (horiz) {par(mai=c(0.5,0.4,0.3,0.7))} else {par(mai=c(0.2,0.8,0.3,0.2))} 
  plot.one.station(station)
  ## plot QQ plot diagnostic
  if (horiz) {par(mai=c(0.8,0.7,0.3,0.5),cex=0.8)} else {{par(mai=c(0.8,0.8,0.2,0.2),cex=1)}}
  stopifnot(exists("stations") & exists("maxima14days"))
  stopifnot(exists("location.scale.transformed.temp") & exists("shape"))
  x <- location.scale.transformed.temp[maxima14days$stn==station]
  qqPlot(x, distribution = "evd", loc=0, scale = 1, shape=shape, pch=20,
         xlab="Model quantiles",ylab="Empirical quantiles")
  par(mfrow=c(1,1))
}


## non-parametric estimators of bivariate extremal coefficients

## cooley et al biometrika 2006
F.madogram.estimator <- function(data){
  data <- as.matrix(data)
  stopifnot(dim(data)[2]==2)
  l <- nrow(data)
  nu <- mean(abs((rank(data[,1])/l)-(rank(data[,2])/l)))/2
  return((1+2*nu)/(1-2*nu)) 
}

## naveau et al biometrika 2009
lambda.madogram.estimator <- function(data){
  data <- as.matrix(data)
  stopifnot(dim(data)[2]==2)
  l <- nrow(data)
  nu <- mean(abs(sqrt(rank(data[,1])/l)-sqrt(rank(data[,2])/l)))/2
  return((2+3*nu)/(2-6*nu)) 
}

## caperaa et al 1997
cfg.estimator <- function(data){
  require(evd)  
  return(2*abvnonpar(data = data, method="cfg", empar=T))
}

