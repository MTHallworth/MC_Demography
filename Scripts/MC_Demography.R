
library(raster)
library(rgdal)
library(sp)
library(RColorBrewer)
library(maptools)
library(rgeos)
library(geosphere)
library(SGAT)
library(TwGeos)
library(MASS)
library(mth)

##############################################################################
#
#
# Light-Level Geolocator Analysis 
#
#
#############################################################################

Americas<-shapefile("Spatial_Layers/Americas.shp")
states<-shapefile("Spatial_Layers/st99_d00.shp")

# Black-throated blue warbler 

BTBWdist<-shapefile("Spatial_Layers/BTBWdist.shp")
BTBWdist<-gUnaryUnion(BTBWdist, id = NULL)

# Read in the data from LUX files 
BTBW_HBEF<-list.files("Data/GLdata/BTBW/LigFiles", pattern = ".lig", full.names = TRUE)

BTBWnames<-list.files("Data/GLdata/BTBW/LigFiles",pattern = ".lig")

# Read just the file names for Bird ID
BirdId<- list.files("Data/GLdata/BTBW/LigFiles", pattern = ".lig")

# Determine the number of birds
nBirds<-length(BirdId)

# Read in the lux file

BTBWdata <- lapply(BTBW_HBEF, readLig, skip=1)  

head(BTBWdata[[1]])

# Set the capture coordinates for each bird #
CapLocs<-array(NA,c(nBirds,2))

# give row names the same as the BTBW names files to keep everything organized.
rownames(CapLocs)<-BTBWnames

# Capture locations
CapLocs[1,]<- cbind(-71.77347,43.93279)
CapLocs[2,]<- cbind(-71.77193,43.93485)
CapLocs[3,]<- cbind(-71.77704,43.93618)
CapLocs[4,]<- cbind(-71.78147,43.93536)


## ----echo=FALSE----------------------------------------------------------
twl <- twlEdit <- vector('list',nBirds)
seed <- as.POSIXct("2015-11-01 04:00", origin  = "1970-01-01", tz = "GMT")

for(i in 1:nBirds){

twl[[i]]  <- findTwilights(tagdata = BTBWdata[[i]], 
                      threshold = 1, 
                      include = seed,
                      dark.min = 240) # minimum dark period in minutes

twlEdit[[i]] <- twilightEdit(twilights = twl[[i]], 
                    window = 4,           # two days before and two days after
                    outlier.mins = 30,    # difference in mins
                    stationary.mins = 15, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)

twlEdit [[i]] <- twilightAdjust(twilights = twlEdit[[i]], 
                      interval = 120) # The unit here is seconds
}

twlEdit[[3]] <- twlEdit[[3]][8:nrow(twlEdit[[3]]),]
twlEdit[[4]] <- twlEdit[[4]][5:nrow(twlEdit[[4]]),]

# Create a vector with the dates known to be at deployment #
calib.dates <- vector('list',nBirds)

for(i in 1:nBirds){
  calib.dates[[i]] <- c(strptime(twlEdit[[i]][1,1],format="%Y-%m-%d"),as.POSIXct("2015-08-15"))
}

calibration.data<-vector('list',nBirds)

for(i in 1:nBirds){
  calibration.data[[i]]<-subset(twlEdit[[i]],twlEdit[[i]]$Twilight>=calib.dates[[i]][1] & twlEdit[[i]]$Twilight<=calib.dates[[i]][2])
}



# Generate empty lists to store data 
sun<-z<-twl_t<-twl_deviation<-fitml<-alpha<-vector('list',nBirds)

# Determine the sun elevation angle - here called the Zenith angle #

for(i in 1:nBirds){
  
  # Calculate solar time from calibration data 
  sun[[i]]  <- solar(calibration.data[[i]][,1])
  
  # Adjust the solar zenith angle for atmospheric refraction
  z[[i]] <- refracted( zenith(sun = sun[[i]],
                              lon = CapLocs[i,1], 
                              lat = CapLocs[i,2]))
  
  twl_t[[i]] <- twilight(tm = calibration.data[[i]][,1],
                         lon = CapLocs[i,1], 
                         lat = CapLocs[i,2], 
                         rise = calibration.data[[i]][,2],
                         zenith = quantile(z[[i]],probs=0.5))
  
  # Determine the difference in minutes from when the sun rose and the geolocator said it rose 
  twl_deviation[[i]] <- ifelse(calibration.data[[i]]$Rise, as.numeric(difftime(calibration.data[[i]][,1], twl_t[[i]], units = "mins")),
                               as.numeric(difftime(twl_t[[i]], calibration.data[[i]][,1], units = "mins")))
  
  # Throw out values less than 0 - These values are not valid 
  twl_deviation[[i]]<-subset(twl_deviation[[i]], subset=twl_deviation[[i]]>=0)
  
  # Describe the distribution of the error 
  fitml[[i]] <- fitdistr(twl_deviation[[i]], "log-Normal")
  
  # save the Twilight model parameters
  alpha[[i]]<- c(fitml[[i]]$estimate[1], fitml[[i]]$estimate[2]) 
}


b<-unlist(twl_deviation)
b[b>200]<-NA

cols<-c("red","blue","green","yellow","orange","purple","brown","gray","black","pink")
seq <- seq(0,60, length = 100)
par(mfrow=c(1,2),mar=c(4,4,0,0))
hist(b, freq = F,
     yaxt="n",
     ylim = c(0, 0.15),
     xlim = c(0, 60),
     breaks=15,
     col="gray",
     main = "",
     xlab = "Twilight error (mins)")
axis(2,las=2)
for(i in 1:nBirds){
  lines(seq, dlnorm(seq, alpha[[i]][1], alpha[[i]][2]), col = cols[i], lwd = 3, lty = 2)
}

#Zenith angle plot
par(bty="l")
plot(median(z[[1]],na.rm=TRUE),xlim=c(1,nBirds),ylim=c(80,100),pch=19,ylab="Zenith Angle",xlab="BTBW",col=cols[1])
segments(1,quantile(z[[1]],probs=0.025),1,quantile(z[[1]],probs=0.975),col=cols[1])
for(i in 2:nBirds){
  par(new = TRUE)
  plot(median(z[[i]],na.rm=TRUE)~i,xlim=c(1,nBirds),ylim=c(80,100),pch=19,yaxt="n",xaxt="n",ylab="",xlab="",col=cols[i])
  segments(i,quantile(z[[i]],probs=0.025),i,quantile(z[[i]],probs=0.975),col=cols[i])
}

# Create empty vectors to store objects #
d.twl<-path<-vector('list',nBirds)

zenith0<-zenith1<-rep(NA,nBirds)

# loop through the birds #
for(i in 1:nBirds){
  # Store the zenith (sun-elevation angle)
  zenith0[i] <-quantile(z[[i]],prob=0.5)
  zenith1[i]<-quantile(z[[i]],prob=0.95)
}

twlEdit[[3]] <- twlEdit[[3]][8:nrow(twlEdit[[3]]),]

# subset the twilight file for dates after the first calibration date (presumably the deployment date)  
# and exclude points that were deleted  
# note we didn't delete any transitions here

tol <- array(0.08, c(nBirds,2))
# Manual Adjustments
tol[1,1] <- 0.12
tol[2,1] <- 0.16
tol[3,1] <- 0.12
tol[4,1] <- 0.12

tol[1,2] <- 0.23
tol[2,2] <- 0.2
tol[3,2] <- 0.142
tol[4,2] <- 0.185


for(i in 1:nBirds){  
  twlEdit[[i]]<-subset(twlEdit[[i]],twlEdit[[i]]$Twilight>=calib.dates[[i]][1] & !Deleted)
  
  path[[i]] <- thresholdPath(twilight = twlEdit[[i]]$Twilight,
                             rise = twlEdit[[i]]$Rise,
                             zenith = zenith0[i],
                             tol = tol[i,])
}

for(i in 1:nBirds){
  print(BirdId[[i]])
  layout(matrix(c(1,3,
                  2,3), 2, 2, byrow = TRUE))
  par(mar=c(2,4,2,0))
  plot(path[[i]]$time, path[[i]]$x[, 2], type = "b", pch = 16, cex = 0.5, ylab = "Latitude", xlab = '',xaxt="n",
       col=ifelse(path[[i]]$time<as.POSIXct("2016-01-01",format="%Y-%m-%d"),"blue","green"))
  abline(h = CapLocs[i,2])
  abline(v = as.POSIXct("2015-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2016-03-20"),col="red",lty=2,lwd=1.5)
  par(mar=c(2,4,2,0))
  plot(path[[i]]$time, path[[i]]$x[, 1], type = "b", pch = 16, cex = 0.5, ylab = "Longitude", xlab = '',
       col=ifelse(path[[i]]$time<as.POSIXct("2016-01-01",format="%Y-%m-%d"),"blue","green"))
  abline(h = CapLocs[i,1])
  abline(v = as.POSIXct("2015-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2016-03-20"),col="red",lty=2,lwd=1.5)
  plot(Americas, col = "grey95",xlim = c(-120,-60),ylim=c(0,40))
  plot(BTBWdist, col = "grey45",border="grey45",add=TRUE)
  box()
  lines(path[[i]]$x, col = ifelse(path[[i]]$time<as.POSIXct("2016-01-01",format="%Y-%m-%d"),"blue","green"))
  points(path[[i]]$x, pch = 16, cex = 0.5, col=ifelse(path[[i]]$time<as.POSIXct("2016-01-01",format="%Y-%m-%d"),"blue","green"))
}


x0 <- z0 <- vector('list',nBirds)

for(i in 1:nBirds){
  # Take the location estimates created above
  x0[[i]]<- path[[i]]$x
  
  # the model also needs the mid-points - generate those here
  z0[[i]]<- trackMidpts(x0[[i]])
}

beta <- c(0.7, 0.08)

fixedx <- vector('list',nBirds)

for(i in 1:nBirds){
  fixedx[[i]]<- rep(FALSE, nrow(x0[[i]]))
  # from capture to Sept 1. 
  fixedx[[i]][c(1:max(which(twlEdit[[i]][,1] < as.POSIXct("2015-09-01",format = "%Y-%m-%d"))),
               (nrow(x0[[i]])-3):nrow(x0[[i]]))] <- TRUE
  
  x0[[i]][fixedx[[i]], 1] <- CapLocs[i,1]
  x0[[i]][fixedx[[i]], 2] <- CapLocs[i,2]
  
  z0[[i]] <- trackMidpts(x0[[i]]) # update z0 positions
}


## Function to construct a land/sea mask
distribution.mask <- function(xlim, ylim, n = 4, land = TRUE, shape) {
  r <- raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1], 
              xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(shape))
  r <- cover(rasterize(shape, shift = c(-360, 0), r, 1, silent = TRUE), 
             rasterize(shape, r, 1, silent = TRUE), rasterize(shape, r, 1, silent = TRUE))
  r <- as.matrix(is.na(r))[nrow(r):1, ]
  if (land) 
    r <- !r
  xbin <- seq(xlim[1], xlim[2], length = ncol(r) + 1)
  ybin <- seq(ylim[1], ylim[2], length = nrow(r) + 1)
  
  function(p) {
    r[cbind(.bincode(p[, 2], ybin), .bincode(p[, 1], xbin))]
  }
}

WGS84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(Americas)<-WGS84
crs(BTBWdist)<-WGS84

xlim=c(-100,-60)
ylim=c(0,55)

## Define mask for Ovenbird distribution
is.dist <- distribution.mask(shape=BTBWdist,
                             xlim = xlim,
                             ylim = ylim,
                             n = 4,
                             land = TRUE)

# Define the log prior for x and z
log.prior <- function(p) {
  f <- is.dist(p)
  ifelse(f | is.na(f), 0, -10)
}


# Define the threshold model - slimilar to above #
model <-  vector('list', nBirds)

for(i in 1:nBirds){
  model[[i]]<- thresholdModel(twilight = twlEdit[[i]]$Twilight,
                              rise = twlEdit[[i]]$Rise,
                              twilight.model = "ModifiedLogNormal",
                              alpha = alpha[[i]],
                              beta = beta,
                              # Here is where we set the constraints for land
                              logp.x = log.prior, logp.z = log.prior, 
                              x0 = x0[[i]],
                              z0 = z0[[i]],
                              zenith = zenith1[i],
                              fixedx = fixedx[[i]])
}

# This defines the error distribution around each location #
proposal.x <- proposal.z <- vector('list',nBirds)

for(i in 1:nBirds){
  proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0[[i]]))
  proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0[[i]]))
}

fit <- xsum <- zsum <- vector('list', nBirds)

for(i in 1:nBirds){
  
  fit[[i]] <- estelleMetropolis(model = model[[i]],
                                proposal.x = proposal.x[[i]],
                                proposal.z = proposal.z[[i]],
                                iters = 2000, # This value sets the number of iterations to run
                                thin = 10,
                                chains = 3)

xsum[[i]] <- locationSummary(fit[[i]]$x)
zsum[[i]] <- locationSummary(fit[[i]]$z)  

### Fine Tuning 

proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0[[i]]))
proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0[[i]]))

fit[[i]] <- estelleMetropolis(model = model[[i]],
                              proposal.x = proposal.x[[i]],
                              proposal.z = proposal.z[[i]],
                              x0 = cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%'),
                              z0 = cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%'),
                              iters=2000, # This value sets the number of iterations to run
                              thin=10,
                              chains=3)

# Final Run

proposal.x[[i]] <- mvnorm(chainCov(fit[[i]]$x),s=0.1)
proposal.z[[i]] <- mvnorm(chainCov(fit[[i]]$z),s=0.1)

xsum[[i]] <- locationSummary(fit[[i]]$x)
zsum[[i]] <- locationSummary(fit[[i]]$z)  

# Note the increase in number of interations - this takes a bit longer to run
fit[[i]] <- estelleMetropolis(model = model[[i]],
                              proposal.x = proposal.x[[i]],
                              proposal.z = proposal.z[[i]],
                              x0 = cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%'),
                              z0 = cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%'),
                              iters=5000,  # This value sets the number of iterations to run
                              thin=10,
                              chains=3)

}
## Inital results 

xlim <- c(-100,-60)
ylim <- c(0,55)
 

# This step makes an empty raster #
r <- raster(res = c(0.25,0.25),
            xmn = xlim[1],
            xmx = xlim[2], 
            ymn = ylim[1],
            ymx = ylim[2])

S <- nonbreed <- breed <- vector('list',nBirds)

for(i in 1:nBirds){
  S[[i]] <- slices(type="intermediate",
                   breaks="day",
                   mcmc=fit[[i]],
                   grid=r,
                   weight = rep(0.5,length(fit[[i]][[1]]$time)))
}

###

DATES <-  tm_breed <-tm_winter<-breed<-winter<- vector('list',nBirds)
Aug31 <- Nov01  <- April1<- rep(NA,nBirds)

# Create a vector of the Dates in the file - here we set rise=true because the data has both rise and sunset

for(i in 1:nBirds){ 
  DATES[[i]] <- S[[i]]$mcmc[[1]]$time[ which( S[[i]]$mcmc[[1]]$rise==TRUE) ]
  
  
  # Here we specify our dates of interest. 
  ReleaseDay<-1
  Aug31[i]<-which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct("2015-08-31"))
  
  
  # Breeding 
  # Get time when on the breeding grounds 
  tm_breed[[i]]<-sliceInterval(S[[i]],k=c(ReleaseDay:Aug31[i]))
  
  
  # "Slice" the data and save all dates between Release date and July 31 2011, and June 1 2012 until capture.
  breed[[i]]<-slice(S[[i]],k=c(ReleaseDay:Aug31[i]))
  
  print(BirdId[i])
  
  plot(BTBWdist,col="gray74")
  
  plot(breed[[i]],useRaster=TRUE,
       axes=FALSE, add=TRUE,
       legend=FALSE,
       col=rev(bpy.colors(50)),
       cex.axis=0.7)
  
  plot(Americas,border="gray",add=TRUE)
  
  # Non-breeding 
  
  Nov01[i] <- which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct("2015-12-01"))
  April1[i] <- which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct("2016-04-01"))
  
  tm_winter[[i]]<- sliceInterval(S[[i]],k=c(Nov01[i]:April1[i]))
  
  # "Slice" data and merge between Nov 01 2014 and Feb 31 2015.
  winter[[i]]<- slice(S[[i]],k=c(Nov01[[i]],April1[[i]]))
  
  #plot(BTBWdist,border="black",ylim=ylim,xlim=xlim)
  plot(winter[[i]],useRaster=TRUE,add=TRUE,
       axes=FALSE,
       legend=FALSE,
       col=rev(bpy.colors(50)),
       cex.axis=0.7)
  
  plot(SpatialPoints(cbind(CapLocs[i,1],CapLocs[i,2])),add=TRUE,pch=19,cex=1.5)
  box()
} 


scaledWinter<-scaledBreeding<-vector('list',nBirds)
for(i in 1:nBirds){
scaledBreeding[[i]]<-breed[[i]]/cellStats(breed[[i]],max)
scaledWinter[[i]]<-winter[[i]]/cellStats(winter[[i]],max)
}

scaledBreeding<-stack(scaledBreeding)
scaledWinter<-stack(scaledWinter)

sumBirdsBreed_btbw<-sum(scaledBreeding,na.rm=TRUE)
sumBirdsWinter_btbw<-sum(scaledWinter,na.rm=TRUE)

sumBirdsBreed_btbw[sumBirdsBreed_btbw==0]<-NA

E<-sumBirdsBreed_btbw/nBirds

prOrigin<-sumBirdsWinter_btbw/4

MC_btbw<-(cellStats(prOrigin,max) - 1/nBirds) / (1 - 1/nBirds) 

########################################################################
#
#
#                               OVENBIRDS
#
#########################################################################

OVENdist<-shapefile("Spatial_Layers/OVENdist.shp")
OVENdist<-gUnaryUnion(OVENdist, id = NULL)

# Read in the data from LUX files 
OVEN_HBEF<-list.files("Data/GLdata/OVEN/LigFiles", pattern = ".lig", full.names = TRUE)

OVENnames<-list.files("Data/GLdata/OVEN/LigFiles",pattern = ".lig")

# Read just the file names for Bird ID
BirdId<- list.files("Data/GLdata/OVEN/LigFiles", pattern = ".lig")

# Determine the number of birds
nBirds<-length(BirdId)

# Loop through all the files and read them in as LUX files #

OVENdata <- lapply(OVEN_HBEF,readLig,skip=1)  

head(OVENdata[[1]])

# Set the capture coordinates for each bird #
CapLocs<-array(NA,c(nBirds,2))

# give row names the same as the OVEN names files to keep everything organized.
rownames(CapLocs)<-OVENnames

# Capture locations
CapLocs<- cbind(-71.73,43.94)

lapply(OVENdata,head)
## ----echo=FALSE----------------------------------------------------------
twl <- twlEdit <- seed <- vector('list',nBirds)
for(i in 1:nBirds){
if(i %in% c(1,3,8)){
seed[[i]] <- as.POSIXct("2010-11-01 04:00", origin  = "1970-01-01", tz = "GMT")
} 
if(i %in% c(2,4:7,9:20)){
seed[[i]] <- as.POSIXct("2012-01-01 04:00", origin = "1970-01-01", tz = "GMT")
}
}


for(i in 1:nBirds){

twl[[i]]  <- findTwilights(tagdata = OVENdata[[i]], 
                      threshold = 1, 
                      include = seed[[i]],
                      dark.min = 240) # minimum dark period in minutes

twlEdit[[i]] <- twilightEdit(twilights = twl[[i]], 
                    window = 4,           # two days before and two days after
                    outlier.mins = 30,    # difference in mins
                    stationary.mins = 15, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)

twlEdit [[i]] <- twilightAdjust(twilights = twlEdit[[i]], 
                      interval = 120) # The unit here is seconds
}

# Save only the twilights that were not deleted #

twlEdit <- lapply(twlEdit, subset, !Deleted)

twlEdit[[1]] <- twlEdit[[1]][4:nrow(twlEdit[[1]]),]

twlEdit[[3]] <- twlEdit[[3]][7:nrow(twlEdit[[3]]),]

twlEdit[[5]] <- twlEdit[[5]][1:(nrow(twlEdit[[5]])-7),]
twlEdit[[6]] <- twlEdit[[6]][6:nrow(twlEdit[[6]]),]
twlEdit[[7]] <- twlEdit[[7]][7:nrow(twlEdit[[7]]),]
twlEdit[[8]] <- twlEdit[[8]][4:nrow(twlEdit[[8]]),]
twlEdit[[9]] <- twlEdit[[9]][1:(nrow(twlEdit[[9]])-4),]
twlEdit[[10]] <- twlEdit[[10]][4:(nrow(twlEdit[[10]])-3),]
twlEdit[[11]] <- twlEdit[[11]][6:(nrow(twlEdit[[11]])-3),]
twlEdit[[12]] <- twlEdit[[12]][4:(nrow(twlEdit[[12]])-14),]
twlEdit[[13]] <- twlEdit[[13]][5:nrow(twlEdit[[13]]),]
twlEdit[[14]] <- twlEdit[[14]][5:(nrow(twlEdit[[14]])-5),]
twlEdit[[15]] <- twlEdit[[15]][4:(nrow(twlEdit[[15]])-5),]
twlEdit[[16]] <- twlEdit[[16]][1:(nrow(twlEdit[[16]])-5),]

# Create a vector with the dates known to be at deployment #
calib.dates <- vector('list',nBirds)

# Set the calibration dates to determine error in twilight times from civil twilight #
# set year for the calibration # a few birds were captured in 2010. 
year<-rep("2011-07-31",20)
year[c(1,3,8)]<-"2010-07-31"

for(i in 1:nBirds){
  calib.dates[[i]] <- c(strptime(twl[[i]][1,1],format="%Y-%m-%d"),as.POSIXct(year[i]))
}

calibration.data<-vector('list',nBirds)

for(i in 1:nBirds){
  calibration.data[[i]]<-subset(twl[[i]],twl[[i]]$Twilight>=calib.dates[[i]][1] & twl[[i]]$Twilight<=calib.dates[[i]][2])
}



# Generate empty lists to store data 
sun<-z<-twl_t<-twl_deviation<-fitml<-alpha<-vector('list',nBirds)

# Determine the sun elevation angle - here called the Zenith angle #

for(i in 1:nBirds){
  
  # Calculate solar time from calibration data 
  sun[[i]]  <- solar(calibration.data[[i]][,1])
  
  # Adjust the solar zenith angle for atmospheric refraction
  z[[i]] <- refracted( zenith(sun = sun[[i]],
                              lon = CapLocs[1], 
                              lat = CapLocs[2]))
  
  twl_t[[i]] <- twilight(tm = calibration.data[[i]][,1],
                         lon = CapLocs[1], 
                         lat = CapLocs[2], 
                         rise = calibration.data[[i]][,2],
                         zenith = quantile(z[[i]],probs=0.5))
  
  # Determine the difference in minutes from when the sun rose and the geolocator said it rose 
  twl_deviation[[i]] <- ifelse(calibration.data[[i]]$Rise, as.numeric(difftime(calibration.data[[i]][,1], twl_t[[i]], units = "mins")),
                               as.numeric(difftime(twl_t[[i]], calibration.data[[i]][,1], units = "mins")))
  
  # Throw out values less than 0 - These values are not valid 
  twl_deviation[[i]]<-subset(twl_deviation[[i]], subset=twl_deviation[[i]]>=0)
  
  # Describe the distribution of the error 
  fitml[[i]] <- fitdistr(twl_deviation[[i]], "log-Normal")
  
  # save the Twilight model parameters
  alpha[[i]]<- c(fitml[[i]]$estimate[1], fitml[[i]]$estimate[2]) 
}

b<-unlist(twl_deviation)
b[b>100]<-NA

cols<-colorRampPalette(c("blue","purple","red"))
seq <- seq(0,60, length = 100)
par(mfrow=c(1,2),mar=c(4,4,0,0))
hist(b, freq = F,
     yaxt="n",
     ylim = c(0, 0.15),
     xlim = c(0, 60),
     breaks=25,
     col="gray",
     main = "",
     xlab = "Twilight error (mins)")
axis(2,las=2)
for(i in 1:nBirds){
  lines(seq, dlnorm(seq, alpha[[i]][1], alpha[[i]][2]), col = cols(nBirds)[as.numeric(cut(i,breaks = seq(0,nBirds,1)))], lwd = 3, lty = 2)
}

#Zenith angle plot
par(bty="l")
plot(median(z[[1]],na.rm=TRUE),
     xlim=c(1,nBirds),
     ylim=c(80,100),
     pch=19,
     ylab="Zenith Angle",
     xlab="OVEN",
     col=cols(nBirds)[as.numeric(cut(1,breaks = seq(0,nBirds,1)))])
segments(1,quantile(z[[1]],probs=0.025),
         1,quantile(z[[1]],probs=0.975),
         col=cols(nBirds)[as.numeric(cut(1,breaks = seq(0,nBirds,1)))])
for(i in 2:nBirds){
  par(new = TRUE)
  plot(median(z[[i]],na.rm=TRUE)~i,
       xlim=c(1,nBirds),
       ylim=c(80,100),
       pch=19,
       yaxt="n",
       xaxt="n",
       ylab="",
       xlab="",
       col=cols(nBirds)[as.numeric(cut(i,breaks = seq(0,nBirds,1)))])
  segments(i,quantile(z[[i]],probs=0.025),
           i,quantile(z[[i]],probs=0.975),
           col = cols(nBirds)[as.numeric(cut(i,breaks = seq(0,nBirds,1)))])
}

# Create empty vectors to store objects #
path <- Zeniths <- vector('list',nBirds)

zenith0<-zenith1<-rep(NA,nBirds)


# loop through the birds #
for(i in 1:nBirds){
  # Store the zenith (sun-elevation angle)
  zenith0[i] <-quantile(z[[i]],prob=0.5)
  zenith1[i]<-quantile(z[[i]],prob=0.95)
Zeniths[[i]] <- rep(zenith0[i],nrow(twlEdit[[i]]))
if(i %in% c(1,3,8)){
Zeniths[[i]][which(twlEdit[[i]][,1] > as.POSIXct("2010-09-15", tz= "GMT") &
                   twlEdit[[i]][,1] < as.POSIXct("2011-05-01", tz = "GMT"))] <- zenith1[i]
}
if(i %in% c(2,4:7,9:20)){
Zeniths[[i]][which(twlEdit[[i]][,1] > as.POSIXct("2011-09-15", tz= "GMT") &
                   twlEdit[[i]][,1] < as.POSIXct("2012-05-01", tz = "GMT"))] <- zenith1[i]
}
}


# subset the twilight file for dates after the first calibration date (presumably the deployment date)  
# and exclude points that were deleted  
# note we didn't delete any transitions here

tol <- array(0.08, c(nBirds,2))

# Manual Adjustments 
# These adjustments were made to maximize the locations on land

tol[1,1] <- 0.26
tol[2,1] <- 0.14
tol[3,1] <- 0.11
tol[4,1] <- 0.125
tol[5,1] <- 0.15
tol[6,1] <- 0.18
tol[7,1] <- 0.14
tol[8,1] <- 0.18
tol[9,1] <- 0.1
tol[10,1] <- 0.08
tol[11,1] <- 0.16
tol[12,1] <- 0.12
tol[13,1] <- 0.14
tol[14,1] <- 0.13
tol[15,1] <- 0.12
tol[16,1] <- 0.1
tol[17,1]
tol[18,1]
tol[19,1]
tol[20,1]


tol[1,2] <- 0.25
tol[2,2] <- 0.23
tol[3,2] <- 0.383  # GL died during non-breeding season
tol[4,2] <- 0.18
tol[5,2] <- 0.2
tol[6,2] <- 0.2
tol[7,2] <- 0.21
tol[8,2] <- 0.24
tol[9,2] <- 0.25
tol[10,2] <- 0.23
tol[11,2] <- 0.26
tol[12,2] <- 0.24
tol[13,2] <- 0.19
tol[14,2] <- 0.2205
tol[15,2] <- 0.1
tol[16,2] <- 0.1
tol[17,2]
tol[18,2]
tol[19,2]
tol[20,2]

for(i in 1:nBirds){  
  
  path[[i]] <- thresholdPath(twilight = twlEdit[[i]]$Twilight,
                             rise = twlEdit[[i]]$Rise,
                             zenith = zenith0[i],
                             tol = tol[i,])
}
cols <- colorRampPalette(c("blue","purple","red"))

#for(i in 1:nBirds){
  i = 16
  print(BirdId[[i]])
  layout(matrix(c(1,3,
                  2,3), 2, 2, byrow = TRUE))
  par(mar=c(2,4,2,0))
  plot(path[[i]]$time, path[[i]]$x[, 2], 
       type = "b", 
       pch = 16, 
       cex = 0.5, 
       ylab = "Latitude",
       xlab = '',
       xaxt="n",
       col=cols(40)[as.numeric(cut(abs(solar(twlEdit[[i]]$Twilight[twlEdit[[i]]$Rise])$sinSolarDec),breaks = seq(0,0.4,0.01)))])
  abline(h = CapLocs[1,2])
  abline(v = as.POSIXct("2015-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2016-03-20"),col="red",lty=2,lwd=1.5)
  par(mar=c(2,4,2,0))
  plot(path[[i]]$time, path[[i]]$x[, 1], 
       type = "b",
       pch = 16,
       cex = 0.5,
       ylab = "Longitude",
       xlab = '',
       col= cols(40)[as.numeric(cut(abs(solar(twlEdit[[i]]$Twilight[twlEdit[[i]]$Rise])$sinSolarDec),breaks = seq(0,0.4,0.01)))])
  abline(h = CapLocs[1,1])
  abline(v = as.POSIXct("2015-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2016-03-20"),col="red",lty=2,lwd=1.5)
  plot(Americas,
       col = "grey95",
       xlim = c(-120,-60),
       ylim=c(0,40))
  plot(OVENdist, 
       col = "grey45",
       border="grey45",
       add=TRUE)
  box()
  lines(path[[i]]$x,
        col = ifelse(path[[i]]$time<as.POSIXct("2016-01-01",format="%Y-%m-%d"),"blue","green"))
  points(path[[i]]$x, 
        pch = 16, cex = 0.5, col=ifelse(path[[i]]$time<as.POSIXct("2016-01-01",format="%Y-%m-%d"),"blue","green"))
#Sys.sleep(3)
#}

x0 <- z0 <- vector('list',nBirds)

for(i in 1:nBirds){
  # Take the location estimates created above
  x0[[i]]<- path[[i]]$x
  
  # the model also needs the mid-points - generate those here
  z0[[i]]<- trackMidpts(x0[[i]])
}

beta <- c(0.7, 0.08)

fixedx <- vector('list',nBirds)

for(i in 1:nBirds){
  fixedx[[i]]<- rep(FALSE, nrow(x0[[i]]))
  
  fixedx[[i]][c(1:10,(nrow(x0[[i]])-3):nrow(x0[[i]]))] <- TRUE
  
  x0[[i]][fixedx[[i]], 1] <- CapLocs[1]
  x0[[i]][fixedx[[i]], 2] <- CapLocs[2]
  
  z0[[i]] <- trackMidpts(x0[[i]]) # update z0 positions
}


## Function to construct a land/sea mask
distribution.mask <- function(xlim, ylim, n = 4, land = TRUE, shape) {
  r <- raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1], 
              xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(shape))
  r <- cover(rasterize(shape, shift = c(-360, 0), r, 1, silent = TRUE), 
             rasterize(shape, r, 1, silent = TRUE), rasterize(shape, r, 1, silent = TRUE))
  r <- as.matrix(is.na(r))[nrow(r):1, ]
  if (land) 
    r <- !r
  xbin <- seq(xlim[1], xlim[2], length = ncol(r) + 1)
  ybin <- seq(ylim[1], ylim[2], length = nrow(r) + 1)
  
  function(p) {
    r[cbind(.bincode(p[, 2], ybin), .bincode(p[, 1], xbin))]
  }
}

WGS84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(Americas)<-WGS84
crs(OVENdist)<-WGS84

xlim=c(-100,-60)
ylim=c(0,55)

## Define mask for Ovenbird distribution
is.dist <- distribution.mask(shape=OVENdist,
                             xlim = xlim,
                             ylim = ylim,
                             n = 4,
                             land = TRUE)

# Define the log prior for x and z
log.prior <- function(p) {
  f <- is.dist(p)
  ifelse(f | is.na(f), 0, -10)
}


# Define the threshold model - slimilar to above #
model <-  vector('list', nBirds)

# Set up changes in alpha and zenith angles here between breeding and non-breeding # 
# first - find the dates you want to change #
# I want the values to change from breeding to non-breeding on Nov1 and change back on April1 # 

November01<-rep("2011-11-01",20)
November01[c(1,3,8)]<-"2010-11-01"
Apr1<-rep("2012-04-01",20)
Apr1[c(1,3,8)]<-"2011-04-01"

which(strptime(path[[4]]$time,format="%Y-%m-%d")=="2012-04-01")

winterZ<-quantile(unlist(z_NB),probs=0.95)

alphaBoth<-zenithBoth<-vector('list',20)

for(i in 1:nBirds){

alphaBoth[[i]]<-cbind(rep(alpha[[i]][1],length(path[[i]]$time)),rep(alpha[[i]][2],length(path[[i]]$time)))
zenithBoth[[i]]<-rep(zenith1[[i]],length(path[[i]]$time))

alphaBoth[[i]][which(strptime(path[[i]]$time,format="%Y-%m-%d")==November01[i])[1]:
               which(strptime(path[[i]]$time,format="%Y-%m-%d")==Apr1[i])[1],1]<-meanNBalpha[1]
alphaBoth[[i]][which(strptime(path[[i]]$time,format="%Y-%m-%d")==November01[i])[1]:
               which(strptime(path[[i]]$time,format="%Y-%m-%d")==Apr1[i])[1],2]<-meanNBalpha[2]

zenithBoth[[i]][which(strptime(path[[i]]$time,format="%Y-%m-%d")==November01[i])[1]:
               which(strptime(path[[i]]$time,format="%Y-%m-%d")==Apr1[i])[1]]<-winterZ
}

for(i in 1:nBirds){
  model[[i]]<- thresholdModel(twilight = d.twl[[i]]$Twilight,
                              rise = d.twl[[i]]$Rise,
                              twilight.model = "ModifiedLogNormal",
                              alpha = alphaBoth[[i]],
                              beta = beta,
                              # Here is where we set the constraints for land
                              logp.x = log.prior, logp.z = log.prior, 
                              x0 = x0[[i]],
                              z0 = z0[[i]],
                              zenith = zenithBoth[[i]],
                              fixedx = fixedx[[i]])
}


# This defines the error distribution around each location #
proposal.x <- proposal.z <- vector('list',nBirds)

for(i in 1:nBirds){
  proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0[[i]]))
  proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0[[i]]))
}

fit <- vector('list', nBirds)

for(i in 1:nBirds){
  
  fit[[i]] <- estelleMetropolis(model = model[[i]],
                                proposal.x = proposal.x[[i]],
                                proposal.z = proposal.z[[i]],
                                iters = 10000, # This value sets the number of iterations to run
                                thin = 1,
                                chains = 3)

### Fine Tuning 

proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0[[i]]))
proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0[[i]]))

fit[[i]] <- estelleMetropolis(model = model[[i]],
                              proposal.x = proposal.x[[i]],
                              proposal.z = proposal.z[[i]],
                              x0 = chainLast(fit[[i]]$x),
                              z0 = chainLast(fit[[i]]$z),
                              iters=10000, # This value sets the number of iterations to run
                              thin=2,
                              chains=3)

# Final Run

proposal.x[[i]] <- mvnorm(chainCov(fit[[i]]$x),s=0.1)
proposal.z[[i]] <- mvnorm(chainCov(fit[[i]]$z),s=0.1)

# Note the increase in number of interations - this takes a bit longer to run
fit[[i]] <- estelleMetropolis(model = model[[i]],
                              proposal.x = proposal.x[[i]],
                              proposal.z = proposal.z[[i]],
                              x0=chainLast(fit[[i]]$x),
                              z0=chainLast(fit[[i]]$z),
                              iters=5000,  # This value sets the number of iterations to run
                              thin=2,
                              chains=3)

}
## Inital results 

# This step makes an empty raster #
r <- raster(nrows=4*diff(ylim),ncols=4*diff(xlim),xmn=xlim[1],xmx=xlim[2],ymn=ylim[1],ymx=ylim[2])

S <- vector('list',nBirds)

for(i in 1:nBirds){
  S[[i]] <- slices(type="intermediate",
                   breaks="day",
                   mcmc=fit[[i]],
                   grid=r)
}

DATES <-  tm_breed <-tm_winter<-breed<-winter<- vector('list',nBirds)
Aug31 <- Nov01  <- April1<- rep(NA,nBirds)

year<-rep("2011-07-31",20)
year[c(1,3,8)]<-"2010-07-31"

yearAug1<-rep("2011-08-01",20)
yearAug1[c(1,3,8)]<-"2010-08-01"
yearNov1<-rep("2011-11-01",20)
yearNov1[c(1,3,8)]<-"2010-11-01"
yearApril1<-rep("2012-04-01",20)
yearApril1[c(1,3,8)]<-"2011-04-01"

# Create a vector of the Dates in the file - here we set rise=true because the data has both rise and sunset

for(i in 1:nBirds){ 
  DATES[[i]] <- S[[i]]$mcmc[[1]]$time[ which( S[[i]]$mcmc[[1]]$rise==TRUE) ]
  
  
  # Here we specify our dates of interest. 
  ReleaseDay<-1
  Aug31[i]<-which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct(yearAug1[i]))
  
  
  # Breeding 
  # Get time when on the breeding grounds 
  tm_breed[[i]]<-sliceInterval(S[[i]],k=c(ReleaseDay:Aug31[i]))
  
  
  # "Slice" the data and save all dates between Release date and July 31 2011, and June 1 2012 until capture.
  breed[[i]]<-slice(S[[i]],k=c(ReleaseDay:Aug31[i]))
  
  print(BirdId[i])
  
  plot(OVENdist,col="gray74")
  
  plot(breed[[i]],useRaster=TRUE,
       axes=FALSE, add=TRUE,
       legend=FALSE,
       col=rev(bpy.colors(50)),
       cex.axis=0.7)
  
  plot(Americas,border="gray",add=TRUE)
  
  # Non-breeding 
  
  Nov01[i] <- which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct(yearNov1[i]))
  April1[i] <- which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct(yearApril1[i]))
  
  tm_winter[[i]]<- sliceInterval(S[[i]],k=c(Nov01[i]:April1[i]))
  
  # "Slice" data and merge between Nov 01 2014 and Feb 31 2015.
  winter[[i]]<- slice(S[[i]],k=c(Nov01[[i]],April1[[i]]))
  
  #plot(OVENdist,border="black",ylim=ylim,xlim=xlim)
  plot(winter[[i]],useRaster=TRUE,add=TRUE,
       axes=FALSE,
       legend=FALSE,
       col=rev(bpy.colors(50)),
       cex.axis=0.7)
  
  plot(SpatialPoints(cbind(CapLocs[1],CapLocs[2])),add=TRUE,pch=19,cex=1.5)
  box()
} 

par(mar=c(0,0,0,0),mfrow=c(2,2))
for(i in 1:4){
plot(winterIslands,col="lightgray")
plot(mask(winter[[i]],winterIslands),add=TRUE,legend=FALSE,col=rev(bpy.colors(10)))
box()
}



scaledWinter<-scaledBreeding<-vector('list',nBirds)
plot(sum(scaledWinter,na.rm=TRUE))

for(i in 1:nBirds){
scaledBreeding[[i]]<-breed[[i]]/cellStats(breed[[i]],max)
scaledWinter[[i]]<-winter[[i]]/cellStats(winter[[i]],max)
}

scaledBreeding<-stack(scaledBreeding)
scaledWinter<-stack(scaledWinter)

sumBirdsBreed<-sum(scaledBreeding,na.rm=TRUE)
sumBirdsWinter<-sum(scaledWinter,na.rm=TRUE)

sumBirdsBreed[sumBirdsBreed==0]<-NA

E<-sumBirdsBreed/nBirds

prOrigin<-sumBirdsWinter/nBirds

MC_oven<-(cellStats(prOrigin,max) - 1/nBirds) / (1 - 1/nBirds) 
cbind(MC_btbw,MC_oven)

OVENprOrigin<-prOrigin

OVENprOrigin[OVENprOrigin==0]<-NA

OVENpts<-rasterToPoints(OVENprOrigin)

BTBWpts<-rasterToPoints(prOrigin)

plot(winterIslands)
plot(Americas,add=TRUE,col="lightgray")
plot(mask(prOrigin,winterIslands),add=TRUE,legend=FALSE)
plot(Americas,add=TRUE)



############################################################################################################
#
#
# Non-breeding location - Rain
#
#
############################################################################################################
MC_BTBW<-raster("Data/BTBW_MC.grd")
MC_OVEN<-raster("Data/OVEN_MC.grd")

BTBWpts<-rasterToPoints(MC_BTBW)
OVENpts<-rasterToPoints(MC_OVEN)

Year<-1998:2015
nYears<-length(Year)

YearList<-RainRasters<-vector('list',nYears)
for(i in 1:nYears){
  YearList[[i]]<-list.files(paste0("Data/TRMM/",Year[i],"/"),full.names=TRUE)
}

for(i in 1:(nYears-1)){
  temp<-vector('list',12)
  for(m in 1:12){
    temp[[m]]<-raster(YearList[[i]][m])
    temp[[m]][temp[[m]]>10000]<-NA
  }
  RainRasters[[i]]<-stack(temp)
}

temp<-vector('list',8)
for(m in 1:12){
  temp[[m]]<-raster(YearList[[18]][m])
  temp[[m]][temp[[m]]>10000]<-NA
}
RainRasters[[18]]<-stack(temp)

Jan<-Feb<-Mar<-Apr<-May<-Jun<-Jul<-Aug<-Sep<-Oct<-Nov<-Dec<-list()
for(i in 1:nYears){
  Jan[[i]]<-RainRasters[[i]][[1]]
  Feb[[i]]<-RainRasters[[i]][[2]]
  Mar[[i]]<-RainRasters[[i]][[3]]
  Apr[[i]]<-RainRasters[[i]][[4]]
  May[[i]]<-RainRasters[[i]][[5]]
  Jun[[i]]<-RainRasters[[i]][[6]]
  Jul[[i]]<-RainRasters[[i]][[7]]
  Aug[[i]]<-RainRasters[[i]][[8]]
  Sep[[i]]<-RainRasters[[i]][[9]]
  Oct[[i]]<-RainRasters[[i]][[10]]
  Nov[[i]]<-RainRasters[[i]][[11]]
  Dec[[i]]<-RainRasters[[i]][[12]]
} 


for(i in 1:16){
  plot((Oct[[i]]+Nov[[i]]+Dec[[i]]+Jan[[i+1]]+Feb[[i]]+Mar[[i+1]]+Apr[[i+1]])/7)
  Sys.sleep(3)
}

world<-shapefile("Spatial_Layers/WorldAsOne.shp")
Americas<-shapefile("Spatial_Layers/Americas.shp")
PR<-subset(Americas,NAME=="Puerto Rico")
Cuba<-subset(Americas,NAME=="Cuba")
Jamaica<-subset(Americas,NAME=="Jamaica")
Hisp<-subset(Americas,(NAME=="Haiti" | NAME=="Dominican Republic"))
Hisp<-rgeos::gUnaryUnion(Hisp)

rain<-array(NA,c(12,19,4))

months<-c(stack(Jan),stack(Feb),stack(Mar),stack(Apr),stack(May),stack(Jun),
          stack(Jul),stack(Aug),stack(Sep),stack(Oct),stack(Nov),stack(Dec))

sampleYear<-sum(Nov[[1]],Dec[[1]],stack(Jan[[1]],Feb[[1]],Mar[[1]],Apr[[1]]))

CUBelev<-getData('alt',country="CUB")
PRelev<-getData('alt',country="PRI")
DOMelev<-getData('alt',country="DOM")
HAIelev<-getData('alt',country="HTI")
JAMelev<-getData('alt',country="JAM")

WIelev<-merge(DOMelev,HAIelev,JAMelev,CUBelev,PRelev)
WLhs<-hillShade(terrain(WIelev,"slope"),terrain(WIelev,'aspect'),angle=45,degree=10)


winterIslands<-gUnion(gUnion(gUnion(Cuba,Jamaica),Hisp),PR)


PrecipBreaks<-seq(0,1800,)


col<-colorRampPalette(c("transparent", # clear - no color 
                        rgb(238,233,233,alpha=175,maxColorValue =255),
                        rgb(135,206,235,alpha=175,maxColorValue = 255), # skyblue
                        rgb(173,216,230,alpha=175,maxColorValue = 255), # lightblue
                        rgb(0,0,255,alpha=175,maxColorValue = 255),     # blue
                        rgb(75,255,255,alpha=175,maxColorValue = 255), # lightcyan
                        rgb(0,255,255,alpha=175,maxColorValue = 255),     # cyan
                        rgb(0,139,139,alpha=175,maxColorValue = 255),     # darkcyan
                        rgb(255,255,0,alpha=175,maxColorValue = 255),     # yellow
                        rgb(255,173,14,alpha=175,maxColorValue = 255),    # yellowish
                        rgb(255,165,0,alpha=175,maxColorValue = 255),     # orange
                        rgb(255,69,0,alpha=175,maxColorValue = 255),     # orangered
                        rgb(255,0,0,alpha=175,maxColorValue = 255),     # red
                        rgb(139,0,0,alpha=175,maxColorValue = 255),     # darkred
                        rgb(176,48,96,alpha=175,maxColorValue = 255)),   # maroon
                       
alpha=TRUE)
empty<-raster(nrow=1800,ncol=1)
values(empty)<-1:1800

crop(sampleYear,winterIslands)
par(mar=c(4,0,0,0),bty="n")
plot(winterIslands,axes=FALSE)
plot(Americas,add=TRUE,col="gray")
plot(sampleYear,add=TRUE,breaks=PrecipBreaks,col=col(length(PrecipBreaks)),legend=FALSE)
plot(empty,legend.only=TRUE,breaks=PrecipBreaks,col=col(length(PrecipBreaks)),horizontal=TRUE,
     legend.width=0.25, legend.shrink=0.75,
     axis.args=list(at=seq(0,1800, 50),
                    labels=seq(0, 1800, 50), 
                    cex.axis=0.8),
     legend.args=list(text='Precipitation (mm)', side=1,line=1.9,font=2))
plot(Americas,border="black",add=TRUE,lwd=3)

plot(Hisp)
plot(WLhs,add=TRUE)

year<-c(1998:2015)

place<-c(PR,Cuba,Jamaica,Hisp)

totalRain<-SDrain<-array(NA,c(12,18,4))

for(i in 1:12){
  for(p in 1:4){
    
    rain[i,1,p]<-extract(mean(stack(months[[i]])),place[[p]],fun=mean)
    
    for(y in 1:18){
      #standardizedRain
      rain[i,y+1,p]<-extract(months[[i]][[y]],place[[p]],fun=mean)-rain[i,1,p]
      #totalRain
      totalRain[i,y,p]<-extract(months[[i]][[y]],place[[p]],fun=mean)
      SDrain[i,y,p]<-extract(months[[i]][[y]],place[[p]],fun=sd)
    }
  }
}

winterRAINstack<-vector('list',17)
for(i in 2:18){
winterRAINstack[[i-1]]<-stack(Nov[[i-1]],Dec[[i-1]],Jan[[i]],Feb[[i]],Mar[[i]],Apr[[i]])
}

BTBWweightRain<-array(NA,c(nrow(BTBWpts),12,17))
OVENweightRain<-array(NA,c(nrow(OVENpts),12,17))
for(i in 1:17){
BTBWweightRain[,,i]<-extract(winterRAINstack[[i]],cbind(BTBWpts[,1:2]))
OVENweightRain[,,i]<-extract(winterRAINstack[[i]],cbind(OVENpts[,1:2]))
}

WeightedRain<-array(NA,c(17,5,2))
rownames(WeightedRain)<-c(1999:2015)
colnames(WeightedRain)<-c("TotalWinter","SD","MinimumMonth","March","April")
for(i in 1:17){
WeightedRain[i,1,2]<-mean(apply(BTBWweightRain[,c(11,12,1,2,3,4),i],1,sum),weight=BTBWpts[,3],na.rm=TRUE)
WeightedRain[i,2,2]<-mean(apply(BTBWweightRain[,c(11,12,1,2,3,4),i],1,sd),weight=BTBWpts[,3],na.rm=TRUE)
WeightedRain[i,3,2]<-mean(BTBWweightRain[,2,i],weight=BTBWpts[,3],na.rm=TRUE)
WeightedRain[i,4,2]<-mean(BTBWweightRain[,3,i],weight=BTBWpts[,3],na.rm=TRUE)
WeightedRain[i,5,2]<-mean(BTBWweightRain[,4,i],weight=BTBWpts[,3],na.rm=TRUE)
WeightedRain[i,1,1]<-mean(apply(OVENweightRain[,c(11,12,1,2,3,4),i],1,sum,na.rm=TRUE),weight=OVENpts[,3])
WeightedRain[i,2,1]<-mean(apply(OVENweightRain[,c(11,12,1,2,3,4),i],1,sd,na.rm=T),weight=OVENpts[,3],na.rm=T)
WeightedRain[i,3,1]<-mean(OVENweightRain[,2,i],weight=OVENpts[,3],na.rm=TRUE)
WeightedRain[i,4,1]<-mean(OVENweightRain[,3,i],weight=OVENpts[,3],na.rm=TRUE)
WeightedRain[i,5,1]<-mean(OVENweightRain[,4,i],weight=OVENpts[,3],na.rm=TRUE)
}

col2rgb("brown")
plot(WeightedRain[,1,1],pch=19,cex=1.25,type="o",ylim=c(500,900))
polygon(c(1:17,17:1),c(WeightedRain[,1,1]-WeightedRain[,2,1],rev(WeightedRain[,1,1]+WeightedRain[,2,1])),col=rgb(1,0,0,0.5),border="lightgray")
polygon(c(1:17,17:1),c(WeightedRain[,1,2]-WeightedRain[,2,2],rev(WeightedRain[,1,2]+WeightedRain[,2,2])),col=rgb(165,42,42,175,maxColorValue=255),border=rgb(165,42,42,175,maxColorValue=255))
par(new=TRUE)
plot(WeightedRain[,1,2],pch=19,cex=1.25,type="o",ylim=c(500,900),axes=F,ylab="",xlab="")
par(new=TRUE)
plot(OVENwinterWeight,pch=19,cex=1.25,type="o",ylim=c(500,900))

PRmean<-apply(totalRain[,,1],1,mean)
PRse<-apply(totalRain[,,1],1,sd)/sqrt(18)
Hispmean<-apply(totalRain[,,4],1,mean)
Hispse<-apply(totalRain[,,1],1,sd)/sqrt(18)
Cubamean<-apply(totalRain[,,2],1,mean)
Cubase<-apply(totalRain[,,2],1,sd)/sqrt(18)
Jammean<-apply(totalRain[,,3],1,mean)
Jamse<-apply(totalRain[,,3],1,sd)/sqrt(18)

barCenters<-barplot(PRmean[c(10:12,1:9)],
                    ylim=c(0,300),ylab="Precipitation (mm)",yaxt="n",
                    names.arg=c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug","Sept"))
axis(2,las=2)
segments(barCenters,(PRmean[c(10:12,1:9)]-PRse[c(10:12,1:9)]),barCenters,(PRmean[c(10:12,1:9)]+PRse[c(10:12,1:9)]))

barCenters<-barplot(Hispmean[c(10:12,1:9)],
                    ylim=c(0,300),ylab="Precipitation (mm)",yaxt="n",
                    names.arg=c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug","Sept"))
axis(2,las=2)
segments(barCenters,(Hispmean[c(10:12,1:9)]-Hispse[c(10:12,1:9)]),barCenters,(Hispmean[c(10:12,1:9)]+Hispse[c(10:12,1:9)]))

barCenters<-barplot(Cubamean[c(10:12,1:9)],
                    ylim=c(0,300),ylab="Precipitation (mm)",yaxt="n",
                    names.arg=c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug","Sept"))
axis(2,las=2)
segments(barCenters,(Cubamean[c(10:12,1:9)]-Cubase[c(10:12,1:9)]),barCenters,(Cubamean[c(10:12,1:9)]+Cubase[c(10:12,1:9)]))

barCenters<-barplot(Jammean[c(10:12,1:9)],
                    ylim=c(0,300),ylab="Precipitation (mm)",yaxt="n",
                    names.arg=c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug","Sept"))
axis(2,las=2)
segments(barCenters,(Jammean[c(10:12,1:9)]-Jamse[c(10:12,1:9)]),barCenters,(Jammean[c(10:12,1:9)]+Jamse[c(10:12,1:9)]))








# Vector of winter rainfall
# Column 1 = PR
#        2 = Cuba
#        3 = Jamaica
#        4 = Hisp

winterRain<-AprilRain<-MarchRain<-array(NA,c(18,4))
for(i in 1:4){
winterRain[,i]<-apply(totalRain[c(1,2,3,4),1:18,i],2,sum)
AprilRain[,i]<-totalRain[4,1:18,i]
MarchRain[,i]<-totalRain[3,1:18,i]
}

par(bty="l")
plot(winterRain[,1],ylim=c(100,700),type="o",pch=19,cex=1.25,xaxt="n",yaxt="n",ylab="Precipitation (mm)",xlab="Non-breeding season")
par(new=TRUE)
plot(winterRain[,2],ylim=c(100,700),type="o",pch=17,cex=1.25,axes=FALSE,ylab="",xlab="")
par(new=TRUE)
plot(winterRain[,3],ylim=c(100,700),type="o",pch=16,cex=1.25,col="red",axes=FALSE,ylab="",xlab="")
par(new=TRUE)
plot(winterRain[,4],ylim=c(100,700),type="o",pch=8,cex=1.25,col="green",axes=FALSE,xlab="",ylab="")
axis(2,las=2)
axis(1,at=1:18,labels=c(1998:2015))

plot(subset(Americas,(NAME=="Cuba" | NAME=="Jamaica" | NAME=="Puerto Rico" | NAME=="Dominican Republic" | NAME=="Haiti")))


# Read in forest Loss in Km2

PRloss<-readRDS("Data/LossYearsPR.rds")
CubaLoss<-readRDS("Data/LossYearsCuba.rds")
JamaicaLoss<-readRDS("Data/LossYearsJam.rds")
HispLoss<-readRDS("Data/LossYearsHisp.rds")

forestloss<-raster("Data/lossyear.tif")
islandsFL<-crop(forestloss,winterIslands)
minFL<-aggregate(islandsFL,fact=4,fun=modal)
crop(OVENprOrigin,islandsFL)
zonal(islandsFL,crop(OVENprOrigin,islandsFL),fun=min)
plot(minFL)


forestloss
lossOVEN<-extract(forestloss,cbind(OVENpts[,c(1:2)]),buffer=5000)

hist(forestloss)

# Create a Plot x Period x Replicate array of abundance for the species of interest #
# this makes an empty 3 dimensional array - Plot x Period x Replicate that is filled with NAs #

VWdata<-read.csv("Data/Valley-wide_bird_census_1999_2015.csv")

VWdata<-VWdata[VWdata$Distance==1,]

replicate<-1:3

VWdata$PLOT<-as.numeric(as.factor(VWdata$Plot))

# Remove data that needs to be deleted # This was fixed in FileMaker #
VWdata<-VWdata[VWdata$Comments!="**PLEASE DELETE**",]
str(VWdata)
# Some of the species codes have new lines for example "BTBW\n" this next code removes that & turns "BTBW\n" into "BTBW"
VWdata$Species<-gsub("[\r\n]", "", VWdata$Species)

# Some rows of the data were entered by accident with no new species - they have blank species#
# remove those rows from the data set #
VWdata<-VWdata[VWdata$Species != "", ]

# Confirm that they were removed by showing the rows that have blank for species #
# should be 0 rows and it is #
VWdata[which(VWdata$Species ==""),]

#########################################################################################
#########################################################################################
#
#         Set up observation covariates 
#
#########################################################################################
#########################################################################################
# load data for plot site covariates #

PlotAttributes<-read.csv("Data/PlotAttributes.csv")
names(PlotAttributes)
PLOTS<-PlotAttributes$SchwarzID

# Create plot level covariates #
Elev<-(PlotAttributes[,3]-mean(PlotAttributes[,3],na.rm=TRUE))/sd(PlotAttributes[,3],na.rm=TRUE)
Elev2<-(PlotAttributes[,3]^2-mean(PlotAttributes[,3]^2,na.rm=TRUE))/sd(PlotAttributes[,3]^2,na.rm=TRUE)
Slope<-(PlotAttributes[,4]-mean(PlotAttributes[,4],na.rm=TRUE))/sd(PlotAttributes[,4],na.rm=TRUE)
Aspect<-(PlotAttributes[,5]-mean(PlotAttributes[,5],na.rm=TRUE))/sd(PlotAttributes[,5],na.rm=TRUE)


# create an empty array (filled with NA) -   # of rows = sites   # of columns = periods   # of dimensions = replicates
spec.mat<- array(NA, c(length(PLOTS),3,17,2))
Date<-Time<-Obs<-array(NA,c(length(PLOTS),3,17))
rownames(spec.mat)<-PLOTS

# Here is how to fill the array with the data that we want # 
# Here is an example for the Ovenbird #
# OVEN 
S.O.I<-c("OVEN","BTBW")

for(s in 1:2){
newrecord<-subset(VWdata,is.na(New.Record) | New.Record==1)
temp.sp <-subset(newrecord, Species == S.O.I[s]) # SET S.O.I number to SPECIES OF INTEREST
temp.sp$YEAR<-as.numeric(temp.sp$YEAR)

for(y in 1:17){
  temp.yr<- subset(temp.sp, YEAR == y+1998)
# for each replicate #
for(i in 1:3){
  temp.rp <- subset(temp.yr, Replicate == i)
    # for each plot #
    for(k in 1:length(PLOTS)){
      temp.st <- subset(temp.rp, PLOT == PLOTS[k])
      temp.abund <- as.numeric(length(temp.st$Species))
      spec.mat[k,i,y,s]<-ifelse(y %in% c(1:4,7:17), temp.abund,NA) 
      # Ordinal Date 
      Date[k,i,y]<-ifelse(y %in% c(1:4,7:17) & nrow(temp.st)>0, 
                          strptime(as.POSIXct(temp.st$Date,format="%m/%d/%Y"), "%Y-%m-%d")$yday+1,
                          NA)
      # Time of survey
      Time[k,i,y]<-ifelse(y %in% c(1:4,7:17) & nrow(temp.st)>0,
                          format(as.POSIXlt(temp.st$Time,format="%H:%M"),"%H:%M:%S"),
                          NA)
      # Observer
      Obs[k,i,y]<-ifelse(y %in% c(1:4,7:17) & nrow(temp.st)>0,
                         temp.st$Observer,
                         NA)
    }  # For plot
  }  # For period
}# For Replicate
}# For Species

Time<-gsub(Time,pattern=":",replacement="")
Time<-substr(Time,1,4)
Time<-structure(as.numeric(Time), dim=c(373,3,17))
Time<-structure(apply(Time,2,scale),dim=c(373,3,17))

Date<-structure(apply(Date,2,scale),dim=c(373,3,17))
Obs<-structure(apply(Obs,2,scale),dim=c(373,3,17))

Date[is.na(Date)]<-0
Time[is.na(Time)]<-0
Obs[is.na(Obs)]<-0

############ Dail - Madsen simplified ###################
cat("
model {
  ##################################################################################################################
  #
  #  Priors
  #
  ##################################################################################################################
  for(s in 1:nspp){

      gam0[s] ~ dnorm(0, 0.01)
      pInt[s] ~ dnorm(0, 0.01)

  for(k in 1:nyears){

      gam1[k,s] ~ dunif(-10, 10)
      Error[k,s]~dnorm(0,0.01)
    } # nyears

      alpha[s]~dnorm(0,0.01)
    
    for (c in 1:ncovs){ 
      beta[s,c]~dnorm(mu.theta,tau.theta)
    } #ncovs
    
    # detection
    for (m in 1:pcovs){
      betaP[s,m]~dnorm(mu.thetaP,tau.thetaP)  
    } #pcovs
  }    #nspp
  
  ### Hyperpriors #######
  mu.theta ~ dnorm(0,0.01)
  tau.theta ~ dgamma(0.001,0.001)
  mu.thetaP ~ dnorm(0,0.01)
  tau.thetaP ~ dgamma(0.001,0.001)
  
  #################################################################################################################
  #
  #  Likelihood
  #
  #################################################################################################################
  for(s in 1:nspp){                            # Species
    for(i in 1:nschwarz) {                       # Schwarz Plot
      N[i,1,s] ~ dpois(lambda[i,1,s])T(0,)
      lambda[i,1,s]<-exp( alpha[s]+
                          beta[s,1]*Elev[i]+
                          beta[s,2]*Elev2[i]+
                          beta[s,3]*Slope[i]+
                          beta[s,4]*Aspect[i]+
                          Error[1,s]) 

      for(j in 1:nreps) {                        # Replicates
        y[i,j,1,s] ~ dbin(p[i,j,1,s], N[i,1,s]) 
        p[i,j,1,s]<-1/(1+exp(-logit.p[i,j,1,s]))
        logit.p[i,j,1,s]<-pInt[s]+
                          betaP[s,1]*time[i,j,1]+
                          betaP[s,2]*date[i,j,1]#+
                          #betaP[s,3]*obsvr[i,j,1]
      } #REPS
    
      for(k in 2:nyears) {                      # Year 2000-2015
        N[i,k,s] ~ dpois(gamma[i,k-1,s])T(0,)
        gamma[i,k-1,s] <- exp(gam0[s]+
                              gam1[k,s]*N[i,k-1,s]+       #PriorYear&Trend
                                beta[s,1]*Elev[i]+
                                beta[s,2]*Elev2[i]+
                                beta[s,3]*Slope[i]+
                                beta[s,4]*Aspect[i]+
                                Error[k,s] )
        
        for(j in 1:nreps){                      # Replicates
          y[i,j,k,s] ~ dbin(p[i,j,k,s], N[i,k,s]) 
          p[i,j,k,s]<-1/(1+exp(-logit.p[i,j,k,s]))
          logit.p[i,j,k,s]<-pInt[s]+
                            betaP[s,1]*time[i,j,k]+
                            betaP[s,2]*date[i,j,k]#+
                           # betaP[s,3]*obsvr[i,j,k]
          
        }    # REPS
      }      # YEARS
    }        # SITES
  }         # SPECIES
 
  #############################################################################################################################
  #############################################################################################################################
  #
  # Derived parameters 
  #
  #############################################################################################################################
  #############################################################################################################################
  #
  # Calculate the population size sampled at HBEF in each year, k.
  #
  for(s in 1:nspp){
  for(k in 1:nyears){
Ntot[k,s]<-sum(N[,k,s])
  }
#for(i in 1:20844){


  #HBEF[i,1,s]<- exp( alpha[s]+
  #              beta[s,1]*HBEFgridCovs[i,3]+ # Elev
  #              beta[s,2]*pow(HBEFgridCovs[i,3],2)+ #Elev2
  #              beta[s,3]*HBEFgridCovs[i,5]+ #Slope
  #              beta[s,4]*HBEFgridCovs[i,4]+ #Aspect
  #              Error[1,s])

  # for(k in 2:nyears) {

  # Derive the population size within each 50x50m grid cell within HBEF
# HBEF[i,k,s]<-exp(gam0[s]+
  #              gam1[k,s]*HBEF[i,k-1,s]+                                                             
  #              beta[s,1]*HBEFgridCovs[i,3]+ # Elev
  #              beta[s,2]*pow(HBEFgridCovs[i,3],2)+ #Elev2
  #              beta[s,3]*HBEFgridCovs[i,5]+ #Slope
  #              beta[s,4]*HBEFgridCovs[i,4]+ #Aspect
  #              Error[k,s] )
  #  }
 # }
  # determine mean detection for each species
    mean.p[s]<- mean(p[,,,s])
  }

  
  
} # END MODEL",
    fill=TRUE,file="MC_demo.txt")


# Function to put into initial values for Nst #
N<-function(x){
  Nst<-array(NA,dim=c(nrow(x),17,2))
  for(s in 1:2){
    Nst[,,s]<-apply(x[,,,s],c(1,3),max,na.rm=TRUE)+3
  }
  # If Nst==NA change value to 3 #
  Nst[Nst==-Inf]<-NA
  Nst[is.na(Nst)]<-3
  return(Nst)
}

inits<-function() list(
  mu.theta=0,
  tau.theta=10,
  mu.thetaP=0,
  tau.thetaP=10,
  N=N(spec.mat))

nchains<-3
nrow(HBEFgridCovs)

HBEFgridCovs<-read.csv("Data/HBEFgridcovs.csv")

win.data<-list(Elev = Elev,
               Elev2 = Elev2,
               Slope = Slope,
               Aspect = Aspect,
               y = spec.mat,
               nreps=3,
               nschwarz=373,
               nspp=2,
               time = Time,
               date = Date,
               #obsvr = Obs,
               ncovs=4,
               pcovs=2,
               nyears=17,
               HBEFgridCovs = HBEFgridCovs)

params<-c("gam1","Ntot","mean.p","beta","betaP")

library(jagsUI)

a<-Sys.time()
M<-jags(model="MC_demo.txt",
        data= win.data,
        parameters.to.save = params,
        seed = 9328,
        inits = inits,
        n.chain=3,
        n.thin = 2,
        n.iter=1000,
        n.burnin=500,
        DIC = FALSE,
        parallel = TRUE)
Sys.time()-a
Sys.time()

#saveRDS(M,"Recruitment.rds")

M<-readRDS("Data/recruitment.rds")

##########################################################################################################################################
#
# NDVI
#
##############################################################################################################################################


NDVIfiles<-list.files("Data/NDVI",full.names=TRUE)
NDVI<-vector('list',length(NDVIfiles))
for(i in 1:17){
NDVI[[i]]<-raster(NDVIfiles[i])
NDVI[[i]][NDVI[[i]]>1]<-NA
NDVI[[i]]<-crop(NDVI[[i]],winterIslands)
}

NDVIyrs<-sdNDVI<-array(NA,c(17,2))
for(i in 1:17){
NDVIyrs[i,1]<-mean(extract(NDVI[[i]],cbind(OVENpts[,c(1:2)])),weight=OVENpts[,3],na.rm=TRUE)
NDVIyrs[i,2]<-mean(extract(NDVI[[i]],cbind(BTBWpts[,c(1:2)])),weight=BTBWpts[,3],na.rm=TRUE)
sdNDVI[i,1]<-sd(extract(NDVI[[i]],cbind(OVENpts[,c(1:2)])),na.rm=TRUE)
sdNDVI[i,2]<-sd(extract(NDVI[[i]],cbind(BTBWpts[,c(1:2)])),na.rm=TRUE)
}

library(jagsUI)
############ Dail - Madsen fully parameterized ###################

cat("
model {
#####################################################
# Model parameters 
  # Ni1 ~ Poisson(Lambda)
  # Sit ~ Binomial(Nit-1,W)
  # Git ~ Poisson(Nit-1*Gamma)
  # Nit = Sit + Git 
  # yijt ~ Binomial(Nit,p)

  # Nit - latent var., abundance at site i in year t
  # Sit - latent var., survivors
  # Git - latent var., recruits
  # Lambda - mean abundance in year t = 1 
  # W - apparent survival 
  # Gamma - recruitment rate
  # p - detection probability 

######################################################
##################################################################################################################
#
#  Priors
#
##################################################################################################################


#SURVIVAL #
# Priors and constraints
for(t in 1:nyears){
phi.spp[t]~dunif(0,1)
}
p.spp~dunif(0,1)

for (i in 1:n){
     for (t in First[i]:(n.occasions-1)){
          phi[i,t]<-phi.spp[year[t]]
          p.sur[i,t]<-p.spp
         }#t
      }#i

# Likelihood
for (i in 1:n){
        z[i,First[i]]<-1
    for (t in (First[i]+1):n.occasions){
        # state process
        z[i,t]~dbern(mu1[i,t])
        mu1[i,t]<-phi[i,t-1]*z[i,t-1]
        # Observation process
        y.sur[i,t]~dbern(mu2[i,t])
        mu2[i,t]<-p.sur[i,t-1]*z[i,t]
      }#t
     }#i

for(s in 1:nspp){
#W[s] ~ dunif(0.3,0.7)  # Apparent Survival 
#Gamma[s] ~ dnorm(0,0.001) # Recruitment Rate
pInt[s] ~ dnorm(0, 0.01)
alpha[s]~dnorm(0,0.01)
gam0[s]~ dnorm(0,0.01)


    gamMin[s]~dnorm(0,0.01)
    gamMar[s]~dnorm(0,0.01)
    gamApr[s]~dnorm(0,0.01)
    gamTotal[s]~dnorm(0,0.01)

for (c in 1:ncovs){ 
      beta[s,c]~dnorm(0,0.01)
    } #ncovs
    
    # detection
    for (m in 1:pcovs){
      betaP[s,m]~dnorm(0,0.01)  
    } #pcovs

  for(k in 1:nyears){
      Error[k,s]~dnorm(0,0.01)
      gam1[k,s]~dnorm(0,0.01)

   #for(r in 1:4){
  #gamApril[k,s,r]~dnorm(0,0.01)
   #  } # r
    } # nyears
  }    #nspp
#################################################################################################################
#
#  Likelihood
#
#################################################################################################################
# Year 1 - 1999
for(s in 1:nspp){                   # Species
    for(i in 1:nschwarz) {          # Schwarz Plot

    N[i,1,s] ~ dpois(lambda[i,1,s])T(0,)
    lambda[i,1,s]<-exp( alpha[s]+
                        beta[s,1]*Elev[i]+
                        beta[s,2]*Elev2[i]+
                        beta[s,3]*Slope[i]+
                        beta[s,4]*Aspect[i]+
                       Error[1,s]) 
    
  for(j in 1:nreps) {                        # Replicates

    y[i,j,1,s] ~ dbin(p[i,j,1,s], N[i,1,s]) 

    p[i,j,1,s]<-1/(1+exp(-logit.p[i,j,1,s]))

    logit.p[i,j,1,s]<-pInt[s]+
                      betaP[s,1]*time[i,j,1]+
                      betaP[s,2]*date[i,j,1]#+
                      #betaP[s,3]*obsvr[i,j,1]
    } #REPS

for(t in 2:nyears){
S[i,t-1,s] ~ dbin(phi.spp[t-1],N[i,t-1,s])    # Survival

G[i,t-1,s] ~ dpois(Gamma[i,t-1,s]) # Recruits

Gamma[i,t-1,s] <- exp(gam0[s]+gam1[t-1,s]*N[i,t-1,s]+
                              gamTotal[s]*totalRain[t-1,s]+
                              gamMin[s]*minRain[t-1,s]+
                              gamMar[s]*marRain[t-1,s]+
                              gamApr[s]*aprRain[t-1,s])  

N[i,t,s] ~ dpois(Npred[i,t,s]) 
Npred[i,t,s] <- S[i,t-1,s] + G[i,t-1,s]

for(j in 1:nreps) {                        # Replicates

    y[i,j,t,s] ~ dbin(p[i,j,t,s], N[i,t,s]) 

    p[i,j,t,s]<-1/(1+exp(-logit.p[i,j,t,s]))

    logit.p[i,j,t,s]<-pInt[s]+
                      betaP[s,1]*time[i,j,t]+
                      betaP[s,2]*date[i,j,t]#+
                      #betaP[s,3]*obsvr[i,j,t]
          } #REPS
       } # t
    } # i
} # s

#### Derived parameters ####
for(s in 1:2) {
for(t in 1:(nyears-1)) {
#RecruitRate[t,s]<-mean(Gamma[,t,s])
TotRecruit[t,s]<-sum(G[,t,s])
Ntot[t,s]<-sum(N[,t,s])
  } # t
} # s 

######## Derived Parameters #############

} # END MODEL",
    fill=TRUE,file="MC_DM_fullParam_IPM.txt")



# Read in annual capture history for adults (M and F) 
AnnualHistory<-read.csv("D:/Google_Drive/BTBW_HBEF/Cost_Of_Repro/Data/CaptureHistory_86_15.csv",header=TRUE)

# Convert capture history into a matrix - ch
ch<-as.matrix(AnnualHistory[,23:39])

# check if any cells have NA
any(is.na(ch)) # returns FALSE 

# sum columns - if 0 inds were not seen #
above0<-which(apply(ch,1,sum)>=1)

ch<-ch[which(apply(ch,1,sum)>=1),]

sex<-AnnualHistory$Sex[above0]

bandNum<-AnnualHistory$AlumBand[above0]

YearBand<-AnnualHistory$FirstYear[above0]

# assume ? and blanks are females #
S<-rep(NA,dim(ch)[1])
for(i in 1:dim(ch)[1]){
if(sex[i]=="F") S[i]<-2
if(sex[i]=="M") S[i]<-1
}

S[is.na(S)]<-2

n.occasions<-dim(ch)[2]
N.inds<-dim(ch)[1]

#remove any individuals that were not seen during that time frame 
## Create a function to find the when the individual was first marked / seen ##
get.first<-function(x) min(which(x!=0))
get.last<-function(x) max(which(x!=0))

First<-apply(ch,1,get.first)
Last<-apply(ch,1,get.last)


# Function to put into initial values for Nst #
N<-function(x){
  Nst<-array(NA,dim=c(nrow(x),17,2))
  for(s in 1:2){
    Nst[,,s]<-apply(x[,,,s],c(1,3),max,na.rm=TRUE)+1
  }
  # If Nst==NA change value to 3 #
  Nst[Nst==-Inf]<-NA
  Nst[is.na(Nst)]<-1
  return(Nst)
}


# Function to change leading zeros before banding to NA #
ch.init<-function(ch,First){
          for (i in 1:dim(ch)[1]){
             ch[i,1:First[i]]<-NA}
              return(ch)
}


inits<-function() list(N=N(spec.mat),
                       z=ch.init(matrix(1,dim(ch)[1],dim(ch)[2]),First),
                       phi.spp=runif(17,0,1),
                       p.spp=runif(1,0,1))


nchains<-3

win.data<-list(Elev = Elev,
               Elev2 = Elev2,
               Slope = Slope,
               Aspect = Aspect,
               y = spec.mat,
               nreps=3,
               nschwarz=373,
               nspp=2,
               time = Time,
               date = Date,
               #obsvr = Obs,
               ncovs=4,
               pcovs=2,
               nyears=17,
               minRain = WeightedRain[,3,],
               totalRain = WeightedRain[,1,],
               marRain = WeightedRain[,4,],
               aprRain = WeightedRain[,5,],
               y.sur=ch,
               First=First,
               n=dim(ch)[1],
               n.occasions=dim(ch)[2],
               sex=S,
                year=c(1:18))

               #HBEFgridCovs = HBEFgridCovs)

params<-c("gam0","gam1","TotRecruit","gamMin","gamMar","gamApr","gamTotal","phi.spp","p.spp","G","S","Ntot","phi.spp")

library(jagsUI)

Sys.time()
a<-Sys.time()
M<-jags(model="MC_DM_fullParam_IPM.txt",
        data= win.data,
        parameters.to.save = params,
        seed = 9328,
        inits = inits,
        n.chain=3,
        n.thin = 1,
        n.iter=2500,
        codaOnly = c("Gamma","p.spp","phi.spp","TotRecruit","gam1","G","S","Ntot"),
        n.burnin=500,
        DIC = FALSE,
        parallel = TRUE)
Sys.time()-a
Sys.time()
#saveRDS(M,"FullyParamRainResults.rds")


1999:2015
str(M$sims.list$TotRecruit)
Recruits<-apply(M$sims.list$TotRecruit,c(2,3),mean)
q97.5<-array(NA,c(16,2,2))
for(i in 1:16){
for(s in 1:2){
q97.5[i,1,s]<-quantile(M$sims.list$TotRecruit[,i,s],probs=0.975)
q97.5[i,2,s]<-quantile(M$sims.list$TotRecruit[,i,s],probs=0.025)
}
}
OVENr<-Recruits[,1]
par(bty='l')
plot(Recruits[,1]~WeightedRain[2:17,5,1],pch=19,cex=2,ylab="Recruits",xlab="April Precipitation (mm)",ylim=c(0,2000),yaxt="n")
segments(WeightedRain[2:17,5,1],q97.5[,1,1],WeightedRain[2:17,5,1],q97.5[,2,1])
text(WeightedRain[2:17,5,1]+1,Recruits[,1]-50,labels=paste0("'",substr(rownames(WeightedRain)[2:17],3,4)))
axis(2,las=2)

plot(Recruits[,2]~WeightedRain[2:17,5,2],pch=19,cex=2,ylab="Recruits",xlab="April Precipitation (mm)",ylim=c(0,3000),yaxt="n")
segments(WeightedRain[2:17,5,2],q97.5[,1,2],WeightedRain[2:17,5,2],q97.5[,2,2])
text(WeightedRain[2:17,5,2]+3,Recruits[,2]-100,labels=paste0("'",substr(rownames(WeightedRain)[2:17],3,4)))
axis(2,las=2)


############################################################################################################################
#
#
#  Plot-level 
#
#
###########################################################################################################################

# OVENBIRD 
OVENdata<-read.csv("Data/OvenCaptures.csv",header=FALSE)
names(OVENdata)<-c("AlumBand","Sex","Date","Age","L","R","ColorBands")

OVENdata<-subset(OVENdata,Age=="ASY" | Age=="SY")
OVENdata$Date<-as.character.Date(OVENdata$Date)
OVENdata$Date<-as.POSIXct(OVENdata$Date,format="%m/%d/%Y")
OVENdata$Year<-format(OVENdata$Date,"%Y")

yr<-c(2009,2010,2011,2012,2013,2014)
OVENpropSY<-rep(NA,length(yr))
for(i in 1:length(yr)){
OVENpropSY[i]<-table(subset(OVENdata$Age,OVENdata$Year==yr[i]))[7]/sum(table(subset(OVENdata$Age,OVENdata$Year==yr[i]))[c(4,7)])  
}

# Black-throated Blue Warbler
BTBWnests<-shapefile("C:/Users/Michael/Dropbox (Smithsonian)/BTBW Project/Spatial_Layers/NestShapefiles/NestLocations_1986_2015.shp")
names(BTBWnests)
BTBWnestData<-BTBWnests@data
BTBWnestData$fledged<-as.numeric(BTBWnestData$fledged)
library(dplyr)

Nestlings<-group_by(BTBWnestData,Aluminum_F,Year)%>%
             summarize(Young=sum(fledged,na.rm=TRUE))

Fecund<-group_by(Nestlings,Year)%>%
             summarize(Fecund=mean(Young,na.rm=TRUE))

Fecund<-as.data.frame(Fecund)

BTBWdata<-read.csv("Data/BtbwCaptures.csv")

BTBWdata<-subset(BTBWdata,AgeInital=="SY" | AgeInital=="ASY")
BTBWdata$Date<-as.POSIXct(BTBWdata$EarliestDate,format="%m/%d/%Y")
BTBWdata$Year<-format(BTBWdata$Date,"%Y")

BTBWdata<-subset(BTBWdata,Year>=1999)

BTBWyr<-1999:2015
BTBWpropSY<-BTBWrecruits<-SYs<-rep(NA,length(BTBWyr))
for(i in 1:length(BTBWyr)){
  BTBWpropSY[i]<-table(subset(BTBWdata$AgeInital,BTBWdata$Year==BTBWyr[i]))[9]/sum(table(subset(BTBWdata$AgeInital,BTBWdata$Year==BTBWyr[i]))[c(4,9)])  
  SYs[i]<-table(subset(BTBWdata$AgeInital,BTBWdata$Year==BTBWyr[i]))[9]
  BTBWrecruits[i]<-sum(table(subset(BTBWdata$AgeInital,BTBWdata$Year==BTBWyr[i])))
}

years<-1999:2015
plotEffort<-vector('list',length(years))
effort<-rep(NA,length(years))
for(i in 1:length(years)){
plotEffort[[i]]<-shapefile(paste0("C:/Users/Michael/Dropbox (Smithsonian)/BTBW Project/Spatial_Layers/PlotEffort/BoundingPlotBoxes/BB_Territories",years[i],".shp"))
effort[i]<-gArea(plotEffort[[i]])/10000
}

plot(BTBWrecruits/effort~Fecund[15:31,2])
cbind(1999:2015,SYs/effort,BTBWrecruits/effort,Fecund[15:31,])


linearModel<-lm((BTBWrecruits/effort)~Fecund[15:31,2])
linearModel$residuals
plot((linearModel$residuals~NDVIyrs[,2]),pch=19,cex=2,yaxt="n",ylab="Residuals",xlab="NDVI")
text(NDVIyrs[,2],linearModel$residuals-0.01,label=paste0("'",substr(BTBWyr,3,4)))
axis(2,las=2)
abline(lm(linearModel$residuals[c(1:3,5:17)]~NDVIyrs[c(1:3,5:17),2]))
abline(lm(linearModel$residuals~NDVIyrs[,2]),lty=2)

summary(lm(linearModel$residuals[c(1:3,5:17)]~NDVIyrs[c(1:3,5:17),2]))
#points(loess.smooth(NDVIyrs[,2],linearModel$residuals),type="l")

par(bty="l")
plot((BTBWrecruits/effort)~NDVIyrs[,2],pch=19,cex=2,yaxt="n",ylab="New recruits / ha",xlab="Weighted NDVI",xlim=c(0.46,0.56),ylim=c(0,0.8))
axis(2,las=2)
abline(lm((BTBWrecruits/effort)~NDVIyrs[,2]))
text(NDVIyrs[,2],(BTBWrecruits/effort)-0.05,label=paste0("'",substr(BTBWyr,3,4)))



plot(BTBWpropSY~NDVIyrs[,1])
abline(lm(BTBWpropSY~NDVIyrs[,1]))


plot(apply(M$sims.list$TotRecruit,c(2,3),mean)[3:16,1]~PRloss,pch=19,cex=2,ylab="",xlab="Forest Loss (km2)",yaxt="n")
axis(2,las=2)
abline(lm(apply(M$sims.list$TotRecruit,c(2,3),mean)[3:16,2]~JamaicaLoss))


plot(apply(M$sims.list$TotRecruit,c(2,3),mean)[3:16,1]~HispLoss,pch=19,cex=2,ylab="",xlab="Forest Loss (km2)",yaxt="n")
axis(2,las=2)
abline(lm(apply(M$sims.list$TotRecruit,c(2,3),mean)[3:16,1]~HispLoss))

par(bty="l")
plot(apply(M$sims.list$TotRecruit,c(2,3),mean)[3:16,2]~NDVIyrs[2:17,2],pch=19,cex=2,ylab="",xlab="",yaxt="n")
axis(2,las=2)
##########################################################################################################################################
#
# NDVI
#
##############################################################################################################################################
M$mean$gamApr

NDVIfiles<-list.files("Data/NDVI",full.names=TRUE)
NDVI<-vector('list',length(NDVIfiles))
for(i in 1:17){
NDVI[[i]]<-raster(NDVIfiles[i])
NDVI[[i]][NDVI[[i]]>1]<-NA
NDVI[[i]]<-crop(NDVI[[i]],winterIslands)
}

NDVIyrs<-sdNDVI<-array(NA,c(17,2))
for(i in 1:17){
NDVIyrs[i,1]<-mean(extract(NDVI[[i]],cbind(OVENpts[,c(1:2)])),weight=OVENpts[,3],na.rm=TRUE)
NDVIyrs[i,2]<-mean(extract(NDVI[[i]],cbind(BTBWpts[,c(1:2)])),weight=BTBWpts[,3],na.rm=TRUE)
sdNDVI[i,1]<-sd(extract(NDVI[[i]],cbind(OVENpts[,c(1:2)])),na.rm=TRUE)
sdNDVI[i,2]<-sd(extract(NDVI[[i]],cbind(BTBWpts[,c(1:2)])),na.rm=TRUE)
}

plot(Recruits[,2]~NDVIyrs[2:17,2])

quantile(M$sims.list$phi.spp[,1],probs=0.975)
par(bty="l")
plot(M$mean$phi.spp~NDVIyrs[,2],pch=19,cex=2,ylab="Survival",xlab="Weighted NDVI",yaxt="n",ylim=c(0,1))
segments(NDVIyrs[,2],M$q2.5$phi.spp,NDVIyrs[,2],M$q97.5$phi.spp)
axis(2,las=2)
text(NDVIyrs[,2],M$mean$phi-0.03,label=paste0("'",substr(BTBWyr,3,4)))
abline(lm(M$mean$phi.spp~NDVIyrs[,2]))

visNDVI<-raster(NDVIfiles[17])
plot(visNDVI)
set.breaks<-seq(0,1,,100)
cols<-colorRampPalette(c("brown","beige","forestgreen"))
plot(winterIslands)
plot(visNDVI,add=TRUE)