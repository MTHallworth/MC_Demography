library(raster)
library(rgdal)
library(sp)
library(RColorBrewer)
library(maptools)
library(rgeos)
library(geosphere)
library(SGAT)
library(BAStag)
library(MASS)

##############################################################################
#
#
# Light-Level Geolocator Analysis 
#
#
#############################################################################

# Black-throated blue warbler 

Americas<-shapefile("Spatial_Layers/Americas.shp")
states<-shapefile("Spatial_Layers/st99_d00.shp")
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
BTBWdata<-vector('list',nBirds)

# Loop through all the files and read them in as LUX files #
for(i in 1:nBirds){
  BTBWdata[[i]] <- readLig(BTBW_HBEF[i],skip=1)  
}

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


# Defining Twilights
tm<-rise<-vector('list',nBirds)

for(i in 1:nBirds){
  tm[[i]] <- seq(from = BTBWdata[[i]][1,2], 
                 to = BTBWdata[[i]][nrow(BTBWdata[[i]]),2], 
                 by = "day")
  
  rise[[i]] <- rep(c(TRUE, FALSE), length(tm[[i]]))
}

# making predicted twilight times given location and zenith #
cal.dat<-vector('list',nBirds)

for(i in 1:nBirds){
  cal.dat[[i]] <- data.frame(Twilight = twilight(rep(tm[[i]], each = 2),
                                                 lon = CapLocs[i,1], 
                                                 lat = CapLocs[i,2], 
                                                 rise = rise[[i]], zenith = 94),
                             Rise = rise[[i]]) 
}

twl<-vector('list',nBirds)
for(i in 1:nBirds){
twl[[i]]<-readRDS(paste0("Data/GLdata/BTBW/twilightFiles/",BirdId[[i]],".rds"))
}

# Create a vector with the dates known to be at deployment #
calib.dates <- vector('list',nBirds)

for(i in 1:nBirds){
  calib.dates[[i]] <- c(strptime(twl[[i]][1,1],format="%Y-%m-%d"),as.POSIXct("2015-08-15"))
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
b[b>400]<-NA

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

# subset the twilight file for dates after the first calibration date (presumably the deployment date)  
# and exclude points that were deleted  
# note we didn't delete any transitions here

for(i in 1:nBirds){  
  d.twl[[i]]<-subset(twl[[i]],twl[[i]]$Twilight>=calib.dates[[i]][1] & !Deleted)
  
  path[[i]] <- thresholdPath(twilight = d.twl[[i]]$Twilight,
                             rise = d.twl[[i]]$Rise,
                             zenith = zenith0[i],
                             tol = 0)
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
  
  fixedx[[i]][c(1:10,(nrow(x0[[i]])-3):nrow(x0[[i]]))] <- TRUE
  
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
  model[[i]]<- thresholdModel(twilight = d.twl[[i]]$Twilight,
                              rise = d.twl[[i]]$Rise,
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




#########################################################################
#   OVENBIRDS
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

# Read in the lux file
OVENdata<-vector('list',nBirds)

# Loop through all the files and read them in as LUX files #
for(i in 1:nBirds){
  OVENdata[[i]] <- readLig(OVEN_HBEF[i],skip=1)  
}

head(OVENdata[[1]])

# Set the capture coordinates for each bird #
CapLocs<-array(NA,c(nBirds,2))

# give row names the same as the OVEN names files to keep everything organized.
rownames(CapLocs)<-OVENnames

# Capture locations
CapLocs<- cbind(-71.73,43.94)

# Defining Twilights
tm<-rise<-vector('list',nBirds)

for(i in 1:nBirds){
  tm[[i]] <- seq(from = OVENdata[[i]][1,2], 
                 to = OVENdata[[i]][nrow(OVENdata[[i]]),2], 
                 by = "day")
  
  rise[[i]] <- rep(c(TRUE, FALSE), length(tm[[i]]))
}

# making predicted twilight times given location and zenith #
cal.dat<-vector('list',nBirds)

for(i in 1:nBirds){
  cal.dat[[i]] <- data.frame(Twilight = twilight(rep(tm[[i]], each = 2),
                                                 lon = CapLocs[1], 
                                                 lat = CapLocs[2], 
                                                 rise = rise[[i]], zenith = 94),
                             Rise = rise[[i]]) 
}

twl<-vector('list',nBirds)
for(i in 1:nBirds){
twl[[i]]<-readRDS(paste0("Data/GLdata/OVEN/twilightFiles/","NH_",BirdId[[i]],".rds"))
}

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
CapLocs_NB<-array(NA,c(16,2))

twl_NB<-vector('list',16)
twl_rds<-list.files(path = "D:/GL_wind",pattern="twl_",full.names=TRUE)
for(i in 1:16){
twl_NB[[i]]<-readRDS(twl_rds[i])
}

# Create an array for the long / lat where birds were captured
CapLocs_NB[4:12,1]<- -77.93
CapLocs_NB[4:12,2]<- 18.04
CapLocs_NB[1:3,1]<- -80.94
CapLocs_NB[1:3,2]<- 25.13
CapLocs_NB[13:16,1]<- -66.86
CapLocs_NB[13:16,2]<- 17.97

plot(OVENdist,col="gray",border="gray")
plot(Americas,add=TRUE)
points(CapLocs_NB,cex=1.5,pch=19)


head(twl_NB[[13]])


end.date_NB<-c(rep("2011-04-15",3),
                   "2011-04-15",
               rep("2010-04-15",2),
               rep("2011-04-15",4),
                   "2010-04-15",
                   "2011-04-15",
               rep("2012-04-15",4))
               
               

# Create a vector with the dates known to be at deployment #
calib.dates_NB <- vector('list',16)

for(i in 1:16){
  calib.dates_NB[[i]] <- c(strptime(twl_NB[[i]][1,1],format="%Y-%m-%d"),as.POSIXct(end.date_NB[i]))
}

calibration.data_NB<-vector('list',16)

for(i in 1:16){
  calibration.data_NB[[i]]<-subset(twl_NB[[i]],twl_NB[[i]]$Twilight>=calib.dates_NB[[i]][1] & twl_NB[[i]]$Twilight<=calib.dates_NB[[i]][2])
}



# Generate empty lists to store data 
sun_NB<-z_NB<-twl_t_NB<-twl_deviation_NB<-fitml_NB<-alpha_NB<-vector('list',16)

# Determine the sun elevation angle - here called the Zenith angle #

for(i in 1:16){
  
  # Calculate solar time from calibration data 
  sun_NB[[i]]  <- solar(calibration.data_NB[[i]][,1])
  
  # Adjust the solar zenith angle for atmospheric refraction
  z_NB[[i]] <- refracted( zenith(sun = sun_NB[[i]],
                              lon = CapLocs_NB[i,1], 
                              lat = CapLocs_NB[i,2]))
  
  twl_t_NB[[i]] <- twilight(tm = calibration.data_NB[[i]][,1],
                         lon = CapLocs_NB[i,1], 
                         lat = CapLocs_NB[i,2], 
                         rise = calibration.data_NB[[i]][,2],
                         zenith = quantile(z_NB[[i]],probs=0.5))
  
  # Determine the difference in minutes from when the sun rose and the geolocator said it rose 
  twl_deviation_NB[[i]] <- ifelse(calibration.data_NB[[i]]$Rise, as.numeric(difftime(calibration.data_NB[[i]][,1], twl_t_NB[[i]], units = "mins")),
                               as.numeric(difftime(twl_t_NB[[i]], calibration.data_NB[[i]][,1], units = "mins")))
  
  # Throw out values less than 0 - These values are not valid 
  twl_deviation_NB[[i]]<-subset(twl_deviation_NB[[i]], subset=twl_deviation_NB[[i]]>=0)
  
  # Describe the distribution of the error 
  fitml_NB[[i]] <- fitdistr(twl_deviation_NB[[i]], "log-Normal")
  
  # save the Twilight model parameters
  alpha_NB[[i]]<- c(fitml_NB[[i]]$estimate[1], fitml_NB[[i]]$estimate[2]) 
}

alpha
alpha_NB

meanBalpha<-c(mean(unlist(lapply(alpha,function(x) x[1]))),mean(unlist(lapply(alpha,function(x) x[2]))))
meanNBalpha<-c(mean(unlist(lapply(alpha_NB,function(x) x[1]))),mean(unlist(lapply(alpha_NB,function(x) x[2]))))
names(meanBalpha)<-names(meanNBalpha)<-names(alpha[[1]])



winterZ<-quantile(unlist(z_NB),probs=0.5)
print("Breeding Zenith Angle")
median(unlist(z))
print("Non-breeding Zenith Angle")
median(unlist(z_NB))


b<-unlist(twl_deviation)
b[b>400]<-NA

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
plot(median(z[[1]],na.rm=TRUE),xlim=c(1,nBirds),ylim=c(80,100),pch=19,ylab="Zenith Angle",xlab="OVEN",col=cols[1])
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

# subset the twilight file for dates after the first calibration date (presumably the deployment date)  
# and exclude points that were deleted  
# note we didn't delete any transitions here

for(i in 1:nBirds){  
  d.twl[[i]]<-subset(twl[[i]],twl[[i]]$Twilight>=calib.dates[[i]][1] & !Deleted)
  
  path[[i]] <- thresholdPath(twilight = d.twl[[i]]$Twilight,
                             rise = d.twl[[i]]$Rise,
                             zenith = zenith1[i],
                             tol = 0)
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


scaledWinter<-scaledBreeding<-vector('list',nBirds)
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

prOrigin<-sumBirdsWinter/4

MC_btbw<-(cellStats(prOrigin,max) - 1/nBirds) / (1 - 1/nBirds) 











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

year<-c(1998:2015)

place<-c(PR,Cuba,Jamaica,Hisp)

totalRain<-array(NA,c(12,18,4))

for(i in 1:12){
  for(p in 1:4){
    
    rain[i,1,p]<-extract(mean(stack(months[[i]])),place[[p]],fun=mean)
    
    for(y in 1:18){
      #standardizedRain
      rain[i,y+1,p]<-extract(months[[i]][[y]],place[[p]],fun=mean)-rain[i,1,p]
      #totalRain
      totalRain[i,y,p]<-extract(months[[i]][[y]],place[[p]],fun=mean)
    }
  }
}



# Vector of winter rainfall
# Column 1 = PR
#        2 = Cuba
#        3 = Jamaica
#        4 = Hisp

winterRain<-AprilRain<-MarchRain<-array(NA,c(18,4))
for(i in 1:4){
winterRain[,i]<-apply(rain[c(1,2,3,4),2:19,i],2,sum)
AprilRain[,i]<-totalRain[4,1:18,i]
MarchRain[,i]<-totalRain[3,1:18,i]
}



# Read in forest Loss in Km2

PRloss<-readRDS("Data/LossYearsPR.rds")
CubaLoss<-readRDS("Data/LossYearsCuba.rds")
JamaicaLoss<-readRDS("Data/LossYearsJam.rds")
HispLoss<-readRDS("Data/LossYearsHisp.rds")

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
for(s in 1:nspp){

W[s] ~ dunif(0.3,0.7)  # Apparent Survival 
#Gamma[s] ~ dnorm(0,0.001) # Recruitment Rate
pInt[s] ~ dnorm(0, 0.01)
alpha[s]~dnorm(0,0.01)
gam0[s]~ dnorm(0,0.01)


for(r in 1:4){
    gamApril[s,r]~dnorm(0,0.01)
    gamTotal[s,r]~dnorm(0,0.01)
    } # r

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
S[i,t-1,s] ~ dbin(W[s],N[i,t-1,s])    # Survival

G[i,t-1,s] ~ dpois(Gamma[i,t-1,s]) # Recruits

Gamma[i,t-1,s] <- exp(gam0[s]+gam1[t-1,s]*N[i,t-1,s]+
                              gamApril[s,1]*AprilRain[t-1,1]+ # PuertoRico
                              gamApril[s,2]*AprilRain[t-1,2]+ # Cuba
                              gamApril[s,3]*AprilRain[t-1,3]+ # Jamaica
                              gamApril[s,4]*AprilRain[t-1,4]+ # Hisp
                              gamTotal[s,1]*totalRain[t-1,1]+ # PuertoRico
                              gamTotal[s,2]*totalRain[t-1,2]+ # Cuba
                              gamTotal[s,3]*totalRain[t-1,3]+ # Jamaica
                              gamTotal[s,4]*totalRain[t-1,4])  # Hisp
  

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
RecruitRate[t,s]<-mean(Gamma[,t,s])
TotRecruit[t,s]<-sum(G[,t,s])
  } # t
} # s 

######## Derived Parameters #############

} # END MODEL",
    fill=TRUE,file="MC_DM_fullParam.txt")

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


inits<-function() list(N=N(spec.mat),
                       W = c(0.5,0.5))


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
               AprilRain = AprilRain[2:18,],
               totalRain = winterRain[2:18,])

               #HBEFgridCovs = HBEFgridCovs)

params<-c("W","Gamma","RecruitRate","gam1","TotRecruit","gamApril","gamTotal")

library(jagsUI)

Sys.time()
a<-Sys.time()
M<-jags(model="MC_DM_fullParam.txt",
        data= win.data,
        parameters.to.save = params,
        seed = 9328,
        inits = inits,
        n.chain=3,
        n.thin = 5,
        n.iter=10000,
        codaOnly = c("Gamma"),
        n.burnin=5000,
        DIC = FALSE,
        parallel = TRUE)
Sys.time()-a
Sys.time()

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
BTBWdata<-read.csv("Data/BtbwCaptures.csv")

BTBWdata<-subset(BTBWdata,AgeInital=="SY" | AgeInital=="ASY")
BTBWdata$Date<-as.POSIXct(BTBWdata$EarliestDate,format="%m/%d/%Y")
BTBWdata$Year<-format(BTBWdata$Date,"%Y")

BTBWdata<-subset(BTBWdata,Year>=1999)

BTBWyr<-1999:2015
BTBWpropSY<-rep(NA,length(BTBWyr))
for(i in 1:length(BTBWyr)){
  BTBWpropSY[i]<-table(subset(BTBWdata$AgeInital,BTBWdata$Year==BTBWyr[i]))[9]/sum(table(subset(BTBWdata$AgeInital,BTBWdata$Year==BTBWyr[i]))[c(4,9)])  
}

par(bty="l",mfrow=c(2,2))
plot(BTBWpropSY~MarchRain[2:18,1],pch=19,ylim=c(0,1),main="Puerto Rico",xlab="Precipitation")
text(x=MarchRain[2:18,1],y=BTBWpropSY-0.025,labels=BTBWyr,cex=0.7)
plot(BTBWpropSY~MarchRain[2:18,2],pch=19,ylim=c(0,1),main="Cuba",xlab="Precipitation")
text(x=MarchRain[2:18,2],y=BTBWpropSY-0.025,labels=BTBWyr,cex=0.7)
plot(BTBWpropSY~MarchRain[2:18,3],pch=19,ylim=c(0,1),main="Jamaica",xlab="Precipitation")
text(x=MarchRain[2:18,3],y=BTBWpropSY-0.025,labels=BTBWyr,cex=0.7)
plot(BTBWpropSY~MarchRain[2:18,4],pch=19,ylim=c(0,1),main="Hispainola",xlab="Precipitation")
text(x=MarchRain[2:18,4],y=BTBWpropSY-0.025,labels=BTBWyr,cex=0.7)