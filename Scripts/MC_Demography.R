library(raster)
library(rgdal)

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

plot(rain[1,2:19,1],pch=19,type="o")

par(bty="l")
plot(apply(rain[c(1,2,3,4),2:19,1],2,sum),pch=19,type="o",ylab="Centered Rainfall",xlab="Year",yaxt="n",xaxt="n",ylim=c(-300,300))
par(new=TRUE)
plot(apply(rain[c(1,2,3,4),2:19,2],2,sum),pch=19,type="o",col="green",ylab="",xlab="Year",yaxt="n",xaxt="n",ylim=c(-300,300))
par(new=TRUE)
plot(apply(rain[c(1,2,3,4),2:19,3],2,sum),pch=19,type="o",col="red",ylab="",xlab="Year",yaxt="n",xaxt="n",ylim=c(-300,300))
par(new=TRUE)
plot(apply(rain[c(1,2,3,4),2:19,4],2,sum),pch=19,type="o",col="blue",ylab="",xlab="Year",yaxt="n",xaxt="n",ylim=c(-300,300))
abline(h=0,lty=2,col="gray")
axis(2,las=2)
axis(1,at=c(1:18),year)

# Vector of winter rainfall
# Column 1 = PR
#        2 = Cuba
#        3 = Jamaica
#        4 = Hisp

winterRain<-AprilRain<-array(NA,c(18,4))
for(i in 1:4){
winterRain[,i]<-apply(rain[c(1,2,3,4),2:19,i],2,sum)
AprilRain[,i]<-rain[4,2:19,i]
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

      gam1[k,s] ~ dnorm(0, 0.01)
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
      N[i,1,s] ~ dpois(lambda[i,1,s])
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
      
      
      for(k in 2:4) {                      # Year 2000-2002
        N[i,k,s] ~ dpois(gamma[i,k-1,s])
        gamma[i,k-1,s] <- exp(gam0[s]+
                              gam1[k,s]*N[i,k-1,s]+   #PriorYear&Trend
                              beta[s,1]*Elev[i]+
                              beta[s,2]*Elev2[i]+
                              beta[s,3]*Slope[i]+
                              beta[s,4]*Aspect[i]+
                              Error[k,s])
        
        for(j in 1:nreps){                      # Replicates
          y[i,j,k,s] ~ dbin(p[i,j,k,s], N[i,k,s]) 
          p[i,j,k,s]<-1/(1+exp(-logit.p[i,j,k,s]))
          logit.p[i,j,k,s]<-pInt[s]+
                            betaP[s,1]*time[i,j,k]+
                            betaP[s,2]*date[i,j,k]#+
                           # betaP[s,3]*obsvr[i,j,k]
          
        }    # REPS
      }      # YEARS
      
      for(k in 5:6) {                      # Year 2003-2004 - years with no counts
        N[i,k,s] ~ dpois(gamma[i,k-1,s])
        gamma[i,k-1,s] <- exp(gam0[s]+
                              gam1[k,s]*mean(N[i,k-1,s])+    #PriorYear&Trend
                              beta[s,1]*Elev[i]+
                              beta[s,2]*Elev2[i]+
                              beta[s,3]*Slope[i]+
                              beta[s,4]*Aspect[i]+
                              Error[k,s])
        
        for(j in 1:nreps){                      # Replicates
          y[i,j,k,s] ~ dbin(p[i,j,k,s], N[i,k,s]) 
          p[i,j,k,s]<-1/(1+exp(-logit.p[i,j,k,s]))
          logit.p[i,j,k,s]<-pInt[s]+
                            betaP[s,1]*time[i,j,k]+
                            betaP[s,2]*date[i,j,k]#+
                            #betaP[s,3]*obsvr[i,j,k]
          
        }    # REPS
      }      # YEARS
      
      
      
      for(k in 7:nyears) {                      # Year 2005-2014
        N[i,k,s] ~ dpois(gamma[i,k-1,s])
        gamma[i,k-1,s] <- exp(gam0[s]+
                              gam1[k,s]*N[i,k-1,s]+                                                             #PriorYear&Trend
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

M<-jags(model="MC_demo.txt",
        data= win.data,
        parameters.to.save = params,
        inits = inits,
        n.chain=3,
        n.thin = 2,
        n.iter=2000,
        n.burnin=500,
        parallel = TRUE)


library(raster)
plot(rasterFromXYZ(cbind(HBEFgridCovs[,1:2],M$mean$HBEF[,17,2])))

nbLoc<-c("Puerto Rico","Cuba","Jamaica","Hispaniola")

par(mfrow=c(2,4))
for(i in 1:4){
    par(bty="l")
   if(i ==1){ plot(apply(totalRain[1:4,2:18,i],2,sum),M$mean$gam1[1:17,1],pch=19,cex=2,
         main=paste(nbLoc[i]),xlab="", ylab="Recruitment")}
  else{plot(apply(totalRain[1:4,2:18,i],2,sum),M$mean$gam1[1:17,1],pch=19,cex=2,
              main=paste(nbLoc[i]),xlab="", ylab="")}
    
    abline(lm(M$mean$gam1[,1]~apply(totalRain[1:4,2:18,i],2,sum)))
}
for(i in 1:4){
  par(bty="l")
  if(i ==1){ plot(apply(totalRain[1:4,2:18,i],2,sum),M$mean$gam1[1:17,2],pch=19,cex=2,
                 xlab="", ylab="Recruitment")}
  else{plot(apply(totalRain[1:4,2:18,i],2,sum),M$mean$gam1[1:17,2],pch=19,cex=2,
            ,xlab="", ylab="")}
  
  abline(lm(M$mean$gam1[,1]~apply(totalRain[1:4,2:18,i],2,sum)))
}
length(apply(totalRain[1:4,2:18,4],2,sum))
traceplot(M,"gam1")
