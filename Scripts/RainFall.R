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

for(i in 1:12){
    for(p in 1:4){
      
rain[i,1,p]<-extract(mean(stack(months[[i]])),place[[p]],fun=mean)

for(y in 1:18){
  rain[i,y+1,p]<-extract(months[[i]][[y]],place[[p]],fun=mean)-rain[i,1,p]
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



# Read in forest loss per year #
rasterOptions(tmpdir="D:/TempRaster")

hi<-raster("lossYear.tif")
HispForest<-crop(hi,Hisp)
JamForest<-crop(hi,Jamaica)
PRForest<-crop(hi,PR)
CubaForst<-crop(hi,Cuba)

reclassArray<-array(0,c(14,2,14))
reclassArray[,1,]<-1:14
for(i in 1:14){
reclassArray[i,2,i]<-1
}



LossYearsHisp<-LossYearsCuba<-LossYearsPR<-LossYearsJam<-rep(NA,14)
a<-Sys.time()
Sys.time()
for(i in 1:14){
LossYearsHisp[i]<-cellStats((area(HispForest)*reclassify(HispForest,reclassArray[,,i])),sum)
LossYearsCuba[i]<-cellStats((area(CubaForst)*reclassify(CubaForst,reclassArray[,,i])),sum)
LossYearsPR[i]<-cellStats((area(PRForest)*reclassify(PRForest,reclassArray[,,i])),sum)
LossYearsJam[i]<-cellStats((area(JamForest)*reclassify(JamForest,reclassArray[,,i])),sum)
}
Sys.time()-a


par(mfrow=c(2,1))
par(bty="l")
plot(1:18,PRrain[1,2:19,1],pch=19,type="o",col="blue",ylim=c(-150,150),yaxt="n",xaxt="n",ylab="Normalized Rain (mm)") # PR
par(new=TRUE)
plot(1:18,PRrain[1,2:19,2],pch=19,type="o",col="red",ylim=c(-150,150),yaxt="n",xaxt="n",ylab="",xlab="") # Cuba
par(new=TRUE)
plot(1:18,PRrain[1,2:19,3],pch=19,type="o",col="black",ylim=c(-150,150),yaxt="n",xaxt="n",ylab="",xlab="") # Hisp
par(new=TRUE)
plot(1:18,PRrain[1,2:19,4],pch=19,type="o",col="green",ylim=c(-150,150),yaxt="n",xaxt="n",ylab="",xlab="") # Jamaica
abline(h=0,col="gray44",lty=2)
legend(1,-100,legend=c("PR","Cuba","Hisp","Jamaica"),pch=rep(19,4),lty=rep(1,4),col=c("blue","red","black","green"),bty="n")
axis(2,las=2)
axis(1,label=Year,at=1:18)
par(bty="l")
plot(1:14,LossYearsPR,pch=19,type="o",col="blue",ylim=c(0,300),yaxt="n",xaxt="n",ylab="Forest Lost km2") # PR
par(new=TRUE)
plot(1:14,LossYearsCuba,pch=19,type="o",col="red",ylim=c(0,300),yaxt="n",xaxt="n",ylab="",xlab="") # Cuba
par(new=TRUE)
plot(1:14,LossYearsJam,pch=19,type="o",col="green",ylim=c(0,300),yaxt="n",xaxt="n",ylab="",xlab="") # Cuba
par(new=TRUE)
plot(1:14,LossYearsHisp,pch=19,type="o",col="black",ylim=c(0,300),yaxt="n",xaxt="n",ylab="",xlab="") # Cuba
axis(1,label=c(2001:2014),at=1:14)
axis(2,las=2)

