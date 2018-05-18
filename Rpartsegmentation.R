#install.packages("itcSegment")
library(itcSegment)
#install.packages("rgdal")
library(rgdal)
library(raster)
#install.packages("rgeos")
#install.packages("sp")
#install.packages("rgdal")
library (rgeos)
library(sp)
library(rgdal)
library(dplyr)




####################################
####################################
####################################
####################################


## load layers
d1 <- raster("C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/Segmentation/pargreen2006.tif")
d2 <- raster("C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/Segmentation/pargreen2007.tif")
d3 <- raster("C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/Segmentation/pargreen2008.tif")
d4 <- raster("C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/Segmentation/pargreen2015.tif")
d5 <- raster("C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/Segmentation/pargreen2016.tif")
d6 <- raster("C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/Segmentation/pargreen2017.tif")
d7 <- raster("C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/Segmentation/pargreen2006-2007-2008.tif")

meanIgnoringZeroes=function(x){mean(x[x>0.6],na.rm=T)}
s<-overlay(d1, d2, d3,fun=meanIgnoringZeroes)
#plot(s)
s2<-overlay(d4, d5, d6,fun=meanIgnoringZeroes)
#plot(s2)

change<-(s2-s)
plot(change)
 
writeRaster(change, file="C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/change4.tif", format = "ascii")
### uncomment to write the raster file
#setwd("C:\\Users\\cristinabarberal\\Documents\\PhD\\Parcels data for selection\\")
#writeRaster(dif, "GreenDIf17ASDF_06.tif", format = "ascii")



######################################
#####################################


#After segmentation we will resterize the segmented layer, stack it with the change layer and then aggregate 
#so we have every pixel with the segmented poligon ID
seg<-readOGR("C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/Segmentation/layer.shp")
#change to UTM
#spTransform(meuse, CRS("+proj=longlat +datum=WGS84"))
#new("CRS"
 #   , projargs = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#)
plot(seg)

seg_xy=gCentroid(seg,byid=T)
points(seg_xy,pch=19,col="red")



ras<-rasterize(seg, change)
plot(ras)
st<-stack(ras,change)
plot(st)
az<-rasterToPoints(st)
azdata<-data.frame(az)


#Eliminate the poligons out of the study area
azdata2<-filter(azdata, azdata$layer.1!=233)
azdata2<-filter(azdata2, azdata2$layer.1!=254)
azdata2<-filter(azdata2, azdata2$layer.1!=256)
str(azdata2)
head(azdata2)
plot(azdata2$layer.1~azdata2$layer.2)
summary(azdata2$layer.1)


#take the max to create an empty data frame to put the new values
cmean<-matrix(data=NA, nrow=252, ncol = 2)
colnames(cmean)<-c("polygon","meanchange")
v<-sort(unique(azdata2$layer.1))
cmean[,1]<-v
head(cmean)
claschange<-data.frame(cmean)
for (i in 1:length(cmean[,1])) {
  claschange$meanchange[i]<-mean(azdata2$layer.2[azdata2$layer.1==i], na.rm=TRUE)
}

plot(claschange$meanchange~claschange$polygon)
claschange$category<-NA
claschange$category[which(claschange$meanchange>-0.09 & claschange$meanchange<=0.01)]=1
claschange$category[which(claschange$meanchange>0.01 & claschange$meanchange<=0.03)]=2
claschange$category[which(claschange$meanchange>0.03 & claschange$meanchange<=0.06)]=3
claschange$category[which(claschange$meanchange>0.06 & claschange$meanchange<=0.08)]=4
claschange$category[which(claschange$meanchange>0.08 & claschange$meanchange<=0.10)]=5
claschange$category[which(claschange$meanchange>0.10 & claschange$meanchange<=0.12)]=6
claschange$category[which(claschange$meanchange>0.12 & claschange$meanchange<=0.14)]=7
claschange$category[which(claschange$meanchange>0.14 & claschange$meanchange<=0.16)]=8
claschange$category[which(claschange$meanchange>0.16 & claschange$meanchange<=0.18)]=9
claschange$category[which(claschange$meanchange>0.18 & claschange$meanchange<=0.19)]=10
for (i in 1:length(azdata2$layer.1)) {
  azdata2$category[i]=claschange$category[which(claschange$polygon==azdata2$layer.1[i])]
}
write.csv(azdata2,file ="C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/NDVIchange5.csv" )

#rasterize our data base
azdata3<-azdata2
azdata3$layer.1<-NULL
azdata3$layer.2<-NULL
se<-azdata3$category
xy<-cbind(azdata3$x, azdata3$y)
xyz<-cbind(xy,se)
ra<-rasterFromXYZ(xyz)
plot(ra)
writeRaster(ra,file='C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/cleanclass.tif')

###############################################################
###############################################################
#Radnomly select our transects
#category poligons
catpols<-rasterToPolygons(ra, dissolve=TRUE)
par(mfrow=c(1,1))
plot(catpols)
plot(ra)

#random sample
point<-sample(seg_xy, 28)
plot(point, add=TRUE, pch=1, col='red')
#write.csv(point,file='C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/randocent.csv')

#extract coordinates
transects<-data.frame(point@coords)
transects<-rbind(transects,c(-80.17088, 7.432097), c(-80.17657, 7.427357))
#emty matrix to allocate the coordinates of the lines going from N/S and E/W pasing by the centroids
transects1<-matrix(NA, ncol = 2, nrow = 120)
tra<-matrix(NA,ncol = 2, nrow = 2)
#loop to create the two points to conect to obtain the lines
for( i in 1:30) {
  transects1[i,1]<-transects$x[i]
  transects1[i,2]<-7.54
  transects1[i+30,1]<-transects$x[i]
  transects1[i+30,2]<-7.42
  transects1[i+60,1]<--80.26
  transects1[i+60,2]<-transects$y[i]
  transects1[i+90,1]<--80.14
  transects1[i+90,2]<-transects$y[i]
}

#loop to create a matrix with the two extreme coordinates to create a Line class list
#Vertical lines
Tall<-vector('list', 60)
for (i in 1:30) {
  tra[1,1]=transects1[i,1]
  tra[1,2]=transects1[i,2]
  tra[2,1]=transects1[i+30,1]
  tra[2,2]=transects1[i+30,2]
  tra2<-SpatialLines(list(Lines(list(Line(tra)), "id")))
  Tall[i]<-tra2
}
#Horizontal lines
for (i in 61:90) {
  tra[1,1]=transects1[i,1]
  tra[1,2]=transects1[i,2]
  tra[2,1]=transects1[i+30,1]
  tra[2,2]=transects1[i+30,2]
  tra2<-SpatialLines(list(Lines(list(Line(tra)), "id")))
  Tall[i-30]<-tra2
}
#Check the result
plot(seg)
plot(point,col='red',add=TRUE)
for (i in 1:90) {
  plot(Tall[[i]], add=TRUE)
}

#save the result to export to Qgis
alllines=do.call("rbind",Tall)

crs(alllines)=crs(seg)
lines_utm=spTransform(x=croplines,CRS="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

df <- data.frame(len = sapply(1:length(lines_utm), function(i) gLength(lines_utm[i, ])))
rownames(df) <- sapply(1:length(lines_utm), function(i) lines_utm@lines[[i]]@ID)

Sldf <- SpatialLinesDataFrame(lines_utm, data = df)




segUTM=spTransform(x=seg,CRS="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")




writeOGR(Sldf,layer="Sldf",driver="ESRI Shapefile",dsn=".",overwrite=T)

all_intersections=gIntersection(Sldf,segUTM,byid=T)
pldf=gBuffer(Sldf,10,byid=T)
writeOGR(pldf,layer="pldf",driver="ESRI Shapefile",dsn=".",overwrite=T)



go=as.data.frame(alllines)

setwd('C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection')
save(Tall,file='C:/Users/cristinabarberal/Documents/PhD/Parcels data for selection/transects.tif')
