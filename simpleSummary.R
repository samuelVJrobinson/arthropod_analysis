# CREATES SUMMARY TABLES/MAPS FROM CLEANED DATA

load('cleanData.Rdata') #Load data

library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
theme_set(theme_classic())

# Map data for simple plots ------------------------------------------

#Goal: Distribution of trapping effort over years, split by trap type

#Shapefile for outline of AB
ABoutline <- st_read("C:\\Users\\Samuel\\Documents\\Shapefiles\\Canada provincial\\alberta.shp")
ABoutline <- ABoutline %>% st_transform(32612) #Transforms to UTM 12

#City names in AB
ABplaces <- st_read("C:\\Users\\Samuel\\Documents\\Shapefiles\\Alberta place names\\cgn_ab_shp_eng.shp") %>% 
  mutate(REL_SCALE=as.numeric(as.character(REL_SCALE))) %>% 
  filter(CONCISE=='CITY') %>% 
  top_n(n=6,wt=REL_SCALE)
ABplaces <- ABplaces %>% st_transform(32612) #Transforms to UTM 12

#Adds coordinate system to site df
siteCoords <- st_as_sf(site,agr='identity',coords=c('lon','lat'))
st_crs(siteCoords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" #Initial projection (latlon) 
siteCoords <- siteCoords %>% st_transform(32612) #Transforms to UTM 12


#Plot using plot.sf
plot(st_geometry(ABoutline))
plot(st_geometry(siteCoords),add=T,pch=19,cex=0.5)
plot(st_geometry(ABplaces),add=T,pch=19,cex=1,col='red')

#Plot using ggplot
ggplot()+
  geom_sf(data=ABoutline,fill='white') +
  geom_sf(data=ABplaces,col='red',size=2) +
  # geom_sf_label(data=ABplaces,aes(label=GEONAME),col='red',size=2)+ #Labels
  geom_sf(data=siteCoords)
  


  
  # trap %>% select(BLID,trapType,deployedhours,startYear) %>% 
  # group_by(BLID) %>%  #Get lat and lon
  # summarize(lon=first(lon),lat=first(lat)) %>% 
  # data.frame()




coordinates(trapCoords)=~lon+lat 
proj <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") #Initial projection (latlon) to set 
proj4string(trapCoords) <- proj #Set projection
newproj <- CRS("+init=epsg:32612") #CRS for UTM 12
trapCoords <- spTransform(trapCoords,newproj) #Transform to UTM 12
# writeOGR(trapCoords,"C:\\Users\\Samuel\\Documents\\Shapefiles\\Field map\\BlueVaneTrap Locations","BVTraps_2015_2016_Continuous",driver="ESRI Shapefile") #Writes trap coordinates to shapefile directory
distMat <- gDistance(trapCoords,byid=T) #Get distance matrix between sites
colnames(distMat) <- rownames(distMat) <- trapCoords$BLID #Rename rows and columns
distMat <- distMat/10000 #Scales distances to 10s of kms
