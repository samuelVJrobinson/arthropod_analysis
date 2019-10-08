# CREATES SUMMARY TABLES/MAPS FROM CLEANED DATA

load('cleanData.Rdata') #Load data

library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
theme_set(theme_bw())

# Map data for simple plots ------------------------------------------

#Goal: Distribution of trapping effort over years, split by trap type

#Get trapping effort for all site combos
trap2 <- trap %>% rowwise() %>% 
  transmute(BLID,type=trapType,year=max(c(startYear,endYear),na.rm=T),hours=deployedhours) %>% 
  ungroup() %>% mutate(type=factor(type),year=factor(year)) 
  
  
  group_by(BLID,type,year) %>% summarize(hours=sum(hours,na.rm=T))



#Types of map projections
targetCRS <- 32612 #UTM 12
# targetCRS <- 3771 #Alberta Albers equal area
# targetCRS <- 102009 #Lambert conformal conic (tilted too much)

#Shapefile for outline of AB
ABoutline <- st_read("C:\\Users\\Samuel\\Documents\\Shapefiles\\Canada provincial\\alberta.shp")
ABoutline <- ABoutline %>% st_transform(targetCRS) #Transforms projection

#City names in AB
ABplaces <- st_read("C:\\Users\\Samuel\\Documents\\Shapefiles\\Alberta place names\\cgn_ab_shp_eng.shp") %>% 
  mutate(REL_SCALE=as.numeric(as.character(REL_SCALE))) %>% 
  filter(CONCISE=='CITY') %>% 
  top_n(n=6,wt=REL_SCALE) #Top 6 cities
ABplaces <- ABplaces %>% st_transform(targetCRS) #Transforms projection

#Adds coordinate system to site df
site <- st_as_sf(site,agr='identity',coords=c('lon','lat'))
st_crs(site) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" #Initial projection (latlon) 
site <- site %>% st_transform(targetCRS) #Transforms projection

# #Plot using plot.sf
# plot(st_geometry(ABoutline))
# plot(st_geometry(site),add=T,pch=19,cex=0.5)
# plot(st_geometry(ABplaces),add=T,pch=19,cex=1,col='red')

#Plot using ggplot
ggplot()+
  geom_sf(data=ABoutline,fill='white') + #Province outline
  geom_sf(data=site,alpha=0.7)+ #Sampled sites
  geom_sf(data=ABplaces,col='red',size=2) + #City locations
  geom_sf_text(data=ABplaces,aes(label=GEONAME),col='red',size=2,nudge_y=2e4)+ #City Labels
  coord_sf(ylim=range((st_coordinates(site))[,2]))+ #Limit to latitudes of sampling
  labs(x=NULL,y=NULL)



coordinates(trapCoords)=~lon+lat 
proj <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") #Initial projection (latlon) to set 
proj4string(trapCoords) <- proj #Set projection
newproj <- CRS("+init=epsg:32612") #CRS for UTM 12
trapCoords <- spTransform(trapCoords,newproj) #Transform to UTM 12
# writeOGR(trapCoords,"C:\\Users\\Samuel\\Documents\\Shapefiles\\Field map\\BlueVaneTrap Locations","BVTraps_2015_2016_Continuous",driver="ESRI Shapefile") #Writes trap coordinates to shapefile directory
distMat <- gDistance(trapCoords,byid=T) #Get distance matrix between sites
colnames(distMat) <- rownames(distMat) <- trapCoords$BLID #Rename rows and columns
distMat <- distMat/10000 #Scales distances to 10s of kms
