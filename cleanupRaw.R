# THIS FILE CLEANS AND SUMMARIZES WETLAND AND INFIELD DATA FROM ECOLOGICS

library(dplyr)
library(tidyr)

#Load raw data
load('rawData.Rdata')


# Clean up site data ------------------------------------------------------
# head(site)

# Clean up trap data ------------------------------------------------------

trap2 <- trap %>% 
  #"Coloured cups" consist of blue, white & yellow cup, so triple effort makes sense
  mutate(trapType=as.character(trapType)) %>% 
  mutate(deployedhours=case_when(trapType=='Coloured Cups' ~ deployedhours*3,
                                 TRUE ~ deployedhours)) 

#TO DO: consider changing all "X Cup" values to "Coloured Cups". Would have to be done in trap (BTID/trapType) and arthropod data (BTID).
  


  
trap2 %>% filter(!grepl('Coloured Cups',trap$trapType))



# Clean up arthropod data -------------------------------------------------
# head(arth)

# Load map data for simple plots ------------------------------------------
library(sf)
# library(sp)
# library(maptools)
# library(rgdal)
# library(rgeos)

#Distribution of trapping effort over years, split by trap type
trapCoords <- trap %>% select(BLID,trapType,deployedhours,startYear) %>% 
  
  group_by(BLID) %>%  #Get lat and lon
  summarize(lon=first(lon),lat=first(lat)) %>% 
  data.frame()




coordinates(trapCoords)=~lon+lat 
proj <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") #Initial projection (latlon) to set
proj4string(trapCoords) <- proj #Set projection
newproj <- CRS("+init=epsg:32612") #CRS for UTM 12
trapCoords <- spTransform(trapCoords,newproj) #Transform to UTM 12
# writeOGR(trapCoords,"C:\\Users\\Samuel\\Documents\\Shapefiles\\Field map\\BlueVaneTrap Locations","BVTraps_2015_2016_Continuous",driver="ESRI Shapefile") #Writes trap coordinates to shapefile directory
distMat <- gDistance(trapCoords,byid=T) #Get distance matrix between sites
colnames(distMat) <- rownames(distMat) <- trapCoords$BLID #Rename rows and columns
distMat <- distMat/10000 #Scales distances to 10s of kms
