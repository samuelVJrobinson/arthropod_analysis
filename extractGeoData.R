#SCRIPT TO EXTRACT GEOGRAPHIC INFORMATION SURROUNDING TRAPPING SITES
#SR WINTER 2020

library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
library(sf)

# Helper functions --------------------------------------------------------

#Function to extract cover within onion rings surrouding site
# centLoc: central point around which to get onion rings
# coverShp: set of polygons representing some cover class
# ringDists: vector of distances to the outside of each ring (inside is radius of last ring). Zero implied at beginning (first ring is a circle).
getOnionRings <- function(centLoc,coverShp,ringDists){
  #Checks
  if(length(ringDists)<1) stop('ringDists must specify 1 or more distances')
  if(any(ringDists<=0)) stop('ringDists must be greater than 0')
  if(nrow(centLoc)!=1) stop('centLoc must have 1 point only')
  #Create list of buffer circles 
  bObjs <- lapply(ringDists,function(x) st_buffer(centLoc,dist=x))
  #Get intersection area between coverShp and consecutive circles
  overlapAreas <- sapply(bObjs,function(x) st_intersection(st_geometry(x),st_geometry(coverShp)) %>% st_area() %>% sum())
  #Subtract area from smaller circle to get area of intersection only within that ring
  overlapAreas[2:length(overlapAreas)] <- overlapAreas[2:length(overlapAreas)]-overlapAreas[1:length(overlapAreas)-1]
  #Get area from each buffer circle
  ringAreas <- sapply(bObjs,st_area)
  #Subtract area from smaller circle to get area within that ring
  ringAreas[2:length(ringAreas)] <- ringAreas[2:length(ringAreas)]-ringAreas[1:length(ringAreas)-1]
  #Assemble into data frame
  a <- data.frame(ringDists=ringDists,overlap=overlapAreas,ringArea=ringAreas)
  return(a)
}
# #Function works
# b1 <- site[2,] #Test centre
# b3 <- seq(20,1000,20) #Test ring dists
# b2 <- st_geometry(st_buffer(b1,dist=max(b3))) %>% #Buffer maximum distance around centLoc
#   st_intersection(st_geometry(allWetlands)) #Get intersection between buffer and coverShp
# 
# getOnionRings(centLoc=b1,coverShp=b2,ringDists=b3) %>% 
#   mutate(prop=overlap/ringArea) %>% 
#   ggplot(aes(ringDists,prop))+geom_point()+geom_line()

# Load shapefiles ---------------------------------------------------------

load('./data/cleanData.Rdata') #Load site, trap, arth data

#Set coordinates
site <- st_as_sf(site,coords=c('lon','lat'),crs=4326) %>% st_transform(3403)
trap <- st_as_sf(trap,coords=c('lonTrap','latTrap'),crs=4326) %>% st_transform(3403)

#Load 
allEphemeral <- st_read('~/Documents/shapefiles/mergedFeatures/allEphemeral.shp',crs=3403) %>% st_geometry()
allGrass <- st_read('~/Documents/shapefiles/mergedFeatures/allGrass.shp',crs=3403) %>% st_geometry()
allShrubs <- st_read('~/Documents/shapefiles/mergedFeatures/allShrubs.shp',crs=3403) %>% st_geometry()
allTrees <- st_read('~/Documents/shapefiles/mergedFeatures/allTrees.shp',crs=3403) %>% st_geometry()
allWetlands <- st_read('~/Documents/shapefiles/mergedFeatures/allWetlands.shp',crs=3403) %>% st_geometry()
allNonCrop <- st_read('~/Documents/shapefiles/mergedFeatures/allNonCrop.shp',crs=3403) %>% st_geometry()


# Get onion-ring distances from trapping points from 2016-2017 ---------------------------

#Create "alternate" BTID in trap df
trap <- trap %>% 
  mutate(altBTID=paste(BLID,dist,startYear,sep='-'))
  # mutate(altBTID=sub('-\\d-','-X-',BTID)) #Doesn't work - contains trap type within replicate

#Get unique trap locations from 2016-2017
trapU <- trap %>% 
  filter(startYear==2016|startYear==2017) %>% 
  select(altBTID,geometry) %>% 
  unique.data.frame() %>% 
  slice(1:10)
  
# ggplot()+
#   geom_sf(data=allNonCrop,fill='green',col='green')+
#   geom_sf(data=trapU,size=0.2)+
#   coord_sf(xlim=st_bbox(trapU)[c(1,3)],ylim=st_bbox(trapU)[c(2,4)])

#Trap locations in list form
trapList <- lapply(1:nrow(trapU),function(x) trapU[x,])
rDists <- seq(20,1000,20)
oRingNonCrop <- lapply(trapList,function(x){ #Extract area of cover at all distances
  b1 <- st_geometry(st_buffer(x,dist=max(rDists))) %>% #Buffer maximum distance around centLoc
    st_intersection(st_geometry(allNonCrop)) #Get intersection between buffer and coverShp
  return(getOnionRings(centLoc=x,coverShp=b1,ringDists=rDists))
})

#Convert to matrix form
oRingNonCropMat <- t(matrix(sapply(oRingNonCrop,function(x) x$overlap),ncol=nrow(trapU),
                          dimnames=list(paste0('noncrop',rDists),paste0('s',1:10))))



# Create new shapefiles for 2018 data --------------------------------
# 
# #Add ditch sites from 2015/2016/2018
# for(i in c(2015,2016,2018)){
#   trapSf %>% filter(startYear==i,trapLoc=='ditch') %>% 
#     select(BLID,replicate,trapType) %>% unique() %>% 
#     st_transform(3403) %>% 
#     st_write(paste0('~/Documents/shapefiles/digitizedFeaturesSR/DitchSites_',i,'.shp'),delete_dsn=TRUE)
# }
# 
# #Add infields from 2018
# trapSf %>% filter(startYear==2018,trapLoc!='ditch') %>%
#   select(BLID,BTID,trapType,trapLoc,dist,distFrom) %>% unique() %>%
#   st_transform(3403) %>%
#   st_write('~/Documents/shapefiles/digitizedFeaturesSR/InfieldSites_2018.shp',delete_dsn=TRUE)
# 
# #Get various shapefiles from GDB (can't be written to outside of ArcMap), transform, and save in another folder
# layers <- st_layers('~/Documents/shapefiles/DigitizedFeatures_AllYears.gdb')$name
# layers <- layers[!grepl('Buffer',layers)] #Get rid of buffer layers
#   
# for(i in layers){
#   j <- st_read(dsn='~/Documents/shapefiles/DigitizedFeatures_AllYears.gdb',layer=i) %>% 
#     st_zm() %>% st_transform(3403)
#   st_write(j,paste0('~/Documents/shapefiles/digitizedFeaturesSR/',i,'.shp'),driver='ESRI Shapefile',delete_layer=TRUE)
# }
# 
# #Only 16 of the ditch sites from 2018 (repeats of previous years in S. Calgary) have been digitized.
# #The remaining ditch sites (41) and infields (27) need wetland layers created around them

#Get nearest wetland distances from infield sites -------------------

#Wetland distances 
wtrap2016 <- trapSf %>% filter(startYear==2016,grepl('W',replicate)) #trap locations
wland2016 <- st_read(dsn='~/Documents/shapefiles/digitizedFeaturesSR/Wetlands2016.shp') %>%  #wetland polygons
  rbind(select(st_read(dsn='~/Documents/shapefiles/digitizedFeaturesSR/Ephemeral2016.shp'),-Notes)) %>%  # include ephemeral wetlands
  mutate(wType=ifelse(grepl('Wetland',FldrPth),'Permanant','Ephemeral'))

#Buffer all unique traps by 200m
trapBuff <- wtrap2016 %>% select(BLID) %>% unique() %>% st_buffer(dist=300)

#Filter wetlands that overlap buffer distance
wland2016 <- wland2016 %>% slice(sort(unique(unlist(st_intersects(trapBuff,wland2016))))) 
rm(trapBuff)

#Identify wetlands nearest to traps
wland2016$isNearest <- 1:nrow(wland2016) %in% st_nearest_feature(wtrap2016,wland2016)

wland2016Points <- st_cast(wland2016,to='POINT') #Cast wetland to point geometry
wland2016Points$isNearest <- 1:nrow(wland2016Points) %in% st_nearest_feature(wtrap2016,wland2016Points)  #Select points nearest to traps

wtrap2016 <- wtrap2016 %>% 
  #Gets distance from nearest wetland, unless replicate is at W0 (dist=0)
  mutate(wlandDist = ifelse(grepl('W0',replicate),0,apply(st_distance(wtrap2016,filter(wland2016Points,isNearest)),1,min)))

ggplot(wtrap2016) + geom_histogram(aes(x=wlandDist)) #Looks good


#Appears to work
coords <- st_coordinates(distinct(wtrap2016))[1,]
#Polygons
ggplot()+
  geom_sf(data=wland2016,aes(fill=wType))+
  geom_sf(data=wtrap2016,col='blue')+
  coord_sf(xlim=c(1,-1)*500+coords[1],ylim=c(1,-1)*500+coords[2])+
  scale_fill_manual(values=c('lightgreen','forestgreen'))

#Points
ggplot()+
  geom_sf(data=wland2016Points,aes(col=isNearest),size=0.5)+
  geom_sf(data=wtrap2016,col='blue',size=0.5)+
  coord_sf(xlim=c(1,-1)*400+coords[1],ylim=c(1,-1)*400+coords[2])+
  scale_colour_manual(values=c('black','red'))
