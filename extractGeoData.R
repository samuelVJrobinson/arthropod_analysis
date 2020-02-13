#SCRIPT TO EXTRACT GEOGRAPHIC INFORMATION SURROUNDING TRAPPING SITES
#SR WINTER 2020

library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
library(sf)
library(beepr)
library(parallel)

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

# #Function to get nearest neighbour distance from 
# getNNdist <- function(centLoc,coverShp){
#   # st_buffer(centLoc,dist=bDist) %>% st_geometry() %>% #Buffer trap to bDist
#   #   st_intersection(coverShp) %>% #Intersection with coverShp
#     st_distance(coverShp,centLoc) %>% #Get nearest distance
#     as.numeric()
# }

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

#Merge "permanent" wetlands with ephemeral - takes ~5s
allWetlands <- st_union(allWetlands,allEphemeral)
# rm(allEphemeral)

#Create ID in trap df to pair with geodata surrounding sites
trap <- trap %>% mutate(ID=paste(BLID,dist,startYear,sep='-'))

#Get unique trap locations from 2016-2017
trapU <- trap %>% 
  #group_by(BLID) %>% #Get rid of sites with only 1 trap
  #mutate(nDist=length(unique(dist))) %>% ungroup() %>%
  # filter(nDist>1) %>% 
  select(ID,geometry) %>%
  unique.data.frame() %>% #Problem: 1 trap (31161-0-2018) has 2 different geometry 
  group_by(ID) %>% mutate(ntraps=1:n()) %>% 
  filter(ntraps==1) %>% select(-ntraps)
# trapU %>% filter(ID=='31161-0-2018') 

#Convert trap locations to list form
trapList <- lapply(1:nrow(trapU),function(x) trapU[x,])

#List of cover class shapefiles
coverList <- list(allEphemeral,allGrass,allShrubs,allTrees,allWetlands,allNonCrop)
names(coverList) <- c('ephemeral','grass','shrubs','trees','wetlands','noncrop')

# Get nearest-neighbour distances from various features -------------------

nnDistMat <- lapply(coverList,function(x){ #Get nearest distance
  as.numeric(st_distance(trapU,x))})
#Convert to matrix
nnDistMat <- do.call('cbind',nnDistMat)
rownames(nnDistMat) <- trapU$ID
colnames(nnDistMat) <- names(coverList)

# Get onion-ring distances from trapping points from 2016-2017 ---------------------------

#Ring distances
rDists <- seq(20,1000,20)

#Total area within each ring
oRingArea <- getOnionRings(centLoc=trapList[[1]],
                           coverShp=st_intersection(st_geometry(st_buffer(trapList[[1]],dist=max(rDists))),st_geometry(allNonCrop)),
                           ringDists=rDists)$ringArea

#Get matrices for each cover class - takes about 10 mins if using the entire dataset
f <- function(shp){
  #Extract non-crop cover within rings
  temp <- lapply(trapList,function(x){ 
    b1 <- st_geometry(st_buffer(x,dist=max(rDists))) %>% #Buffer maximum distance around centLoc
      st_intersection(st_geometry(shp)) #Get intersection between buffer and coverShp
    b2 <- getOnionRings(centLoc=x,coverShp=b1,ringDists=rDists)
    return(b2)
  })
  #Convert to matrix form
  tempMat <- t(matrix(sapply(temp,function(x) x$overlap),ncol=nrow(trapU),
                      dimnames=list(paste0('d',rDists),trapU$ID)))
  return(tempMat)
}

# oRingMat <- lapply(coverList,f) #Sequential version
oRingMat <- mclapply(coverList,f,mc.cores=10) #Parallel version (mcapply doesn't work on Windows)
beep(2)

# Check results
par(mfrow=c(3,2))
for(j in 1:length(oRingMat)){ #Looks OK
  plot(0,0,type='n',ylab='% cover',xlab='Distance',ylim=c(0,max(apply(oRingMat[[j]],2,max)/oRingArea)),
       xlim=range(rDists),main=names(oRingMat)[j])
  for(i in 1:nrow(oRingMat[[j]])) lines(rDists,oRingMat[[j]][i,]/oRingArea,col=alpha('black',0.1))
  lines(rDists,apply(oRingMat[[j]],2,mean)/oRingArea,lwd=2,col='red')
  abline(h=0,col='red',lty='dashed'); abline(h=max(apply(oRingMat[[j]],2,max)/oRingArea),col='red',lty='dashed');
}
par(mfrow=c(1,1))

# Save results ------------------------------------------------------------

save(nnDistMat,oRingMat,oRingArea,file='./data/geoData.Rdata') #Saves to geodata file

# Other geoprocessing code snippets -----------------------

# Create new shapefiles for 2018 data
 
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
