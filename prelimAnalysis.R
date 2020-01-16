#PRELIMINARY DATA DISPLAY AND ANALYSIS - MAINLY FOR MESSING AROUND WITH THE DATA AND SEEING WHAT WON'T WORK

#Load data and libraries
load('./data/cleanData.Rdata')
library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
library(sf)

# Basic maps of sites/trap distribution -----------------------------------

#Calgary-centred projection from Paul - seems to be problems when using this
abProjString <- '+proj=tmerc +lat_0=0 +lon_0=-114 +k=0.9999 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'

#Provincial map
ab <- st_read('~/Documents/shapefiles/provBorders/provBorders.shp') %>% #Read in provincial borders
  filter(PRNAME=='Alberta') %>% #Alberta only
  st_transform(3403) #Transform to NAD83, UTM 11N (EPSG:2153), or try EPSG:3403
  # st_transform(abProjString) #Custom projection from Paul
# st_write(ab,'~/Documents/shapefiles/provBorders/AB_only.shp',driver='ESRI Shapefile',delete_layer=T)

#Cities
cit <- st_read('~/Documents/shapefiles/ABlocations/alberta_location.shp') %>% 
  filter(PLACE=='city'|NAME=='Vermilion') %>% 
  filter(!is.na(NAME),NAME!='Airdrie',NAME!='Delia',NAME!='Cold Lake') %>% 
  st_transform(3403) 
# st_write(cit,'~/Documents/shapefiles/ABcities',driver='ESRI Shapefile',delete_layer=T)

#City labels of interest
citLabs <- cit %>% 
  filter(NAME %in% c('Grande Prairie','Edmonton','Red Deer','Calgary','Brooks','Lethbridge','Vermilion')) %>%
  # mutate(NAME=gsub(' ','\n',NAME)) %>% 
  select(NAME,geometry) %>% 
  mutate(xnudge=c(7e4,0,-5e4,0,-5e4,4e4,-2e4),ynudge=c(0,2e4,0,2e4,0,0,2e4)) #Distances to nudge city labels

#Add coordinate system to site df
siteSf <- st_as_sf(site,coords=c('lon','lat'),crs=4326) %>% 
  st_transform(3403)

#Add coordinate system to trap df
trapSf <- trap %>% left_join(select(site,BLID:lon),by='BLID') %>% 
  rename('latSite'='lat','lonSite'='lon') %>% 
  #Replace trap lat/lon with site lat/lon if not recorded
  mutate(latTrap=ifelse(latTrap==0,latSite,latTrap),lonTrap=ifelse(lonTrap==0,lonSite,lonTrap)) %>% 
  st_as_sf(coords=c('lonTrap','latTrap'),crs=4326) %>% 
  st_transform(3403)

#Bounding box for 2016:2018 sites (excludes Grande Prairie)
bounds <- st_bbox(filter(trapSf,startYear!=2015))
  
#Overall Site distributions
(p1 <- ggplot(siteSf)+ geom_sf(size=0.3)+
  geom_sf_text(data=citLabs,aes(label=NAME),col='red',size=3,
               nudge_y=citLabs$ynudge,nudge_x=citLabs$xnudge)+
  geom_sf(data=citLabs,col='red',size=1)+
  geom_sf(data=ab,fill=NA) + 
  coord_sf(xlim=st_bbox(siteSf)[c(1,3)],ylim=st_bbox(siteSf)[c(2,4)])+
  ggtitle('Overall site distributions')+labs(x=NULL,y=NULL))
# ggsave('./figures/maps/overallMap.png',p1)

#Trapping type distribution
(p1 <- trapSf %>% mutate(startYear=factor(startYear)) %>% 
  ggplot()+ geom_sf(size=0.3)+
  # geom_sf_text(data=citLabs,aes(label=NAME),col='red',size=3,
               # nudge_y=citLabs$ynudge,nudge_x=citLabs$xnudge)+
  facet_grid(trapType~startYear)+
  geom_sf(data=ab,fill=NA)+
  geom_sf(data=cit,col='red',size=1)+
  coord_sf(xlim=bounds[c(1,3)],ylim=bounds[c(2,4)])+
  ggtitle('Trap type distributions')+
  theme(axis.text=element_text(size=5),axis.text.x=element_text(angle=90,vjust=0.3,hjust=1)))
# ggsave('./figures/maps/trapTypeMap.png',p1)

#Non-ditch trapping (none from 2015)
nonDitch <- trapSf %>% mutate(startYear=factor(startYear)) %>% 
  filter(trapLoc!='ditch')

(p1 <- ggplot(nonDitch)+ geom_sf(size=0.3)+
  facet_grid(trapType~startYear)+
  geom_sf(data=ab,fill=NA)+
  geom_sf(data=cit,col='red',size=1)+
  # geom_sf_text(data=citLabs,aes(label=NAME),col='red',size=3,nudge_y=citLabs$ynudge,nudge_x=citLabs$xnudge)+
  coord_sf(xlim=bounds[c(1,3)],ylim=bounds[c(2,4)]) + 
  ggtitle('Wetland & in-field trap distributions')+
  theme(axis.text=element_text(size=5),axis.text.x=element_text(angle=90,vjust=0.3,hjust=1)))
# ggsave('./figures/maps/nonDitchMap.png',p1)

#Maps of arthropod groups (arthOrder) used to check completeness
arthNumbers <- arth %>% filter(!mismatchBTID,!mismatchBLID) %>% 
  group_by(BLID,year,arthOrder) %>% summarize(n=n()) %>% 
  left_join(select(siteSf,BLID),by='BLID') 
yearNumbers <- arthNumbers %>% 
  group_by(BLID,year) %>% summarize(n=sum(n)) %>% 
  left_join(select(siteSf,BLID),by='BLID') 
  

#Arth maps
(p1 <- ggplot()+ geom_sf(data=yearNumbers,aes(geometry=geometry),size=0.3,shape=1,alpha=0.3)+
  geom_sf(data=arthNumbers,aes(geometry=geometry),size=0.3)+
  geom_sf(data=cit,col='red',size=1)+ #Cities
  facet_grid(arthOrder~year)+
  coord_sf(xlim=bounds[c(1,3)],ylim=bounds[c(2,4)])+
  ggtitle('Sampling and occurance distribution')+
  theme(axis.text=element_text(size=5),axis.text.x=element_text(angle=90,vjust=0.3,hjust=1)))
  # ggsave('./figures/maps/arthDistMap.png',p1,scale=1.5)

rm(p1)
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
#   select(BLID,replicate,trapType) %>% unique() %>% 
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

#Get wetland distances from infield sites -------------------

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

# Question 1: Spillover into crops ----------------------------------------

# How does insect abundance change with distance into cropland?






