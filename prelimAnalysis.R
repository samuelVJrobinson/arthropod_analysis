#PRELIMINARY DATA DISPLAY AND ANALYSIS - MAINLY FOR MESSING AROUND WITH THE DATA AND SEEING WHAT WON'T WORK

#Load data and libraries
load('./data/cleanData.Rdata')
library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA)) #Better for maps
library(sf)

# Basic maps of sites/trap distribution -----------------------------------

#Provincial map
ab <- st_read('~/Documents/shapefiles/provBorders') %>% #Read in provincial borders
  filter(PRNAME=='Alberta') %>% #Alberta only
  st_transform(3403) #Transform to NAD83, UTM 11N (EPSG:2153), or try EPSG:3403  

#Cities
cit <- st_read('~/Documents/shapefiles/ABcities') %>% filter(PLACE=='city') %>% 
  filter(!is.na(NAME),NAME!='Airdrie',NAME!='Delia',NAME!='Cold Lake') %>% 
  st_transform(3403) 
# st_write(cit,'~/Documents/shapefiles/ABcities',driver='ESRI Shapefile')

#City labels of interest
citLabs <- cit %>% 
  filter(NAME %in% c('Grande Prairie','Edmonton','Red Deer','Calgary','Brooks','Lethbridge')) %>%
  # mutate(NAME=gsub(' ','\n',NAME)) %>% 
  select(NAME,geometry) %>% 
  mutate(xnudge=c(8e4,0,-6e4,0,-6e4),ynudge=c(0,2e4,0,2e4,0))

#Distances to nudge city labels

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
  
#Overall Site distributions
ggplot(siteSf)+
  geom_sf(size=1)+
  geom_sf_text(data=citLabs,aes(label=NAME),col='red',size=3,
               nudge_y=citLabs$ynudge,nudge_x=citLabs$xnudge)+
  geom_sf(data=citLabs,col='red',size=1)+
  geom_sf(data=ab,fill=NA) + 
  coord_sf(xlim=st_bbox(siteSf)[c(1,3)],ylim=st_bbox(siteSf)[c(2,4)])+
  ggtitle('Overall site distributions')

#Trapping type distribution
trapSf %>% mutate(startYear=factor(startYear)) %>% 
  ggplot()+
  geom_sf(size=1,shape=4)+
  # geom_sf_text(data=citLabs,aes(label=NAME),col='red',size=3,
               # nudge_y=citLabs$ynudge,nudge_x=citLabs$xnudge)+
  facet_grid(trapType~startYear)+
  geom_sf(data=ab,fill=NA)+
  geom_sf(data=cit,col='red',size=1)+
  coord_sf(xlim=st_bbox(siteSf)[c(1,3)],ylim=st_bbox(siteSf)[c(2,4)])

#Non-ditch trapping (none from 2015)
nonDitch <- trapSf %>% mutate(startYear=factor(startYear)) %>% 
  filter(trapLoc!='ditch')

ggplot(nonDitch)+ geom_sf(shape=4)+
  facet_grid(trapType~startYear)+
  geom_sf(data=ab,fill=NA)+
  geom_sf(data=cit,col='red',size=1)+
  # geom_sf_text(data=citLabs,aes(label=NAME),col='red',size=3,nudge_y=citLabs$ynudge,nudge_x=citLabs$xnudge)+
  coord_sf(xlim=st_bbox(nonDitch)[c(1,3)],ylim=st_bbox(nonDitch)[c(2,4)])

# Question 1: Spillover into crops ----------------------------------------

# How does insect abundance change with distance into cropland?






