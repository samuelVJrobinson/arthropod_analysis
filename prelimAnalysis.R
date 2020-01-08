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
  filter(!is.na(NAME),NAME!='Airdrie',NAME!='Delia',NAME!='Cold Lake',NAME!='Brooks') %>% 
  st_transform(3403) 
# st_write(cit,'~/Documents/shapefiles/ABcities',driver='ESRI Shapefile')

#Add coordinate system to site df
siteSf <- st_as_sf(site,coords=c('lon','lat'),crs=4326) %>% 
  st_transform(3403)

#Overall Site distributions
ggplot(siteSf)+
  geom_sf(size=1)+
  # geom_sf_text(data=cit,aes(label=NAME),col='red',size=3,nudge_y=20000)+
  geom_sf(data=cit,col='red',size=1)+
  geom_sf(data=ab,fill=NA) + 
  coord_sf(xlim=st_bbox(siteSf)[c(1,3)],ylim=st_bbox(siteSf)[c(2,4)])

#Trapping type distribution
trap %>% group_by(BLID,startYear,trapType) %>% 
  summarize() %>% ungroup() %>% 
  left_join(select(siteSf,BLID),by='BLID') %>% 
  mutate(startYear=factor(startYear)) %>% 
  ggplot(aes(geometry=geometry))+
  geom_sf(size=1,shape=4)+
  facet_grid(trapType~startYear)+
  geom_sf(data=ab,fill=NA)+
  geom_sf(data=cit,col='red',size=1)+
  coord_sf(xlim=st_bbox(siteSf)[c(1,3)],ylim=st_bbox(siteSf)[c(2,4)])

#Non-ditch trapping (none from 2015)
trap %>% filter(trapLoc!='ditch') %>% 
  group_by(BLID,startYear,trapType) %>%
  summarize() %>% ungroup() %>% 
  left_join(select(siteSf,BLID),by='BLID') %>% 
  mutate(startYear=factor(startYear)) %>% 
  ggplot(aes(geometry=geometry))+
  geom_sf(shape=4)+
  facet_grid(trapType~startYear)+
  geom_sf(data=ab,fill=NA)+
  geom_sf(data=cit,col='red',size=1)+
  coord_sf(xlim=st_bbox(siteSf)[c(1,3)],ylim=st_bbox(siteSf)[c(2,4)])
  

