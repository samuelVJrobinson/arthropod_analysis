#PRELIMINARY DATA DISPLAY AND ANALYSIS - MAINLY FOR MESSING AROUND WITH THE DATA AND SEEING WHAT WON'T WORK

#Load data and libraries
load('./data/cleanData.Rdata')
library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
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
  mutate(xnudge=c(8e4,0,-6e4,0,-6e4,4e4),ynudge=c(0,2e4,0,2e4,0,-2e4)) #Distances to nudge city labels

#Bounding box for 2016:2018 sites (excludes Grande Prairie)
bounds <- st_bbox(filter(trapSf,startYear!=2015))

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
p1 <- ggplot(siteSf)+ geom_sf(size=0.3)+
  geom_sf_text(data=citLabs,aes(label=NAME),col='red',size=3,
               nudge_y=citLabs$ynudge,nudge_x=citLabs$xnudge)+
  geom_sf(data=citLabs,col='red',size=1)+
  geom_sf(data=ab,fill=NA) + 
  coord_sf(xlim=st_bbox(siteSf)[c(1,3)],ylim=st_bbox(siteSf)[c(2,4)])+
  ggtitle('Overall site distributions')+labs(x=NULL,y=NULL)
ggsave('./figures/maps/overallMap.png',p1)

#Trapping type distribution
p1 <- trapSf %>% mutate(startYear=factor(startYear)) %>% 
  ggplot()+ geom_sf(size=0.3)+
  # geom_sf_text(data=citLabs,aes(label=NAME),col='red',size=3,
               # nudge_y=citLabs$ynudge,nudge_x=citLabs$xnudge)+
  facet_grid(trapType~startYear)+
  geom_sf(data=ab,fill=NA)+
  geom_sf(data=cit,col='red',size=1)+
  coord_sf(xlim=bounds[c(1,3)],ylim=bounds[c(2,4)])+
  ggtitle('Trap type distributions')+
  theme(axis.text=element_text(size=5),axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))
ggsave('./figures/maps/trapTypeMap.png',p1)

#Non-ditch trapping (none from 2015)
nonDitch <- trapSf %>% mutate(startYear=factor(startYear)) %>% 
  filter(trapLoc!='ditch')

p1 <- ggplot(nonDitch)+ geom_sf(size=0.3)+
  facet_grid(trapType~startYear)+
  geom_sf(data=ab,fill=NA)+
  geom_sf(data=cit,col='red',size=1)+
  # geom_sf_text(data=citLabs,aes(label=NAME),col='red',size=3,nudge_y=citLabs$ynudge,nudge_x=citLabs$xnudge)+
  coord_sf(xlim=bounds[c(1,3)],ylim=bounds[c(2,4)]) + 
  ggtitle('Wetland & in-field trap distributions')+
  theme(axis.text=element_text(size=5),axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))
ggsave('./figures/maps/nonDitchMap.png',p1)

#Maps of arthropod groups (arthOrder) used to check completeness
arthNumbers <- arth %>% filter(!mismatchBTID,!mismatchBLID) %>% 
  group_by(BLID,year,arthOrder) %>% summarize(n=n()) %>% 
  left_join(select(siteSf,BLID),by='BLID') 
yearNumbers <- arthNumbers %>% 
  group_by(BLID,year) %>% summarize(n=sum(n)) %>% 
  left_join(select(siteSf,BLID),by='BLID') 

#Arth maps
p1 <- ggplot()+ geom_sf(data=yearNumbers,aes(geometry=geometry),size=0.3,shape=1,alpha=0.3)+
  geom_sf(data=arthNumbers,aes(geometry=geometry),size=0.3)+
  geom_sf(data=cit,col='red',size=1)+ #Cities
  facet_grid(arthOrder~year)+
  coord_sf(xlim=bounds[c(1,3)],ylim=bounds[c(2,4)])+
  ggtitle('Sampling and occurance distribution')+
  theme(axis.text=element_text(size=5),axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))
ggsave('./figures/maps/arthDistMap.png',p1,scale=1.5)


# Create new shapefiles for 2018 data --------------------------------

#Add ditch sites from 2015/2016/2018
for(i in c(2015,2016,2018)){
  trapSf %>% filter(startYear==i,trapLoc=='ditch') %>% 
    select(BLID,replicate,trapType) %>% unique() %>% 
    st_transform(3403) %>% 
    st_write(paste0('~/Documents/shapefiles/DitchSites_',i,'.shp'),delete_dsn=TRUE)
}
#Add infields from 2018
trapSf %>% filter(startYear==2018,trapLoc!='ditch') %>% 
  select(BLID,replicate,trapType) %>% unique() %>% 
  st_transform(3403) %>% 
  st_write('~/Documents/shapefiles/InfieldSite_2018.shp',delete_dsn=TRUE)

#Only 16 of the ditch sites from 2018 (repeats of previous years in S. Calgary) have been digitized.
#The remaining ditch sites (41) and infields (27) need wetland layers created around them

# Question 1: Spillover into crops ----------------------------------------

# How does insect abundance change with distance into cropland?






