#PRELIMINARY DATA DISPLAY AND ANALYSIS - MAINLY FOR MESSING AROUND WITH THE DATA AND SEEING WHAT WON'T WORK
#RUN ON GALPERN LAB MACHINE
#SR WINTER 2020

#Load data and libraries
load('./data/cleanData.Rdata')
library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
library(sf)

# Abundance tables --------------------------------------------------------

#What spp of beetles are present in pitfall traps?
arth %>% filter(grepl('PF',BTID),arthOrder=='Coleoptera') %>%
  mutate(genSpp=paste(genus,species,sep=' ')) %>% group_by(family,genus,species) %>% 
  summarize(count=n()) %>% arrange(desc(count)) %>% 
  unite(genSpp,genus,species,sep=' ') %>% 
  mutate(genSpp=factor(genSpp,levels=genSpp)) %>% 
  filter(genSpp!='Pterostichus melanarius') %>% #Strip out P. melanarius
  ggplot(aes(x=genSpp,y=count))+geom_col()+
  # scale_y_sqrt()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
#Tons of Pterostichus melanarius. 

#Genera
arth %>% filter(grepl('PF',BTID),arthOrder=='Coleoptera') %>% filter(species!='melanarius') %>% 
  group_by(genus) %>% summarize(count=n()) %>% data.frame() %>% 
  arrange(desc(count)) %>% mutate(genus=factor(genus,levels=genus)) %>% 
  ggplot(aes(x=genus,y=count))+geom_col()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#What spp of spiders are present in pitfall traps?
arth %>% filter(grepl('PF',BTID),arthOrder=='Araneae') %>%
  mutate(genSpp=paste(genus,species,sep=' ')) %>% group_by(family,genus,species) %>% 
  summarize(count=n()) %>% ungroup() %>% arrange(desc(count)) %>%
  unite(genSpp,genus,species,sep=' ') %>% 
  mutate(genSpp=factor(genSpp,levels=genSpp)) %>% 
  ggplot(aes(x=genSpp,y=count))+geom_col()+ 
  # scale_y_sqrt()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
#Lots of Pardosa distincta and moesta. Could try all lycosids or Pardosa at once, but this might be a stretch

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
  #Replace trap lat/lon with site lat/lon if not recorded
  mutate(latTrap=ifelse(latTrap==0,latSite,latTrap),lonTrap=ifelse(lonTrap==0,lonSite,lonTrap)) %>% 
  st_as_sf(coords=c('lon','lat'),crs=4326) %>% 
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

#Pterostichus sampling (2017-2018)
arth %>% select(BTID:BLID,year,genus,species) %>% filter(genus=='Pterostichus',species=='melanarius') %>% 
  group_by(BTID) %>% summarize(count=n()) %>%
  right_join(filter(trapSf,startYear==2017|startYear==2018),by='BTID') %>% 
  mutate(count=ifelse(is.na(count),0,count),deployedWeeks=deployedhours/(24*7)) %>%
  filter(!is.na(deployedWeeks)) %>% 
  # ggplot()+geom_histogram(aes(x=deployedWeeks))
  ggplot()+geom_sf(aes(geometry=geometry,size=count/deployedWeeks),alpha=0.1)+
  facet_grid(trapType~startYear)

#Opiliones sampling
arth %>% select(BTID:BLID,year,genus,species) %>% filter(genus=='Phalangium',species=='opilio') %>% 
  group_by(BTID) %>% summarize(count=n()) %>%
  right_join(filter(trapSf,startYear==2017|startYear==2018),by='BTID') %>% 
  mutate(count=ifelse(is.na(count),0,count),deployedWeeks=deployedhours/(24*7)) %>%
  filter(!is.na(deployedWeeks)) %>% 
  # ggplot()+geom_histogram(aes(x=deployedWeeks))
  ggplot()+geom_sf(aes(geometry=geometry,size=count/deployedWeeks))+
  facet_grid(trapType~startYear)

#Pardosa sampling (2017 only)
arth %>% select(BTID:BLID,year,genus,species) %>% filter(genus=='Pardosa',species=='distincta') %>% 
  group_by(BTID) %>% summarize(count=n()) %>%
  right_join(filter(trapSf,startYear==2017),by='BTID') %>% 
  mutate(count=ifelse(is.na(count),0,count),deployedWeeks=deployedhours/(24*7)) %>%
  filter(!is.na(deployedWeeks)) %>% 
  # ggplot()+geom_histogram(aes(x=deployedWeeks))
  ggplot()+geom_sf(aes(geometry=geometry,size=count/deployedWeeks,col=count/deployedWeeks))+
  facet_grid(trapType~startYear)

  
  

