#FUNCTIONAL REGRESSION ANALYSIS OF IN-FIELD INSECT ABUNDANCE USING ONION-RING LANDSCAPE DATA
# SR WINTER 2020

library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
library(mgcv)
library(sf)
library(gstat)

# Load everything ---------------------------------------------------------
load('./data/cleanData.Rdata') #Load site, trap, arth data
load('./data/geoData.Rdata') #Load NNdist and oring data

# Pterostichus ----------------------------------

#What spp of beetles are present?
arth %>% filter(grepl('PF',BTID),arthOrder=='Coleoptera') %>%
  mutate(genSpp=paste(genus,species,sep=' ')) %>% group_by(family,genus,species) %>% 
  summarize(n=n()) %>% arrange(genus,desc(n)) %>% data.frame()

arth %>% filter(grepl('PF',BTID),arthOrder=='Coleoptera') %>%
  group_by(genus) %>% summarize(n=n()) %>% data.frame() %>% 
  arrange(desc(n))

tempArth <- arth %>% filter(genus=='Pterostichus',species=='melanarius') %>% group_by(BTID) %>% summarize(n=n())

#Only using pitfall traps from 2017
tempTrap <- trap %>% filter(startYear==2017,grepl('PF',BTID)) %>%
  # filter(!grepl('PF',BTID)) %>% #Some ditch sites don't have cover properly digitized at further distances, but it's OK for now
  select(BLID,BTID,pass,contains('julian'),deployedhours,trapLoc,distFrom,dist,lonTrap,latTrap,lonSite:ID) %>% 
  mutate(trapdays=deployedhours/24) %>% select(-deployedhours) %>% 
  left_join(rownames_to_column(data.frame(nnDistMat),'ID'),by='ID') %>% 
  left_join(tempArth,by='BTID') %>% mutate(n=ifelse(is.na(n),0,n)) %>% filter(!is.na(grass)) %>% 
  arrange(BLID,dist,pass) %>% mutate(BLID=factor(BLID)) %>% 
  mutate_at(vars(ephemeral:noncrop),function(x) ifelse(x>1500,1500,x)) %>% 
  rename('count'='n') %>% 
  #Converts 0 dist trapLoc to distFrom - pitfall at 0 m are "in" feature
  mutate(trapLoc=factor(ifelse(dist==0 & distFrom!='control',distFrom,trapLoc))) %>% 
  #Get UTM coordinates for traps
  mutate(lonTrap2=lonTrap,latTrap2=latTrap) %>% #Duplicate columns
  st_as_sf(coords=c('lonTrap2','latTrap2'),crs=4326) %>% 
  st_transform(3403) %>% mutate(easting=st_coordinates(.)[,1],northing=st_coordinates(.)[,2]) %>% 
  mutate_at(vars(easting,northing),scale) #Scale coordinates to be a reasonable range

#Arrange oRing cover matrix to correspond with correct rows
tempORing <- lapply(oRingMat,function(x) x[match(tempTrap$ID,rownames(x)),])
#Total area of each ring in matrix form
tempRingArea <- matrix(rep((pi*seq(20,1000,20)^2)-(pi*(seq(20,1000,20)-20)^2),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
#Matrix of distance values
tempRingMat <- matrix(rep(as.numeric(gsub('d','',colnames(tempORing[[1]]))),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
grassMat <- tempORing$grass/10000 #hectares grass cover within each ring
noncropMat <- tempORing$noncrop/10000 #Noncrop cover (very similar to grass)
grassMatPerc <- grassMat/tempRingArea #prop cover within each ring
noncropMatPerc <- noncropMat/tempRingArea #prop cover within each ring

#Model of landscape effect (nnDist only)
#Appears that distance to grass has a positive effect, and distance to noncrop a negative effect 
#Spatio-temporal tensor doesn't fit as well as straightforward random effects. Investigate this later.
#When using trap location, canola clearly has the most 
#Is this a "cultural" species?
mod1 <- gam(count~offset(log(trapdays))+s(endjulian)+
              trapLoc+
              s(BLID,bs='re')
              # s(grass)+s(wetlands)+
            ,
            data=tempTrap,family='nb',method='REML')
gam.check(mod1)
summary(mod1)  
plot(mod1,scheme=2,pages=1,all.terms=F,residuals=F)

# #Trick to do multiple comparisons with GAM
# library(multcomp)
# modMat1 <- glht(lm(count~trapLoc,data=tempTrap),linfct=mcp(trapLoc='Tukey'))$linfct
# modMat2 <- model.matrix(mod1)[c(1:nrow(modMat1)),]
# modMat2[,] <- 0
# modMat2[c(1:nrow(modMat1)),c(1:ncol(modMat1))] <- modMat1
# rownames(modMat2) <- rownames(modMat1)
# modTest <- glht(mod1,linfct=modMat2)
# summary(modTest)
# plot(modTest)
# cld(modTest) #Doesn't work

#Random BLID intercept with value closest to zero
midBLID <- levels(tempTrap$BLID)[which.min(abs(coef(mod1)[grepl('BLID',names(coef(mod1)))]))]
#Effect of trapLoc
p1 <- with(tempTrap,expand.grid(BLID=unique(BLID[BLID==midBLID]),trapdays=7,
                 endjulian=mean(endjulian),grass=mean(grass),
                 wetlands=mean(wetlands),
                 trapLoc=unique(trapLoc))) %>% 
  mutate(pred=predict(mod1,newdata=.),se=predict(mod1,newdata=.,se.fit=T)$se.fit) %>% 
  mutate(upr=pred+se,lwr=pred-se) %>% 
  mutate_at(vars(pred,upr,lwr),exp) %>%
  ggplot(aes(x=trapLoc))+geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Catches/week',title='P. melanarius ~ offset + s(time) + s(BLID,type=re) + trapLoc')
ggsave('./figures/01_P_melanarius_trapLoc.png',p1,height=6,width=6)

p2 <- tempTrap %>% select(BLID,lonSite,latSite) %>% st_drop_geometry() %>% 
  # mutate(trapLoc=ifelse(trapLoc=='ditch',trapLoc,'inField')) %>% 
  distinct() %>%
  mutate(ranef=coef(mod1)[grepl('BLID',names(coef(mod1)))]) %>%
  st_as_sf(coords=c('lonSite','latSite'),crs=4326) %>% 
  ggplot(aes(col=ranef))+
  geom_sf(aes(size=abs(ranef)*0.4),alpha=0.5,show.legend=T)+
  # facet_wrap(~trapLoc,ncol=1)+
  scale_colour_gradient(low='blue',high='red')+
  scale_size(guide='none')+
  labs(title='P. melanarius random intercepts')
ggsave('./figures/01_P_melanarius_intercepts.png',p2,height=6,width=8)


#Model of landscape effect (ring composition)
mod2 <- gam(count~offset(log(trapdays))+s(endjulian)+s(tempRingMat,by=grassMat)+s(BLID,bs='re'),data=tempTrap,family='nb',method='ML')

#Looks like distance is the only important factor in this
AIC(mod1,mod2,update(mod2,.~.-s(tempRingMat,by=grassMat)+s(tempRingMat,by=grassMatPerc)),update(mod2,.~.+s(grass)))
summary(mod2)
plot(mod2,scheme=2,pages=1)
gam.check(mod2)

plot(mod1$model$count,fitted(mod1),ylab='Predicted',xlab='Actual',pch=19,cex=0.5)
points(mod2$model$count,fitted(mod2),pch=19,cex=0.5,col='red')
abline(0,1,lty='dashed')
legend('bottomright',c('Distance','Ring comp'),fill=c('black','red'))

#Summary for Pterostichus melanarius: appears to be more individuals at the centre of the fields

#Messing around with data display
ggplot(tempTrap,aes(x=midjulian,y=count+1,col=distFrom))+
  geom_point(aes(shape=trapLoc))+
  facet_wrap(~BLID)+
  geom_smooth(method='gam',se=F,aes(lty=trapLoc))+
  scale_y_log10()+
  theme(legend.position='bottom')+
  labs(x='DOY',col='Feature')

ggplot(tempTrap,aes(x=dist,y=count+1,col=distFrom))+
  geom_point(aes(shape=trapLoc))+
  facet_wrap(~BLID)+
  geom_smooth(method='gam',se=F,aes(lty=trapLoc))+
  scale_y_log10()+
  # theme(legend.position='bottom')+
  labs(x='Dist',col='Distance from',lty='Trap location',shape='Trap location')

# #This model doesn't make any sense, because distance means different things depending on whether trap was in "normal" field (control) or in a field near a pivot corner or wetland. Better to use above model
# mod3 <- gam(count~offset(log(trapdays))+s(endjulian)+
#               dist*distFrom+trapLoc+
#               s(BLID,bs='re'),data=tempTrap,family='nb')
# summary(mod3) #Once again, distance from the edge seems to be the biggest fixed effect predictor for P.melanarius. Larg
# plot(mod3,all.terms=T,pages=1)
# 
# #Random intercept plots - doesn't look too bad
# tempTrap %>% select(BLID,lonSite,latSite) %>% distinct() %>% 
#   mutate(ranef=coef(mod3)[grepl('BLID',names(coef(mod3)))]) %>% 
#   st_as_sf(coords=c('lonSite','latSite'),crs=4326) %>% st_transform(3403) %>% 
#   ggplot()+geom_sf(aes(size=abs(ranef),col=ranef))+
#   scale_colour_gradient(low='blue',high='red')
# 
# #Distance-only variogram - need
# tempTrap %>% select(BLID,lonTrap,latTrap) %>%
#   mutate(resid=resid(mod3)) %>%
#   st_as_sf(coords=c('lonTrap','latTrap'),crs=4326) %>%
#   st_transform(3403) %>%
#   variogram(resid~1,data=.) %>% plot()
# 
# #Spatio-temporal variogram for residuals - need to look into this more
# library(sp)
# library(spacetime)
#  
# temp <- tempTrap %>% data.frame()
# coordinates(temp) <- ~lonTrap+latTrap
# proj4string(temp) <- CRS("+init=epsg:4326") #Convert to SpatialPointsDataFrame
# temp <- spTransform(temp,CRS("+init=epsg:3403")) #Convert to UTM
# temp <- SpatialPoints(temp@coords,CRS("+init=epsg:3403"))
# 
# timeDF <- STIDF(sp=temp,
#                 time=as.POSIXct(paste(tempTrap$endjulian,'2017',sep='-'),format='%j-%Y'),
#                 data=data.frame(resid=resid(mod3)))
# stplot(timeDF)
# stVar <- variogramST(resid~1,data=timeDF,tunit='days',assumeRegular=T) #Fit spatiotemporal variogram
# plot(stVar,map=T)


# Wolf spiders -----------------------------------------------------------------

tempArth <- arth %>% filter(family=='Lycosidae') %>% group_by(BTID) %>% summarize(n=n())

#Only using pitfall traps from 2017
tempTrap <- trap %>% filter(startYear==2017,grepl('PF',BTID)) %>%
  # filter(!grepl('PF',BTID)) %>% #Some ditch sites don't have cover properly digitized at further distances, but it's OK for now
  select(BLID,BTID,pass,contains('julian'),deployedhours,trapLoc,distFrom,dist,lonTrap,latTrap,lonSite:ID) %>% 
  mutate(trapdays=deployedhours/24) %>% select(-deployedhours) %>% 
  left_join(rownames_to_column(data.frame(nnDistMat),'ID'),by='ID') %>% 
  left_join(tempArth,by='BTID') %>% mutate(n=ifelse(is.na(n),0,n)) %>% filter(!is.na(grass)) %>% 
  arrange(BLID,dist,pass) %>% mutate(BLID=factor(BLID)) %>% 
  mutate_at(vars(ephemeral:noncrop),function(x) ifelse(x>1500,1500,x)) %>% 
  rename('count'='n') %>% 
  #Converts 0 dist trapLoc to distFrom - pitfall at 0 m are "in" feature
  mutate(trapLoc=ifelse(dist==0 & distFrom!='control',distFrom,trapLoc)) %>% 
  #Get UTM coordinates for traps
  mutate(lonTrap2=lonTrap,latTrap2=latTrap) %>% #Duplicate columns
  st_as_sf(coords=c('lonTrap2','latTrap2'),crs=4326) %>% 
  st_transform(3403) %>% mutate(easting=st_coordinates(.)[,1],northing=st_coordinates(.)[,2]) %>% 
  mutate_at(vars(easting,northing),scale) #Scale coordinates to be a reasonable range

#Arrange oRing cover matrix to correspond with correct rows
tempORing <- lapply(oRingMat,function(x) x[match(tempTrap$ID,rownames(x)),])
#Total area of each ring in matrix form
tempRingArea <- matrix(rep((pi*seq(20,1000,20)^2)-(pi*(seq(20,1000,20)-20)^2),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
#Matrix of distance values
tempRingMat <- matrix(rep(as.numeric(gsub('d','',colnames(tempORing[[1]]))),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
grassMat <- tempORing$grass/10000 #hectares grass cover within each ring
noncropMat <- tempORing$noncrop/10000 #Noncrop cover (very similar to grass)
grassMatPerc <- grassMat/tempRingArea #prop cover within each ring
noncropMatPerc <- noncropMat/tempRingArea #prop cover within each ring

mod1 <- gam(count~offset(log(trapdays))+s(endjulian)+
              s(BLID,bs='re')+
              # s(grass)+s(wetlands)+
              trapLoc,
            data=tempTrap,family='nb',method='REML')
gam.check(mod1)
summary(mod1)  
plot(mod1,scheme=2,pages=1,all.terms=F,residuals=F)

#Random BLID intercept with value closest to zero
midBLID <- levels(tempTrap$BLID)[which.min(abs(coef(mod1)[grepl('BLID',names(coef(mod1)))]))]
#Effect of trapLoc
p1 <- with(tempTrap,expand.grid(BLID=unique(BLID[BLID==midBLID]),trapdays=7,
                          endjulian=mean(endjulian),grass=mean(grass),
                          wetlands=mean(wetlands),
                          trapLoc=unique(trapLoc))) %>% 
  mutate(pred=predict(mod1,newdata=.),se=predict(mod1,newdata=.,se.fit=T)$se.fit) %>% 
  mutate(upr=pred+se,lwr=pred-se) %>% 
  mutate_at(vars(pred,upr,lwr),exp) %>%
  ggplot(aes(x=trapLoc))+geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Catches/week',title='Lycosidae ~ offset + s(time) + s(BLID,type=re) + trapLoc')
ggsave('./figures/02_Lycosidae_trapLoc.png',p1,height=6,width=6)

p2 <- tempTrap %>% select(BLID,lonSite,latSite) %>% st_drop_geometry() %>% 
  # mutate(trapLoc=ifelse(trapLoc=='ditch',trapLoc,'inField')) %>% 
  distinct() %>%
  mutate(ranef=coef(mod1)[grepl('BLID',names(coef(mod1)))]) %>%
  st_as_sf(coords=c('lonSite','latSite'),crs=4326) %>% 
  ggplot()+geom_sf(aes(col=ranef,size=abs(ranef)),alpha=0.5,show.legend=F)+
  # facet_wrap(~trapLoc,ncol=1)+
  scale_colour_gradient(low='blue',high='red')+
  labs(title='Lycosidae random intercepts')
ggsave('./figures/02_Lycosidae_intercepts.png',p2,height=6,width=8)


# Harvestmen -------------------------------------------------------------

tempArth <- arth %>% filter(arthOrder=='Opiliones') %>% group_by(BTID) %>% summarize(n=n())

#Only using pitfall traps from 2017
tempTrap <- trap %>% filter(startYear==2017,grepl('PF',BTID)) %>%
  # filter(!grepl('PF',BTID)) %>% #Some ditch sites don't have cover properly digitized at further distances, but it's OK for now
  select(BLID,BTID,pass,contains('julian'),deployedhours,trapLoc,distFrom,dist,lonTrap,latTrap,lonSite:ID) %>% 
  mutate(trapdays=deployedhours/24) %>% select(-deployedhours) %>% 
  left_join(rownames_to_column(data.frame(nnDistMat),'ID'),by='ID') %>% 
  left_join(tempArth,by='BTID') %>% mutate(n=ifelse(is.na(n),0,n)) %>% filter(!is.na(grass)) %>% 
  arrange(BLID,dist,pass) %>% mutate(BLID=factor(BLID)) %>% 
  mutate_at(vars(ephemeral:noncrop),function(x) ifelse(x>1500,1500,x)) %>% 
  rename('count'='n') %>% 
  #Converts 0 dist trapLoc to distFrom - pitfall at 0 m are "in" feature
  mutate(trapLoc=ifelse(dist==0 & distFrom!='control',distFrom,trapLoc)) %>% 
  #Get UTM coordinates for traps
  mutate(lonTrap2=lonTrap,latTrap2=latTrap) %>% #Duplicate columns
  st_as_sf(coords=c('lonTrap2','latTrap2'),crs=4326) %>% 
  st_transform(3403) %>% mutate(easting=st_coordinates(.)[,1],northing=st_coordinates(.)[,2]) %>% 
  mutate_at(vars(easting,northing),scale) #Scale coordinates to be a reasonable range

#Arrange oRing cover matrix to correspond with correct rows
tempORing <- lapply(oRingMat,function(x) x[match(tempTrap$ID,rownames(x)),])
#Total area of each ring in matrix form
tempRingArea <- matrix(rep((pi*seq(20,1000,20)^2)-(pi*(seq(20,1000,20)-20)^2),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
#Matrix of distance values
tempRingMat <- matrix(rep(as.numeric(gsub('d','',colnames(tempORing[[1]]))),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
grassMat <- tempORing$grass/10000 #hectares grass cover within each ring
noncropMat <- tempORing$noncrop/10000 #Noncrop cover (very similar to grass)
grassMatPerc <- grassMat/tempRingArea #prop cover within each ring
noncropMatPerc <- noncropMat/tempRingArea #prop cover within each ring

mod1 <- gam(count~offset(log(trapdays))+s(endjulian)+
              s(BLID,bs='re')+
              # s(grass)+s(wetlands)+
              trapLoc,
            data=tempTrap,family='nb',method='REML')
gam.check(mod1)
summary(mod1)  
plot(mod1,scheme=2,pages=1,all.terms=F,residuals=F)

#Random BLID intercept with value closest to zero
midBLID <- levels(tempTrap$BLID)[which.min(abs(coef(mod1)[grepl('BLID',names(coef(mod1)))]))]
#Effect of trapLoc
p1 <- with(tempTrap,expand.grid(BLID=unique(BLID[BLID==midBLID]),trapdays=7,
                          endjulian=mean(endjulian),grass=mean(grass),
                          wetlands=mean(wetlands),
                          trapLoc=unique(trapLoc))) %>% 
  mutate(pred=predict(mod1,newdata=.),se=predict(mod1,newdata=.,se.fit=T)$se.fit) %>% 
  mutate(upr=pred+se,lwr=pred-se) %>% 
  mutate_at(vars(pred,upr,lwr),exp) %>%
  ggplot(aes(x=trapLoc))+geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Catches/week',title='Opiliones ~ offset + s(time) + s(BLID,type=re) + trapLoc')
ggsave('./figures/03_Opiliones_trapLoc.png',p1,height=6,width=6)

p2 <- tempTrap %>% select(BLID,lonSite,latSite) %>% st_drop_geometry() %>% 
  # mutate(trapLoc=ifelse(trapLoc=='ditch',trapLoc,'inField')) %>% 
  distinct() %>%
  mutate(ranef=coef(mod1)[grepl('BLID',names(coef(mod1)))]) %>%
  st_as_sf(coords=c('lonSite','latSite'),crs=4326) %>% 
  ggplot()+geom_sf(aes(col=ranef,size=abs(ranef)),alpha=0.5,show.legend=F)+
  # facet_wrap(~trapLoc,ncol=1)+
  scale_colour_gradient(low='blue',high='red')+
  labs(title='Opiliones random intercepts')
ggsave('./figures/03_Opiliones_intercepts.png',p2,height=6,width=8)

