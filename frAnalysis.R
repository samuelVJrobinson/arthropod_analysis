#FUNCTIONAL REGRESSION ANALYSIS OF IN-FIELD INSECT ABUNDANCE USING ONION-RING LANDSCAPE DATA
# SR WINTER 2020

library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
library(mgcv)
library(sf)

# Load everything ---------------------------------------------------------
load('./data/cleanData.Rdata') #Load site, trap, arth data
load('./data/geoData.Rdata') #Load NNdist and oring data

# Pterostichus ----------------------------------

tempArth <- arth %>% filter(genus=='Pterostichus',species=='melanarius') %>% group_by(BTID) %>% summarize(n=n())

#Only using pitfall traps from 2017
tempTrap <- trap %>% filter(startYear==2017,grepl('PF',BTID)) %>%
  filter(!grepl('DPF',BTID)) %>% #Some ditch sites don't have cover properly digitized
  select(BLID,BTID,pass,contains('julian'),deployedhours,trapLoc,distFrom,dist,lonTrap,latTrap,lonSite:ID) %>% 
  mutate(trapdays=deployedhours/24) %>% select(-deployedhours) %>% 
  left_join(rownames_to_column(data.frame(nnDistMat),'ID'),by='ID') %>% 
  left_join(tempArth,by='BTID') %>% mutate(n=ifelse(is.na(n),0,n)) %>% filter(!is.na(grass)) %>% 
  arrange(BLID,dist,pass) %>% mutate(BLID=factor(BLID)) %>% 
  mutate_at(vars(ephemeral:noncrop),function(x) ifelse(x>1500,1500,x)) %>% 
  rename('count'='n')

trap %>% filter(startYear==2017,grepl('PF',BTID)) %>% select(BTID) %>% data.frame()

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
#Is this a "cultural" species?
mod1 <- gam(count~offset(log(trapdays))+s(endjulian)+s(grass)+s(BLID,bs='re'),data=tempTrap,family='nb',method='ML')
gam.check(mod1)
summary(mod1)  
# par(mfrow=c(2,1))
plot(mod1,scheme=2,pages=1)

#Map of site intercepts (data scale)
tempTrap %>% select(BLID,lonSite,latSite,endjulian,trapdays:noncrop) %>% 
  mutate_at(vars(-BLID:-latSite),mean) %>% distinct() %>% 
  mutate(pred=predict(mod1,newdata=.)) %>% 
  st_as_sf(coords=c('lonSite','latSite'),crs=4326) %>% 
  ggplot()+geom_sf(aes(col=pred))

#Effect of grass
with(tempTrap,expand.grid(BLID=unique(BLID[BLID=='26168']),trapdays=7,endjulian=200,grass=0:400)) %>% 
  mutate(pred=predict(mod1,newdata=.),se=predict(mod1,newdata=.,se.fit=T)$se.fit) %>% 
  mutate(upr=pred+1.96*se,lwr=pred-1.96*se) %>% 
  mutate_at(vars(pred,upr,lwr),exp) %>%
  ggplot(aes(x=grass))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred))+labs(x='Distance from grass',y='Catches/week',title='Pterosticus')

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

mod3 <- gam(count~offset(log(trapdays))+s(endjulian)+
              dist+trapLoc+grass+
              s(BLID,bs='re'),data=tempTrap,family='nb')
summary(mod3)
plot(mod3,pages=1,all.terms=T,se=F)


with(tempTrap,expand.grid(BLID=unique(BLID[BLID=='26168']),trapdays=7,endjulian=c(180,200,220),dist=0:200)) %>% 
  mutate(pred=predict(mod3,newdata=.),se=predict(mod3,newdata=.,se.fit=T)$se.fit) %>% 
  mutate(upr=pred+1.96*se,lwr=pred-1.96*se) %>% 
  mutate_at(vars(pred,upr,lwr),exp) %>%
  mutate(doy=factor(endjulian,labels=c('early','mid','late'))) %>% 
  ggplot(aes(x=dist))+geom_ribbon(aes(ymax=upr,ymin=lwr,fill=doy),alpha=0.3)+
  geom_line(aes(y=pred,col=doy))+labs(x='Distance',y='Catches/week',title='Pterosticus')
