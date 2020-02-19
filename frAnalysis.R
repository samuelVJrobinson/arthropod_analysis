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

tempArth <- arth %>% filter(genus=='Pterostichus') %>% group_by(BTID) %>% summarize(n=n())

#Only using pitfall traps from 2017
tempTrap <- trap %>% filter(startYear==2017,grepl('PF',BTID)) %>%
  filter(!grepl('DPF',BTID)) %>% #Some ditch sites don't have cover properly digitized
  select(BLID,BTID,pass,contains('julian'),deployedhours,trapLoc,distFrom,dist,lonTrap,latTrap,lonSite:ID) %>% 
  mutate(trapdays=deployedhours/24) %>% select(-deployedhours) %>% 
  left_join(rownames_to_column(data.frame(nnDistMat),'ID'),by='ID') %>% 
  left_join(tempArth,by='BTID') %>% mutate(n=ifelse(is.na(n),0,n)) %>% filter(!is.na(grass)) %>% 
  arrange(BLID,dist,pass) %>% mutate(BLID=factor(BLID)) %>% 
  rename('count'='n')

trap %>% filter(startYear==2017,grepl('PF',BTID)) %>% select(BTID) %>% data.frame()


#Arrange oRing cover matrix to correspond with correct rows
tempORing <- lapply(oRingMat,function(x) x[match(tempTrap$ID,rownames(x)),])

#Total area of each ring in matrix form
tempRingArea <- matrix(rep((pi*seq(20,1000,20)^2)-(pi*(seq(20,1000,20)-20)^2),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))

#Matrix of distance values
tempRingMat <- matrix(rep(as.numeric(gsub('d','',colnames(tempORing[[1]]))),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))

#Model of landscape effect (nnDist only)
#Appears that distance to grass has a positive effect, and distance to noncrop a negative effect 
#Are these "cultural" species?
mod1 <- gam(count~offset(log(trapdays))+
              s(endjulian)+
              # s(wetlands)+
              # s(trees)+
              s(grass)+
              # s(shrubs)+
              s(BLID,bs='re'),
            data=tempTrap,family='nb',method='ML')

AIC(mod1,
    update(mod1,.~.-s(grass)),
    update(mod1,.~.-s(wetlands)))
# gam.check(mod1)
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
tempTrap %>% select(BLID,lonSite,latSite,endjulian,trapdays:noncrop) %>% 
  mutate(trapdays=7) %>% 
  mutate_at(vars(-BLID:-latSite,-trapdays,-grass),mean) %>% distinct() %>% 
  filter(BLID=='26168') %>% 
  mutate(pred=predict(mod1,newdata=.),se=predict(mod1,newdata=.,se.fit=T)$se.fit) %>% 
  mutate(upr=pred+1.96*se,lwr=pred-1.96*se) %>% 
  mutate_at(vars(pred,upr,lwr),exp) %>%
  ggplot(aes(x=grass))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred))



#Model of landscape (oRingDist)
grassMat <- tempORing$grass #m2 cover within each ring
noncropMat <- tempORing$noncrop
grassMatPerc <- grassMat/tempRingArea #prop cover within each ring
noncropMatPerc <- noncropMat/tempRingArea #prop cover within each ring

mod2 <- gam(count~offset(log(trapdays))+s(endjulian)+
              s(tempRingMat,by=grassMatPerc)+
              # s(tempRingMat,by=noncropMatPerc)+
              s(BLID,bs='re'),data=tempTrap,family='nb',method='ML')

AIC(mod2,
    update(mod2,.~.-s(tempRingMat,by=grassMatPerc)),
    update(mod2,.~.-s(grass)))

summary(mod2)
plot(mod2,scheme=2,pages=1)



