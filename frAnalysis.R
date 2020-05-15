#FUNCTIONAL REGRESSION ANALYSIS OF IN-FIELD INSECT ABUNDANCE USING ONION-RING LANDSCAPE DATA
# SR WINTER 2020

library(tidyverse)
library(ggpubr)
theme_set(theme_classic())
maptheme <- theme(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
library(mgcv)
library(sf)
library(gstat)
library(beepr)

# Load everything ---------------------------------------------------------
load('./data/cleanData.Rdata') #Load site, trap, arth data
load('./data/geoData.Rdata') #Load NNdist and oring data
load('./data/geoDataAAFC.Rdata') #Load oring data extracted from AAFC data

source('helperFunctions.R') #Load helper functions

#Arrange oRing data
oRingMat2 <- lapply(oRingMat2,function(x) x[,-1]) #Remove distance=0 measurements (already in trapLoc category)

#Convert to proportion area within each ring
oRingMat2Prop <-  lapply(oRingMat2,function(x) x/Reduce('+',oRingMat2))

#Add total proportion of "noncrop" (not sure about "Pasture", but this is likely untilled/unplanted during that year)
nonCropClasses <- c('Grassland','Pasture','Wetland','Urban','Shrubland','Forest','Water','Barren','Fallow')
oRingMat2Prop$NonCrop <- Reduce('+',oRingMat2Prop[names(oRingMat2Prop) %in% nonCropClasses])

#Combine classes that are collinear/represent similar things

#Trees + Shrubs
oRingMat2Prop$TreeShrub <- oRingMat2Prop$Shrubland + oRingMat2Prop$Forest
# oRingMat2Prop$Shrubland <- NULL; oRingMat2Prop$Forest <- NULL
# #Water + Wetland
# oRingMat2Prop$WetlandWater <- oRingMat2Prop$Wetland + oRingMat2Prop$Water
# oRingMat2Prop$Wetland <- NULL; oRingMat2Prop$Water <- NULL

#Create model formulas
modFormulas <- 'count~offset(log(trapdays))+s(day,k=10,bs=basisFun)+s(E,N,k=50,bs=basisFun)' #Temporal + spatial
# modFormulas <- paste0(modFormulas,'+ti(N,E,day,k=5,bs=basisFun)') # Add spatiotemporal interaction
modFormulas[c(2,3,4)] <- paste0(modFormulas[1],'+trapLoc-1')
# mod3Vars <- c('Grassland','Canola','Pasture','Wetland','TreeShrub','Pulses','Urban') #Variables for mod3
mod3Vars <- c('Grassland','Cereal','Canola','Pasture','Wetland','Forest','Shrubland','Pulses','Urban') #Variables for mod3
#Cereal may be causing problems in estimation - collinear with canola
for(i in 1:length(mod3Vars)){ #Add in specified terms (s and ti)
  #Main effects + interaction (s + ti)
  modFormulas[3] <- paste0(modFormulas[3],"+ s(distMat,by=",mod3Vars[i],",bs=basisFun)")
  modFormulas[3] <- paste0(modFormulas[3],"+ s(endDayMat,by=",mod3Vars[i],",bs=basisFun)")
  modFormulas[3] <- paste0(modFormulas[3],"+ ti(distMat,endDayMat,by=",mod3Vars[i],",bs=basisFun)")
  # #Full interaction (te)
  # modFormulas[3] <- paste0(modFormulas[3],"+ te(distMat,endDayMat,by=",mod3Vars[i],",bs=basisFun)")
}
modFormulas[4] <- paste0(modFormulas[2],'+s(distMat,by=NonCrop,bs=basisFun)+s(endDayMat,by=NonCrop,bs=basisFun)+ti(distMat,endDayMat,by=NonCrop,bs=basisFun)')

# Pterostichus melanarius ----------------------------------

#What spp of beetles are present?
arth %>% filter(grepl('PF',BTID),arthOrder=='Coleoptera') %>%
  mutate(genSpp=paste(genus,species,sep=' ')) %>% group_by(family,genus,species) %>% 
  summarize(n=n()) %>% arrange(genus,desc(n)) %>% data.frame()

arth %>% filter(grepl('PF',BTID),arthOrder=='Coleoptera') %>%
  group_by(genus) %>% summarize(n=n()) %>% data.frame() %>% 
  arrange(desc(n))

#Select only P. melanarius
tempArth <- arth %>% filter(genus=='Pterostichus',species=='melanarius') %>% group_by(BTID) %>% summarize(n=n())

PteMelMod <- runMods(tempArth,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
# save(PteMelMod,file='./data/PteMelMod.Rdata')
load('./data/PteMelMod.Rdata')

#Check models
attach(PteMelMod)
# detach(PteMelMod)
AIC(mod1,mod2,mod3,mod4) #Best model has separate land cover types + SpatioTemporal effect

#Model 3 - landscape effects matter quite a bit
summary(mod3); AIC(mod3)
anova(mod3)

#Variogram of residuals - no pattern
tempTrap %>% mutate(resid=resid(mod3)) %>% 
  # group_by(ID) %>% summarize(resid=sum(abs(resid))) %>% ungroup() %>% 
  as_Spatial() %>% variogram(resid~1,.) %>% plot()

#Check k values
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F,seWithMean=T)
par(mfrow=c(3,1)); for(i in 12:14) plot(mod3,scheme=2,shade=T,rug=F,seWithMean=T,select=i); par(mfrow=c(1,1))

#Check for concurvity between smoothers
concurvity(mod3) #High concurvity
checkMC <- concurvity(mod3,full=F)$worst
termNames <- rownames(checkMC)
# termNames <- gsub('s\\(distMat(,endDayMat)?\\)','s',gsub('\\:oRingMat2Prop\\$',':',rownames(checkMC)))
# termNames <- gsub('ti\\(distMat(,endDayMat)?\\)','ti',termNames)
# termNames <- gsub('s\\(endjulian\\)','s:time',gsub('s\\(easting,northing\\)','s:space',termNames))
# termNames <- gsub('ti\\(northing,easting,endjulian\\)','ti:spacetime',termNames)
matrixplot(checkMC,mar=c(1, 10, 10, 1),termNames)
matrixplot(checkMC[4:nrow(checkMC),4:nrow(checkMC)],mar=c(1, 10, 10, 1),termNames[4:nrow(checkMC)]) #Only landscape terms
abline(h=seq(0-1/length(4:nrow(checkMC))/2,1+1/length(4:nrow(checkMC))/2,length.out=1+length(4:nrow(checkMC))/3)) #Lines to separate terms
abline(v=seq(0-1/length(4:nrow(checkMC))/2,1+1/length(4:nrow(checkMC))/2,length.out=1+length(4:nrow(checkMC))/3))

#High concurvity seems to be coming from similar terms (e.g. s(dist):Canola & s(day):Canola), but not really between terms

#Problem terms:
#Cereal ~ Canola (bad), 
#Wetland ~ Grassland (bad)
#Urban ~ Grassland,Wetland,Canola (not as bad)

#Check for correlation in smoothers
matrixplot(abs(cov2cor(sp.vcov(mod3))),c(names(mod3$sp),'scale'),mar=c(1,10,10,1))

varcomp <- gam.vcomp(mod3) #Large amount of variation explained by easting

varcomp %>% data.frame() %>% rownames_to_column('term') %>%
  mutate(term=factor(term,levels=term)) %>% 
  ggplot(aes(x=term,y=std.dev))+geom_col()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#Problem: need partial effects plots for function regression of space, time, or both.
# First two are doable, but interactions are tricky
#Solution: get predictions for new data frame
# Generate prediction matrix using PredictMat(smootherTerm,newdata)
# Multiply by coefs to get predictions
# If ti() is significant, add predictions together to get overall predictions
#New problem: how to get SEs? Simulation? Figure out how Wood does it.
#Solution: se = sqrt(New prediction matrix %*% Covariance matrix for predictors * New prediction matrix)

#Temporal smoother
p1 <- data.frame(day=min(datList$day):max(datList$day)) %>% 
  smoothPred(mod3,whichSmooth=1) %>% 
  ggplot(aes(x=day,y=pred))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+
  labs(x='Day of Year',y='Activity')

#Spatial smoother
p2 <- with(datList,expand.grid(E=seq(from=min(E)-5,to=max(E)+5,length.out=100),N=seq(from=min(N)-5,to=max(N)+5,length.out=100))) %>%
  smoothPred(mod3,whichSmooth=2) %>% 
  filter(!exclude.too.far(E,N,datList$E,datList$N,0.1)) %>% 
  ggplot(aes(x=E,y=N))+geom_raster(aes(fill=pred))+
  geom_point(data=tempTrap,aes(x=easting,y=northing))+
  scale_fill_gradient(low='blue',high='red')+
  labs(x='Easting (km)',y='Northing (km)',fill='Activity')+
  theme(legend.position = c(0.85, 0.2),legend.background = element_rect(size=0.5,linetype="solid",colour="black"))+
  maptheme

raneffPlot <- ggarrange(p1,p2,labels=letters[1:2]) #Plot spatial/temporal effects
ggsave('./figures/Pterostichus_melanarius_raneff.png',raneffPlot,width=8,height=4,scale=2)

#Plot significant landscape effects

#Days to display on plots (early,mid,late)
dispDays <- data.frame(doy=c(173,203,232)) %>% 
  mutate(date=c('June 20','July 20','August 20')) #Actually June/July 22, but close enough...

#Trap location 
p1 <- data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod3,type='terms',terms='trapLoc',se.fit=T)) %>% 
  rename('pred'='pred.trapLoc','se'='pred.trapLoc.1') %>% distinct() %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96) %>% 
  mutate(trapLoc=factor(trapLoc,labels=firstUpper(levels(trapLoc)))) %>% 
  ggplot(aes(x=trapLoc,y=pred))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Activity')

p2 <- data.frame(distMat=seq(30,1500,30),Pasture=1,y=NA) %>% #Pasture distance: p=0.018
  smoothPred(mod3,whichSmooth=9) %>% rename(x=distMat) %>% 
  effectPlot(leg=F)+labs(x='Distance (m)',y='Pasture effect')

p3 <- data.frame(endDayMat=149:241,Wetland=1,y=NA) %>% #Wetland time: p=0.00160
  smoothPred(mod3,whichSmooth=13) %>% rename(x=endDayMat) %>% 
  mutate(x=as.Date(paste0(x,'-2017'),format='%j-%Y')) %>% 
  effectPlot(leg=F)+labs(x='Time of year',y='Wetland effect')

#Tree/shrub effect
p4 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),TreeShrub=1) %>% 
  smoothPred(mod3,whichSmooth=15:17) %>% rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot()+labs(x='Distance (m)',y='Tree/Shrub effect')

#Pulses: p=0.026
p5 <- data.frame(distMat=seq(30,1500,30),Pulses=1,y=NA) %>% 
  smoothPred(mod3,whichSmooth=18) %>% rename(x=distMat) %>% 
  effectPlot(leg=F)+labs(x='Distance (m)',y='Pulse effect')

#Urban effect
p6 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Urban=1) %>% 
  smoothPred(mod3,whichSmooth=21:23) %>% 
  mutate(endDayMat=factor(endDayMat)) %>% 
  rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot()+labs(x='Distance (m)',y='Urban effect')

#Plot landscape effects
fixeffPlot <- ggarrange(p1,p2,p3,p4,p5,p6,
                        labels=letters[1:6],nrow=2,ncol=3,
                        legend='bottom',common.legend=T) 
ggsave('./figures/Pterostichus_melanarius_fixeff.png',fixeffPlot,width=8,height=4,scale=2)
rm(p1,p2,p3,p4,p5,p6,raneffPlot,fixeffPlot) #Cleanup

#Looks like the ring model of landscape does better. Important landscape features seem to be:
# Urban (spatial + temporal), Pasture, Pulses (weak), Tree/Shrubs (weak)
detach(PteMelMod)

# Pardosa distincta (wolf spider) -----------------------------------------------------------------

#What spp of spiders are present?
arth %>% filter(grepl('PF',BTID),arthOrder=='Araneae') %>%
  mutate(genSpp=paste(genus,species,sep=' ')) %>% group_by(family,genus,species) %>% 
  summarize(n=n()) %>% arrange(family,genus,desc(n)) %>% data.frame()
#Lots of Pardosa distincta and P. moesta. Could try all lycosids at once, but this might be a stretch

# #Select only wolf spiders
# tempArth <- arth %>% filter(family=='Lycosidae') %>% group_by(BTID) %>% summarize(n=n())

#Select only Pardosa distincta
tempArth <- arth %>% filter(genus=='Pardosa',species=='distincta') %>% group_by(BTID) %>% summarize(n=n())

# #Takes way longer to run. 5-10 mins +
# ParDisMod <- runMods(tempArth,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
# save(ParDisMod,file='./data/ParDisMod.Rdata')
load('./data/ParDisMod.Rdata')

#Check models
attach(ParDisMod)

AIC(mod1,mod2,mod3,mod4) 

#Model 3 
summary(mod3); AIC(mod3)
plot(mod3,pages=1,scheme=2,rug=F,shade=T,all.terms=T)

#Check for multicollinearity
checkMC <- concurvity(mod3,full=F)$worst
#Looks OK. Some correlation between grassland and wetland
termNames <- gsub('(te|s)\\(distMat(,endDayMat)?\\)','f',gsub('\\:oRingMat2Prop\\$',':',rownames(checkMC)))
matrixplot(checkMC,termNames,mar=c(1,6,6,1))

#Smoothing term estimates
matrixplot(abs(cov2cor(sp.vcov(mod3))),c(names(mod3$sp),'scale'),numSize=1,mar=c(1,8,8,1))

matrixplot(checkMC,mar=c(1, 10, 10, 1),termNames)
matrixplot(checkMC[4:nrow(checkMC),4:nrow(checkMC)],mar=c(1, 10, 10, 1),termNames[4:nrow(checkMC)]) #Only landscape terms

#High concurvity seems to be coming from similar terms (e.g. s(dist):Canola & s(day):Canola), but not really between terms
plot(mod3$sp)
round(diag(sp.vcov(mod3)))


#Variance component
varcomp <- gam.vcomp(mod3) #Large amount of variation explained by easting
round(varcomp,3)

#Temporal smoother
p1 <- data.frame(day=min(datList$day):max(datList$day)) %>% 
  smoothPred(mod3,whichSmooth=1) %>% 
  ggplot(aes(x=day,y=pred))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+
  labs(x='Day of Year',y='Activity')

#Spatial smoother
p2 <- with(datList,expand.grid(E=seq(from=min(E),to=max(E),length.out=100),N=seq(from=min(N),to=max(N),length.out=100))) %>% 
  smoothPred(mod3,whichSmooth=2) %>% 
  filter(!exclude.too.far(E,N,datList$E,datList$N,0.1)) %>% 
  ggplot(aes(x=E,y=N))+geom_raster(aes(fill=pred))+
  geom_point(data=tempTrap,aes(x=easting,y=northing))+
  scale_fill_gradient(low='blue',high='red')+
  labs(x='Easting (km)',y='Northing (km)',fill='Activity')+
  theme(legend.position = c(0.85, 0.2),legend.background = element_rect(size=0.5,linetype="solid",colour="black"))+
  maptheme

raneffPlot <- ggarrange(p1,p2,labels=letters[1:2]) #Plot spatial/temporal effects
ggsave('./figures/Pardosa_distincta_raneff.png',raneffPlot,width=8,height=4,scale=2)

#Plot important landscape effects

#Days to display on plots (early,mid,late)
dispDays <- data.frame(doy=c(173,203,232)) %>% 
  mutate(date=c('June 20','July 20','August 20')) #Actually June/July 22, but close enough...

#Trap location - far less in canola
p1 <- data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod3,type='terms',terms='trapLoc',se.fit=T)) %>% 
  rename('pred'='pred.trapLoc','se'='pred.trapLoc.1') %>% distinct() %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96) %>% 
  mutate(trapLoc=factor(trapLoc,labels=firstUpper(levels(trapLoc)))) %>% 
  ggplot(aes(x=trapLoc,y=pred))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Activity')

#Grassland
p2 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Grassland=1) %>% 
  smoothPred(mod3,whichSmooth=3:5) %>% rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot()+labs(x='Distance (m)',y='Grassland effect')

#Canola
p3 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Canola=1) %>% 
  smoothPred(mod3,whichSmooth=6:8) %>% rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot()+labs(x='Distance (m)',y='Canola effect')

#Pulses
p4 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Pulses=1) %>% 
  smoothPred(mod3,whichSmooth=18:20) %>% rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot()+labs(x='Distance (m)',y='Pulse effect')

#Urban effect
p5 <- expand.grid(distMat=seq(30,1500,30),Urban=1,y=NA) %>% 
  smoothPred(mod3,whichSmooth=21) %>% 
  rename(x=distMat) %>% 
  effectPlot(leg=F)+labs(x='Distance (m)',y='Urban effect')

#Plot landscape effects
fixeffPlot <- ggarrange(p1,p2,p3,p4,p5,
                        labels=letters[1:6],nrow=2,ncol=3,
                        legend='bottom',common.legend=T) 
ggsave('./figures/Pardosa_distincta_fixeff.png',fixeffPlot,width=8,height=4,scale=2)
rm(p1,p2,p3,p4,p5,raneffPlot,fixeffPlot) #Cleanup
detach(ParDisMod)


# Harvestmen -------------------------------------------------------------

#Select only harvestmen
tempArth <- arth %>% filter(arthOrder=='Opiliones') %>% group_by(BTID) %>% summarize(n=n())

#Takes way longer to run. 5-10 mins +
# OpilioMod <- runMods(tempArth,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
# save(OpilioMod,file='./data/OpilioMod.Rdata')
load('./data/OpilioMod.Rdata')

#Check models
attach(OpilioMod)
# detach(OpilioMod)

AIC(mod1,mod2,mod3,mod4)  #Very little info beyond geography here

#Model 3 - none of the fixed terms are super important
#Model 4 - noncrop land
summary(mod3); AIC(mod3)
anova(mod3)

plot(mod3,pages=1,scheme=2,rug=F,shade=T,all.terms=T)

#Check for multicollinearity
checkMC <- concurvity(mod3,full=F)$worst
#Looks OK. Some correlation between grassland and wetland
termNames <- gsub('(te|s)\\(distMat(,endDayMat)?\\)','f',gsub('\\:oRingMat2Prop\\$',':',rownames(checkMC)))
matrixplot(checkMC,termNames,mar=c(1,6,6,1))

#Smoothing term estimates
matrixplot(abs(cov2cor(sp.vcov(mod3))),c(names(mod3$sp),'scale'),numSize=1,mar=c(1,8,8,1))

plot(mod3$sp)
round(diag(sp.vcov(mod2)))

#Variance component
varcomp <- gam.vcomp(mod3) 
round(varcomp,3)

#Temporal smoother
p1 <- data.frame(day=min(datList$day):max(datList$day)) %>% 
  smoothPred(mod3,whichSmooth=1) %>% 
  ggplot(aes(x=day,y=pred))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+
  labs(x='Day of Year',y='Activity')

#Spatial smoother
p2 <- with(datList,expand.grid(E=seq(from=min(E),to=max(E),length.out=100),N=seq(from=min(N),to=max(N),length.out=100))) %>% 
  smoothPred(mod3,whichSmooth=2) %>% 
  filter(!exclude.too.far(E,N,datList$E,datList$N,0.1)) %>% 
  ggplot(aes(x=E,y=N))+geom_raster(aes(fill=pred))+
  geom_point(data=tempTrap,aes(x=easting,y=northing))+
  scale_fill_gradient(low='blue',high='red')+
  labs(x='Easting (km)',y='Northing (km)',fill='Activity')+
  theme(legend.position = c(0.85, 0.2),legend.background = element_rect(size=0.5,linetype="solid",colour="black"))+
  maptheme

raneffPlot <- ggarrange(p1,p2,labels=letters[1:2]) #Plot spatial/temporal effects
ggsave('./figures/Opiliones_raneff.png',raneffPlot,width=8,height=4,scale=2)

#Plot important landscape effects

#Days to display on plots (early,mid,late)
dispDays <- data.frame(doy=c(173,203,232)) %>% 
  mutate(date=c('June 20','July 20','August 20')) #Actually June/July 22, but close enough...

#Trap location - less in canola, more in wetland/pivot
p1 <- data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod3,type='terms',terms='trapLoc',se.fit=T)) %>% 
  rename('pred'='pred.trapLoc','se'='pred.trapLoc.1') %>% distinct() %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96) %>% 
  mutate(trapLoc=factor(trapLoc,labels=firstUpper(levels(trapLoc)))) %>% 
  ggplot(aes(x=trapLoc,y=pred))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Activity')

#Pasture
p2 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Pasture=1) %>% 
  smoothPred(mod3,whichSmooth=9:11) %>% rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot()+labs(x='Distance (m)',y='Pasture effect')

#Plot landscape effects
fixeffPlot <- ggarrange(p1,p2,
                        labels=letters[1:6],nrow=1,ncol=2,
                        legend='bottom',common.legend=T) 
ggsave('./figures/Opiliones_fixeff.png',fixeffPlot,width=8,height=4,scale=2)
rm(p1,p2,raneffPlot,fixeffPlot) #Cleanup

detach(ParDisMod)