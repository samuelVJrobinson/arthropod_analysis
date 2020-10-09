# FUNCTIONAL REGRESSION ANALYSIS OF IN-FIELD INSECT ABUNDANCE USING ONION-RING LANDSCAPE DATA
# SR WINTER 2020

library(tidyverse)
library(ggpubr)
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

#Combine classes that are concurved/represent similar things
#Trees + Shrubs
oRingMat2Prop$TreeShrub <- oRingMat2Prop$Shrubland + oRingMat2Prop$Forest
# oRingMat2Prop$Shrubland <- NULL; oRingMat2Prop$Forest <- NULL

#Grassland + Wetland
oRingMat2Prop$GrassWetland <- oRingMat2Prop$Wetland + oRingMat2Prop$Grassland
# oRingMat2Prop$Wetland <- NULL; oRingMat2Prop$Grassland <- NULL

#Create model formulas
modFormulas <- 'count~offset(log(trapdays))+s(day,k=10,bs=basisFun)+s(E,N,k=50,bs=basisFun)' #Temporal + spatial
# modFormulas <- paste0(modFormulas,'+ti(N,E,day,k=5,bs=basisFun)') # Add spatiotemporal interaction
modFormulas[c(2,3,4)] <- paste0(modFormulas[1],'+trapLoc-1')
mod3Vars <- c('GrassWetland','Canola','Pasture','TreeShrub','Pulses','Urban') #Variables for mod3
# mod3Vars <- c('GrassWetland','Cereal','Canola','Pasture','TreeShrub','Pulses','Flax','Urban') #As above, but includes cereal, flax
# mod3Vars <- c('Grassland','Cereal','Canola','Pasture','Pulses','Wetland','Urban','Shrubland','Flax','Forest') #Expanded variables for mod3 - adding Water causes problems
#Cereal may be causing problems in estimation - collinear with canola

#Assemble formulas:
for(i in 1:length(mod3Vars)){ #Add in specified terms (s and ti)
  #Main effects + interaction (s + ti)
  modFormulas[3] <- paste0(modFormulas[3],"+ s(distMat,by=",mod3Vars[i],",bs=basisFun)")
  modFormulas[3] <- paste0(modFormulas[3],"+ s(endDayMat,by=",mod3Vars[i],",bs=basisFun)")
  modFormulas[3] <- paste0(modFormulas[3],"+ ti(distMat,endDayMat,by=",mod3Vars[i],",bs=basisFun)")
  # #Full interaction (te)
  # modFormulas[3] <- paste0(modFormulas[3],"+ te(distMat,endDayMat,by=",mod3Vars[i],",bs=basisFun)")
}
modFormulas[4] <- paste0(modFormulas[2],'+s(distMat,by=NonCrop,bs=basisFun)+s(endDayMat,by=NonCrop,bs=basisFun)+ti(distMat,endDayMat,by=NonCrop,bs=basisFun)')

#Parameters for multiplot figures
#Landscape-level parameters
landscapeFigX <- 8 #width
landscapeFigY <- 6 #height
landscapeLabel <- 20 #label size (letters a-f) 
#Spatiotemporal smoother
stFigX <- 6 #width
stFigY <- 3 #height
stLegPosX <- 0.85 #Legend position X
stLegPosY <- 0.3 #Legend position Y
theme_set(theme_classic()) #Classic theme
maptheme <- theme(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better theme for maps
smoothCols <- c('blue','purple','red') #Colours for smoothing lines

# Pterostichus melanarius ----------------------------------

#Select only P. melanarius
tempArth <- arth %>% filter(genus=='Pterostichus',species=='melanarius',year==2017) %>% group_by(BTID) %>% summarize(n=n())

# PteMelMod <- runMods(tempArth,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts',doublePenalize=FALSE); beep(1)
# save(PteMelMod,file='./data/PteMelMod.Rdata')
load('./data/PteMelMod.Rdata')

#Tried running this with double-penalization instead of shrinkage. Results were similar, but some of the shrinkage 
#Also tried running this with Cereal included, despite concurvity. Fits, but causes a bunch of weird instability in other results.

#Check models
attach(PteMelMod)

# detach(PteMelMod)
AIC(mod1,mod2,mod3,mod4) #Best model has separate land cover types + SpatioTemporal effect

#Model 3 - landscape effects matter quite a bit
summary(mod3); AIC(mod3)

# #Checking whether flax adds anything
# fitMethod <- 'REML'; fitFam <- 'nb'; basisFun <- 'ts'; doublePenalize <- F
# mod3a <- update(mod3,.~. - s(distMat, by = Flax, bs = basisFun) - s(endDayMat, by = Flax, bs = basisFun) - ti(distMat, endDayMat, by = Flax, bs = basisFun))
# summary(mod3a)
# AIC(mod3,mod3a) #AIC is better in model with flax, but this doesn't mean the effect is "real". I think it's better to remove it, since only 1 site was near flax

#Variogram of residuals - no pattern
tempTrap %>% mutate(resid=resid(mod3)) %>% 
  # group_by(ID) %>% summarize(resid=sum(abs(resid))) %>% ungroup() %>% 
  as_Spatial() %>% variogram(resid~1,.) %>% plot()

#Check k values and residual distribution
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F,seWithMean=T)

#Check for concurvity between smoothers
concurvity(mod3) #High concurvity
checkMC <- concurvity(mod3,full=F)$estimate
termNames <- rownames(checkMC)
termNames[1] <- 'Trap location'
# termNames <- gsub('s\\(distMat(,endDayMat)?\\)','s',gsub('\\:oRingMat2Prop\\$',':',rownames(checkMC)))
# termNames <- gsub('ti\\(distMat(,endDayMat)?\\)','ti',termNames)
# termNames <- gsub('s\\(endjulian\\)','s:time',gsub('s\\(easting,northing\\)','s:space',termNames))
# termNames <- gsub('ti\\(northing,easting,endjulian\\)','ti:spacetime',termNames)
termNames <- gsub('endDayMat','Time',termNames)
termNames <- gsub('distMat','Distance',termNames)
termNames <- gsub('GrassWetland','Grassland',termNames)
termNames <- gsub('TreeShrub','Woodland',termNames)

N <- 1:nrow(checkMC) #Plot all smooth terms
# N <- 4:nrow(checkMC) #Plot with only landscape smoothers
lN <- length(N)
png('./figures/coverCorPlots/concurvityEstimate_reduced.png',1200,1000,pointsize=20)
matrixplot(checkMC[N,N],mar=c(1, 10, 10, 1),termNames[N])
abline(h=seq(0-1/lN/2,1+1/lN/2,length.out=1+lN/3)) #Lines to separate terms
abline(v=seq(0-1/lN/2,1+1/lN/2,length.out=1+lN/3))
dev.off(); rm(N,lN)
#High concurvity seems to be coming from similar terms (e.g. s(dist):Canola & s(day):Canola), but not really between terms

#Problem terms:
#Cereal ~ Canola (bad), 
#Wetland ~ Grassland (bad)
#Urban ~ Grassland,Wetland,Canola (not as bad)

#Check for correlation in smoothers - seems to be OK for reduced set, but starts acting strangely if Water included
matrixplot(abs(cov2cor(sp.vcov(mod3))),c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75)

#Hessian
matrixplot(mod3$outer.info$hess,c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75)

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
  theme(legend.position = c(stLegPosX, stLegPosY),legend.background = element_rect(size=0.5,linetype="solid",colour="black"))+
  maptheme

raneffPlot <- ggarrange(p1,p2,labels=letters[1:2],font.label=list(size=landscapeLabel)) #Plot spatial/temporal effects
ggsave('./figures/Pterostichus_melanarius_raneff.png',raneffPlot,width=stFigX,height=stFigY,scale=1.2)

#Plot significant landscape effects

#Days to display on plots (early,mid,late)
dispDays <- data.frame(doy=c(173,203,232)) %>% 
  mutate(date=c('June 20','July 20','August 20')) #Actually June/July 22, but close enough...

#Get order of terms to plot
cbind(1:length(mod3$smooth),sapply(mod3$smooth,function(x) x$label))

#Trap location 
p1 <- data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod3,type='terms',terms='trapLoc',se.fit=T)) %>% 
  rename('pred'='pred.trapLoc','se'='pred.trapLoc.1') %>% distinct() %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96,trapLoc=gsub('ditch','road\nmargin',trapLoc)) %>% 
  mutate(trapLoc=gsub('native','grassland',trapLoc),trapLoc=gsub('pivot','field\nedge',trapLoc)) %>%
  mutate(trapLoc=factor(trapLoc,labels=sort(firstUpper(trapLoc)))) %>% 
  ggplot(aes(x=trapLoc,y=pred))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Activity')

#Grass/wetland effect (space + time)
p2 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),GrassWetland=1) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('GrassWetland',sapply(mod3$smooth,function(x) x$label)))) %>% 
  rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot(leg=F,cols=smoothCols)+labs(x='Distance (m)',y='Grassland effect')

#Canola effect (space + time)
p3 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Canola=1) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('Canola',sapply(mod3$smooth,function(x) x$label)))) %>% 
  rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot(leg=T,cols=smoothCols)+labs(x='Distance (m)',y='Canola effect')

#Pulse effect (time)
p4 <- data.frame(endDayMat=149:241,Pulses=1,y=NA) %>% 
  smoothPred(mod3,whichSmooth=16) %>% rename(x=endDayMat) %>% 
  mutate(x=as.Date(paste0(x,'-2017'),format='%j-%Y')) %>% 
  effectPlot(leg=F)+labs(x='Time of year',y='Pulse effect')

#Plot landscape effects
fixeffPlot <- ggarrange(p1,p2,p3,p4,
                        labels=letters[1:5],nrow=2,ncol=2,
                        legend='bottom',common.legend=T,font.label=list(size=landscapeLabel)) 
ggsave('./figures/Pterostichus_melanarius_fixeff.png',fixeffPlot,width=landscapeFigX,height=landscapeFigY,scale=1)
rm(p1,p2,p3,p4,raneffPlot,fixeffPlot) #Cleanup

#Looks like the ring model of landscape does better. Important landscape features seem to be:
# Grassland/Wetland, Canola, Pulses

#Trying out some distance plots.
#Need to identify sites at differing distances into fields

#Figure out which sites have the higest amount of pasture
PteMelMod$datList$Pasture %>% apply(.,1,sum) %>%
  data.frame(pasture=.,ID=PteMelMod$tempTrap$ID)

#Choose only 26133 at pass 4
chooseMe <- which(with(PteMelMod$tempTrap,BLID=='26133'&pass==4))

dat <- lapply(PteMelMod$datList,function(x) if(is.matrix(x)) x[chooseMe,] else x[chooseMe])

PteMelMod$tempTrap$BTID

detach(PteMelMod)

# Pardosa distincta (wolf spider) -----------------------------------------------------------------

#Select only Pardosa distincta
tempArth <- arth %>% filter(genus=='Pardosa',species=='distincta',year==2017) %>% group_by(BTID) %>% summarize(n=n())
# tempArthF <- arth %>% filter(genus=='Pardosa',species=='distincta',year==2017,sex=='F') %>% group_by(BTID) %>% summarize(n=n())
# tempArthM <- arth %>% filter(genus=='Pardosa',species=='distincta',year==2017,sex=='M') %>% group_by(BTID) %>% summarize(n=n())

# Takes way longer to run. 5-10 mins +
# ParDisMod <- runMods(tempArth,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1) #Full model
# ParDisModF <- runMods(tempArthF,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1) #Female model
# ParDisModM <- runMods(tempArthM,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1) #Male model

# save(ParDisMod,file='./data/ParDisMod.Rdata')
load('./data/ParDisMod.Rdata')

#Check models
attach(ParDisMod)
AIC(mod1,mod2,mod3,mod4) 

#Model 3 
summary(mod3); AIC(mod3)
plot(mod3,pages=1,scheme=2,rug=F,shade=T,all.terms=T)

#Check for multicollinearity (same as above)

#Check residuals
sumres <- with(datList,data.frame(trapLoc,day,E,N)) %>% 
  mutate(EN=factor(paste(E,N,sep='_'))) %>% 
  mutate(res=residuals(mod3)) %>% group_by(EN) %>% 
  summarize(E=first(E),N=first(N),res=sum(res))

ggplot(sumres,aes(x=E,y=N,size=res,col=res))+geom_point(alpha=0.5)+
  scale_colour_gradient(low='red',high='blue')

gam(res~s(E,N),data=sumres) %>% plot(.,scheme=2,rug=F)

#Check k values and residual distribution
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F,seWithMean=T)

#Smoothing term estimates - problem here with distance term for TreeShrub, but this appears to be caused by the large number of terms. When "Urban" is dropped (weak effects), smoothing estimate for TreeShrub is fine, and no other effects change.
matrixplot(abs(cov2cor(sp.vcov(mod3))),c('scale',names(mod3$sp)),numSize=1,mar=c(1,8,8,1))

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
  theme(legend.position = c(stLegPosX, stLegPosY),legend.background = element_rect(size=0.5,linetype="solid",colour="black"))+
  maptheme

raneffPlot <- ggarrange(p1,p2,labels=letters[1:2]) #Plot spatial/temporal effects
ggsave('./figures/Pardosa_distincta_raneff.png',raneffPlot,width=stFigX,height=stFigY,scale=1.2)

#Plot important landscape effects

#Get order of terms to plot
data.frame(round(anova(mod3)$s.table,3)) %>% rownames_to_column('smoother')

#Days to display on plots (early,mid,late)
dispDays <- data.frame(doy=c(173,203,232)) %>% 
  mutate(date=c('June 20','July 20','August 20')) #Actually June/July 22, but close enough...

#Trap location - far less in canola
p1 <- data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod3,type='terms',terms='trapLoc',se.fit=T)) %>% 
  rename('pred'='pred.trapLoc','se'='pred.trapLoc.1') %>% distinct() %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96,trapLoc=gsub('ditch','road\nmargin',trapLoc)) %>% 
  mutate(trapLoc=gsub('native','grassland',trapLoc),trapLoc=gsub('pivot','field\nedge',trapLoc)) %>%
  mutate(trapLoc=factor(trapLoc,labels=sort(firstUpper(trapLoc)))) %>% 
  ggplot(aes(x=trapLoc,y=pred))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Activity')

# #Canola (space + time: weak)
# p2 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Canola=1) %>% 
#   smoothPred(mod3,whichSmooth=which(grepl('Canola',sapply(mod3$smooth,function(x) x$label)))) %>% 
#   rename(x=distMat,y=endDayMat) %>% 
#   mutate(y=factor(y,labels=dispDays$date)) %>% 
#   effectPlot(leg=T,cols=smoothCols)+labs(x='Distance (m)',y='Canola effect')

#Pasture (space)
p3 <- expand.grid(distMat=seq(30,1500,30),Pasture=1,y=NA) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('Pasture',sapply(mod3$smooth,function(x) x$label)))[1]) %>% 
  rename(x=distMat) %>% 
  effectPlot(leg=F)+labs(x='Distance (m)',y='Pasture effect')

#Tree/shrub (time)
p4 <- data.frame(endDayMat=149:241,TreeShrub=1,y=NA) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('TreeShrub',sapply(mod3$smooth,function(x) x$label)))[2]) %>% 
  rename(x=endDayMat) %>% 
  mutate(x=as.Date(paste0(x,'-2017'),format='%j-%Y')) %>% 
  effectPlot(leg=F)+labs(x='Time of year',y='Woodland effect')

#Plot landscape effects
fixeffPlot <- ggarrange(p1,p3,p4,
                        labels=letters[1:3],nrow=2,ncol=2,
                        legend='bottom',common.legend=T,font.label=list(size=landscapeLabel)) 
ggsave('./figures/Pardosa_distincta_fixeff.png',fixeffPlot,width=landscapeFigX,height=landscapeFigY,scale=1)
rm(p1,p2,p3,p4,raneffPlot,fixeffPlot) #Cleanup
detach(ParDisMod)

# Pardosa moesta (wolf spider) --------------------------------------------

#Select only Pardosa moesta
tempArth <- arth %>% filter(genus=='Pardosa',species=='moesta',year==2017) %>% group_by(BTID) %>% summarize(n=n())

# Takes way longer to run. 5-10 mins +
# ParMoeMod <- runMods(tempArth,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
# save(ParMoeMod,file='./data/ParMoeMod.Rdata')
load('./data/ParMoeMod.Rdata')

#Check models
attach(ParMoeMod)

AIC(mod1,mod2,mod3,mod4) 

#Model 3 
summary(mod3); AIC(mod3)
plot(mod3,pages=1,scheme=2,rug=F,shade=T,all.terms=T)

#Check for multicollinearity (same as above)

#Check k values and residual distribution
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F,seWithMean=T)

#Smoothing term estimates
matrixplot(abs(cov2cor(sp.vcov(mod3))),c('scale',names(mod3$sp)),numSize=1,mar=c(1,9,9,1))

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
  theme(legend.position = c(stLegPosX, stLegPosY),legend.background = element_rect(size=0.5,linetype="solid",colour="black"))+
  maptheme

raneffPlot <- ggarrange(p1,p2,labels=letters[1:2]) #Plot spatial/temporal effects
ggsave('./figures/Pardosa_moesta_raneff.png',raneffPlot,width=stFigX,height=stFigY,scale=1.2)

#Plot important landscape effects

#Get order of terms to plot
data.frame(round(anova(mod3)$s.table,3)) %>% rownames_to_column('smoother')

#Days to display on plots (early,mid,late)
dispDays <- data.frame(doy=c(173,203,232)) %>% 
  mutate(date=c('June 20','July 20','August 20')) #Actually June/July 22, but close enough...

#Trap location - far less in canola
p1 <- data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod3,type='terms',terms='trapLoc',se.fit=T)) %>% 
  rename('pred'='pred.trapLoc','se'='pred.trapLoc.1') %>% distinct() %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96,trapLoc=gsub('ditch','road\nmargin',trapLoc)) %>% 
  mutate(trapLoc=gsub('native','grassland',trapLoc),trapLoc=gsub('pivot','field\nedge',trapLoc)) %>%
  mutate(trapLoc=factor(trapLoc,labels=sort(firstUpper(trapLoc)))) %>% 
  ggplot(aes(x=trapLoc,y=pred))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Activity')

#Grass/Wetland (space)
p2 <- expand.grid(distMat=seq(30,1500,30),GrassWetland=1,y=NA) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('GrassWetland',sapply(mod3$smooth,function(x) x$label)))[1]) %>% 
  rename(x=distMat) %>% 
  effectPlot(leg=F)+labs(x='Distance (m)',y='Grassland effect')

#Canola (space + time)
p3 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Canola=1) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('Canola',sapply(mod3$smooth,function(x) x$label)))) %>%
  rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot(cols=smoothCols)+labs(x='Distance (m)',y='Canola effect')

#Pulses (space)
p4 <- expand.grid(distMat=seq(30,1500,30),Pulses=1,y=NA) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('Pulses',sapply(mod3$smooth,function(x) x$label)))[1]) %>%
  rename(x=distMat) %>% 
  effectPlot(leg=F)+labs(x='Distance (m)',y='Pulse effect')

#Urban effect (space + time)
p5 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Urban=1) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('Urban',sapply(mod3$smooth,function(x) x$label)))) %>%
  rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot(cols=smoothCols)+labs(x='Distance (m)',y='Road margin effect')

#Plot landscape effects
fixeffPlot <- ggarrange(p1,p2,p3,p4,p5,
                        labels=letters[1:6],nrow=2,ncol=3,
                        legend='bottom',common.legend=T,font.label=list(size=landscapeLabel)) 
ggsave('./figures/Pardosa_moesta_fixeff.png',fixeffPlot,width=landscapeFigX,height=landscapeFigY,scale=1)
rm(p1,p2,p3,p4,p5,raneffPlot,fixeffPlot) #Cleanup
detach(ParMoeMod)

# Harvestmen -------------------------------------------------------------

#Select only harvestmen
tempArth <- arth %>% filter(arthOrder=='Opiliones',year==2017) %>% group_by(BTID) %>% summarize(n=n())

#Takes way longer to run. 5-10 mins +
# OpilioMod <- runMods(tempArth,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
# save(OpilioMod,file='./data/OpilioMod.Rdata')
load('./data/OpilioMod.Rdata')

#Check models
attach(OpilioMod)
# detach(OpilioMod)

AIC(mod1,mod2,mod3,mod4) #Much less info from landscape level

#Model 3 - none of the fixed terms are super important
#Model 4 - noncrop land
summary(mod3); AIC(mod3)
anova(mod3)
summary(mod4)

plot(mod3,pages=1,scheme=2,rug=F,shade=T,all.terms=T)

#Check for multicollinearity - same as above

#Check k values and residual distribution
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F,seWithMean=T)

#Smoothing term estimates
matrixplot(abs(cov2cor(sp.vcov(mod3))),c('scale',names(mod3$sp)),numSize=1,mar=c(1,8,8,1))

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
  theme(legend.position = c(stLegPosX, stLegPosY), legend.background = element_rect(size=0.5,linetype="solid",colour="black"))+
  maptheme

raneffPlot <- ggarrange(p1,p2,labels=letters[1:2]) #Plot spatial/temporal effects
ggsave('./figures/Opiliones_raneff.png',raneffPlot,width=stFigX,height=stFigY,scale=1.2)

#Plot important landscape effects

#Get order of terms to plot
data.frame(round(anova(mod3)$s.table,3)) %>% rownames_to_column('smoother')

#Days to display on plots (early,mid,late)
dispDays <- data.frame(doy=c(173,203,232)) %>% 
  mutate(date=c('June 20','July 20','August 20')) #Actually June/July 22, but close enough...

#Trap location - less in canola, more in wetland/pivot
p1 <- data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod3,type='terms',terms='trapLoc',se.fit=T)) %>% 
  rename('pred'='pred.trapLoc','se'='pred.trapLoc.1') %>% distinct() %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96,trapLoc=gsub('ditch','road\nmargin',trapLoc)) %>% 
  mutate(trapLoc=gsub('native','grassland',trapLoc),trapLoc=gsub('pivot','field\nedge',trapLoc)) %>%
  mutate(trapLoc=factor(trapLoc,labels=sort(firstUpper(trapLoc)))) %>% 
  ggplot(aes(x=trapLoc,y=pred))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Activity')

#Grassland/wetland (space + time)
p2 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),GrassWetland=1) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('GrassWetland',sapply(mod3$smooth,function(x) x$label)))) %>%
  rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot(cols=smoothCols)+labs(x='Distance (m)',y='Grassland effect')

#Tree/Shrub (space + time)
p3 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),TreeShrub=1) %>% 
  smoothPred(mod3,whichSmooth=which(grepl('TreeShrub',sapply(mod3$smooth,function(x) x$label)))) %>%
  rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot(cols=smoothCols)+labs(x='Distance (m)',y='Woodland effect')

#Plot landscape effects
fixeffPlot <- ggarrange(p1,p2,p3,
                        labels=letters[1:4],nrow=2,ncol=2,
                        legend='bottom',common.legend=T) 
ggsave('./figures/Opiliones_fixeff.png',fixeffPlot,width=landscapeFigX,height=landscapeFigY,scale=1)
rm(p1,p2,p3,raneffPlot,fixeffPlot) #Cleanup


#Make plots for mod4 (non-crop only)

#Temporal smoother
p1 <- data.frame(day=min(datList$day):max(datList$day)) %>% 
  smoothPred(mod4,whichSmooth=1) %>% 
  ggplot(aes(x=day,y=pred))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+
  labs(x='Day of Year',y='Activity')

#Spatial smoother
p2 <- with(datList,expand.grid(E=seq(from=min(E),to=max(E),length.out=100),N=seq(from=min(N),to=max(N),length.out=100))) %>% 
  smoothPred(mod4,whichSmooth=2) %>% 
  filter(!exclude.too.far(E,N,datList$E,datList$N,0.1)) %>% 
  ggplot(aes(x=E,y=N))+geom_raster(aes(fill=pred))+
  geom_point(data=tempTrap,aes(x=easting,y=northing))+
  scale_fill_gradient(low='blue',high='red')+
  labs(x='Easting (km)',y='Northing (km)',fill='Activity')+
  theme(legend.position = c(stLegPosX, stLegPosY), legend.background = element_rect(size=0.5,linetype="solid",colour="black"))+
  maptheme

#Trap location - less in canola, more in wetland/pivot
p3 <- data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod4,type='terms',terms='trapLoc',se.fit=T)) %>% 
  rename('pred'='pred.trapLoc','se'='pred.trapLoc.1') %>% distinct() %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96) %>% 
  mutate(trapLoc=factor(trapLoc,labels=firstUpper(levels(trapLoc)))) %>% 
  ggplot(aes(x=trapLoc,y=pred))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
  labs(x='Trap location',y='Activity')

#Noncrop
p4 <- expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),NonCrop=1) %>% 
  smoothPred(mod4,whichSmooth=3:5) %>% rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  effectPlot()+labs(x='Distance (m)',y='Noncrop effect')

ggarrange(p1,p2,p3,p4,labels=letters[1:4],nrow=2,ncol=2,legend='right',common.legend=F) 

detach(OpilioMod)
