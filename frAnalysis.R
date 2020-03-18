#FUNCTIONAL REGRESSION ANALYSIS OF IN-FIELD INSECT ABUNDANCE USING ONION-RING LANDSCAPE DATA
# SR WINTER 2020

library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
library(mgcv)
library(sf)
library(gstat)
library(beepr)

# Helper functions --------------------------------------------------------

#Function to display matrix (m) with row/column names (nam)
matrixplot <- function(m,nam){
  #Dark colours indicate high multicollinearity
  image(m,axes=F,col=gray.colors(12,start=0.2,end=1,rev=T))
  mtext(text=nam, side=2, line=0.3, at=seq(0,1,length.out=nrow(m)), las=1, cex=0.8)
  mtext(text=nam, side=3, line=0.3, at=seq(0,1,length.out=nrow(m)), las=2, cex=0.8)
  
  # #shadowtext trick from:
  # #https://stackoverflow.com/questions/25631216/r-plots-is-there-any-way-to-draw-border-shadow-or-buffer-around-text-labels/39002911
  # theta <-  seq(0, 2*pi, length.out=60) #Vector of angles to use
  # r <- 0.15 #Proportion of character size to extend
  # x0 <- r*strwidth('0')
  # y0 <- r*strheight('0')

  for(i in 1:nrow(m)){ #Add row text
    (rpos <- (i-1)/(nrow(m)-1)) #Row position
    # for(j in 1:ncol(m)){
    j <- 1:ncol(m)
    (cpos <- (j-1)/(ncol(m)-1)) #Column position
    # #White background text
    # for(k in 1:length(theta)){
    #   text(cpos+cos(k)*x0,rpos+sin(k)*y0,round(m[j,i],2),col='white')
    # }
    #Actual text
    text(cpos,rpos,round(m[j,i],2),col='darkred')
    # }
  }
}

# Load everything ---------------------------------------------------------
load('./data/cleanData.Rdata') #Load site, trap, arth data
load('./data/geoData.Rdata') #Load NNdist and oring data
load('./data/geoDataAAFC.Rdata') #Load oring data extracted from AAFC data

# Pterostichus ----------------------------------

#What spp of beetles are present?
arth %>% filter(grepl('PF',BTID),arthOrder=='Coleoptera') %>%
  mutate(genSpp=paste(genus,species,sep=' ')) %>% group_by(family,genus,species) %>% 
  summarize(n=n()) %>% arrange(genus,desc(n)) %>% data.frame()

arth %>% filter(grepl('PF',BTID),arthOrder=='Coleoptera') %>%
  group_by(genus) %>% summarize(n=n()) %>% data.frame() %>% 
  arrange(desc(n))

#Select only P. melanarius
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

#Model of landscape effect (trap location only)
#Appears that distance to grass has a positive effect, and distance to noncrop a negative effect 
#Spatio-temporal tensor doesn't fit as well as straightforward random effects. Investigate this later.
#When using trap location, canola clearly has the most 
#Is this a "cultural" species?
mod1 <- gam(count~offset(log(trapdays))+s(endjulian)+
              trapLoc+
              # s(BLID,bs='re')
              s(BLID,bs='re')
              # s(grass)+s(wetlands)+
            ,
            data=tempTrap,family='nb',method='ML')
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

p2 <- tempTrap %>% #Trap table
  select(BLID,lonSite,latSite) %>% st_drop_geometry() %>% #Get BLID and location data
  distinct() %>% #Get only distinct values
  mutate(ranef=coef(mod1)[grepl('BLID',names(coef(mod1)))]) %>% #Get random intercepts from mod1
  st_as_sf(coords=c('lonSite','latSite'),crs=4326) %>% #Convert df to spatial object
  ggplot(aes(col=ranef))+ #Plot intercepts, using intercept as colour
  geom_sf(aes(size=abs(ranef)*0.4),alpha=0.5,show.legend=T)+
  scale_colour_gradient(low='blue',high='red')+
  scale_size(guide='none')+
  labs(title='P. melanarius random intercepts')
ggsave('./figures/01_P_melanarius_intercepts.png',p2,height=6,width=8)


#Model of landscape effect (ring composition from shapefiles)

#Arrange oRing cover matrix to correspond with correct rows
tempORing <- lapply(oRingMat,function(x) x[match(tempTrap$ID,rownames(x)),])
#Total area of each ring in matrix form
tempRingArea <- matrix(rep((pi*seq(20,1000,20)^2)-(pi*(seq(20,1000,20)-20)^2),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
#Matrix of distance values
distMat <- matrix(rep(as.numeric(gsub('d','',colnames(tempORing[[1]]))),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
grassMat <- tempORing$grass/10000 #hectares grass cover within each ring
noncropMat <- tempORing$noncrop/10000 #Noncrop cover (very similar to grass)
grassMatPerc <- grassMat/tempRingArea #prop cover within each ring
noncropMatPerc <- noncropMat/tempRingArea #prop cover within each ring

mod2 <- gam(count~offset(log(trapdays))+s(endjulian)+
              s(distMat,by=grassMat)+
              s(BLID,bs='re'),
            data=tempTrap,family='nb',method='ML')

#Model of landscape effect (ring composition from AAFC rasters)
#Proportion area within each ring
oRingMat2Prop <-  lapply(oRingMat2,function(x) x/Reduce('+',oRingMat2))

#Arrange matrices from oRingMat2 to correspond with rows of tempTrap
oRingMat2Prop <- lapply(oRingMat2Prop,function(x){
  x %>% as.data.frame() %>% rownames_to_column('ID') %>% 
    left_join(st_drop_geometry(select(tempTrap,ID)),by='ID') %>% 
    select(-ID) %>% as.matrix() })

#Simplify to proportion noncrop
nonCropClasses <- c('Grassland','Pasture','Wetland','Urban','Shrubland','Forest','Water','Barren')
oRingNoncropProp <- Reduce('+',oRingMat2Prop[names(oRingMat2Prop) %in% nonCropClasses])

#Selects top 10 cover classes from full oRing set
oRingMat2 <- oRingMat2[1:10]
oRingMat2Prop <- oRingMat2Prop[1:10]

#Distance matrix to use in functional regression
distMat <- matrix(rep(as.numeric(gsub('d','',colnames(oRingMat2Prop[[1]]))),
    each=nrow(oRingMat2Prop[[1]])),ncol=ncol(oRingMat2Prop[[1]]))

#Entire set of cover classes - "kitchen sink model"
#Looks like this model doesn't do any better than earlier models

# r <- seq(1,20,1)
# mod3reml <- sapply(r,function(x){
mod3 <- gam(count~offset(log(trapdays))+
              s(endjulian,bs='gp',k=20,m=c(2,6))+
              # trapLoc+
              s(distMat,by=oRingMat2Prop$Grassland)+
              s(distMat,by=oRingMat2Prop$Canola)+
              # s(distMat,by=oRingMat2Prop$Cereal)+
              s(distMat,by=oRingMat2Prop$Pasture)+
              s(distMat,by=oRingMat2Prop$Wetland)+
              # s(distMat,by=oRingMat2Prop$Pulses)+
              # s(distMat,by=oRingMat2Prop$Urban)+
              # s(distMat,by=oRingMat2Prop$Shrubland)+
              # s(distMat,by=oRingMat2Prop$Flax)+
              # s(distMat,by=oRingMat2Prop$Forest)+
              # s(easting,northing,k=40), #Spline basis function
              s(easting,northing,bs='gp',k=30,m=c(2,0.1)), #Gaussian process basis function (Matern correlation)
            # te(easting,northing,bs='gp',m=c(3)), #Gaussian process basis function
            # s(BLID,bs='re'),
            data=tempTrap,family='nb')
#     return(mod3$gcv)
# })
# beep(1)
# plot(r,mod3reml)
gam.check(mod3)
plot(mod3,pages=1,scale=0,scheme=2,resid=F)
plot(mod3,select=6,scheme=2,cex=3)
summary(mod3)

weirdPoints <- c(11891,12000,12754,13825,14025,17882,20007,25196) %>% 
  factor(.,levels=levels(tempTrap$BLID))

tempTrap %>% mutate(resid=resid(mod3)) %>% 
  # filter(BLID %in% weirdPoints) %>% 
  ggplot(aes(endjulian,resid))+geom_point()+
  facet_wrap(~BLID)+geom_hline(yintercept=0)

#Check for multicollinearity
checkMC <- concurvity(mod3,full=F)$estimate
#Looks OK. Come correlation between grassland and wetland
matrixplot(checkMC,gsub('s\\(distMat\\)\\:oRingMat2Prop\\$','f:',rownames(checkMC)))




#Seems like some of these categories are important, but possible issues with collinearity?

#Noncrop only:
mod4 <- gam(count~offset(log(trapdays))+s(endjulian)+
              trapLoc+ #Trap location
              s(distMat,by=oRingNoncropProp)+ #Proportion noncrop
              s(BLID,bs='re'), #Site random intercept
            data=tempTrap,family='nb',method='ML')
summary(mod4)
plot(mod4,pages=1,scale=0)

#How do the 4 models compare?
AIC(mod1,mod2,mod3,mod4)

#Looks like distance is the only important factor in this
AIC(mod1,mod2,update(mod2,.~.-s(distMat,by=grassMat)+s(distMat,by=grassMatPerc)),update(mod2,.~.+s(grass)))
summary(mod2)
plot(mod2,scheme=2,pages=1)
gam.check(mod2)

plot(mod1$model$count,fitted(mod1),ylab='Predicted',xlab='Actual',pch=19,cex=0.5)
points(mod2$model$count,fitted(mod2),pch=19,cex=0.5,col='red')
abline(0,1,lty='dashed')
legend('bottomright',c('Distance','Ring comp'),fill=c('black','red'))

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

plot(mod1,scale=0)




#Summary for Pterostichus melanarius: appears to be more individuals at the centre of the fields

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
distMat <- matrix(rep(as.numeric(gsub('d','',colnames(tempORing[[1]]))),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
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
distMat <- matrix(rep(as.numeric(gsub('d','',colnames(tempORing[[1]]))),each=nrow(tempORing[[1]])),ncol=ncol(tempORing[[1]]))
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

