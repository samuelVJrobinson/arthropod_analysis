#FUNCTIONAL REGRESSION ANALYSIS OF IN-FIELD INSECT ABUNDANCE USING ONION-RING LANDSCAPE DATA
# SR WINTER 2020

library(tidyverse)
theme_set(theme_classic())
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps
library(mgcv)
library(sf)
library(gstat)
library(beepr)

# Functions --------------------------------------------------------

#Function to display matrix (m) with row/column names (nam)
matrixplot <- function(m,nam=NULL,mar=NULL,showNums=T){
  if(ncol(m)!=nrow(m)) stop('Matrix not square')
  if(!is.null(nam)){
    if(ncol(m)!=length(nam)) stop('Length of name vector != matrix dimension')
  } else {
    nam <- rep('',ncol(m))
  }
  n <- ncol(m) #Matrix dimensions
  savepar <- par()
  if(is.null(mar)) mar <- c(1,5,5,1)
  par(mfrow=c(1,1),mar=mar)
  #Dark colours indicate high multicollinearity
  image(m[1:n,n:1],axes=F,col=gray.colors(12,start=0.2,end=1,rev=T))
  mtext(text=nam[n:1], side=2, line=0.3, at=seq(0,1,length.out=nrow(m)), las=1, cex=0.8)
  mtext(text=nam, side=3, line=0.3, at=seq(0,1,length.out=nrow(m)), las=2, cex=0.8)
  
  # #shadowtext trick from:
  # #https://stackoverflow.com/questions/25631216/r-plots-is-there-any-way-to-draw-border-shadow-or-buffer-around-text-labels/39002911
  # theta <-  seq(0, 2*pi, length.out=60) #Vector of angles to use
  # r <- 0.15 #Proportion of character size to extend
  # x0 <- r*strwidth('0')
  # y0 <- r*strheight('0')

  if(showNums){
    #Add row text
    for(i in 1:n){ #Column of matrix to access
      cpos <- (i-1)/(n-1) #Col position in figure
      j <- 1:n #Rows of matrix to access
      rpos <- ((j-1)/(n-1))[n:1] #Column position in figure
      # #White background text
      # for(k in 1:length(theta)){
      #   text(cpos+cos(k)*x0,rpos+sin(k)*y0,round(m[j,i],2),col='white')
      # }
      text(cpos,rpos,round(m[j,i],2),col='darkred') #Text
      # }
    }
  }
  
  suppressWarnings(par(savepar))
}

#Takes a df of arthropod counts at each site (unique to each spp), then runs 4 landscape-level models of abundance
#mod3Vars = vector of landscape cover types contained in oRingMat2
#mod4Var = name of "noncrop" cover type to use
runMods <- function(tempArth,trap,nnDistMat,oRingMat2,Kvals=NULL,mod3Vars=NULL,mod4Var=NULL){
  if(is.null(Kvals)) Kvals <- rep(list(c(10,30,5)),4) #Make list of default knot values for spatio-temporal smoothers
  
  #Only use pitfall traps from 2017
  tempTrap <- trap %>% filter(startYear==2017,grepl('PF',BTID)) %>%
    # filter(!grepl('PF',BTID)) %>% #Some ditch sites don't have cover properly digitized at further distances, but it's OK for now
    select(BLID,BTID,pass,contains('julian'),deployedhours,trapLoc,distFrom,dist,lonTrap,latTrap,lonSite:ID) %>%
    mutate(trapdays=deployedhours/24) %>% select(-deployedhours) %>% 
    left_join(rownames_to_column(data.frame(nnDistMat),'ID'),by='ID') %>% 
    left_join(tempArth,by='BTID') %>% mutate(n=ifelse(is.na(n),0,n)) %>% filter(!is.na(grass)) %>% 
    arrange(BLID,dist,pass) %>% mutate(BLID=factor(BLID)) %>% 
    mutate_at(vars(ephemeral:noncrop),function(x) ifelse(x>1500,1500,x)) %>% 
    rename('count'='n') %>% 
    #Converts 0 dist trapLoc to distFrom (pitfall traps at 0 m are "inside of" feature)
    mutate(trapLoc=factor(ifelse(dist==0 & distFrom!='control',distFrom,trapLoc))) %>% 
    #Get UTM coordinates for traps
    mutate(lonTrap2=lonTrap,latTrap2=latTrap) %>% #Duplicate columns
    st_as_sf(coords=c('lonTrap2','latTrap2'),crs=4326) %>% 
    st_transform(3403) %>% mutate(easting=st_coordinates(.)[,1],northing=st_coordinates(.)[,2]) %>% 
    mutate(easting=(easting-mean(easting))/1000,northing=(northing-mean(northing))/1000) #Center and scale coordinates to km
  
  #Arrange matrices from oRingMat2 to correspond with rows of tempTrap
  oRingMat2 <- lapply(oRingMat2,function(x){
    x %>% as.data.frame() %>% rownames_to_column('ID') %>% 
      left_join(st_drop_geometry(select(tempTrap,ID)),by='ID') %>% 
      select(-ID) %>% as.matrix() })
  
  #Distance matrix to use in functional regression
  distMat <- matrix(rep(as.numeric(gsub('d','',colnames(oRingMat2[[1]]))),each=nrow(oRingMat2[[1]])),
                    ncol=ncol(oRingMat2[[1]]))
  #Matrix of end ("julian") days
  endDayMat <- matrix(rep(tempTrap$endjulian,times=ncol(oRingMat2[[1]])),
                      ncol=ncol(oRingMat2[[1]]))
  
  #Make list containing appropriate data
  datList <- with(tempTrap,list(count=count,trapLoc=trapLoc,trapdays=trapdays,day=endjulian,
                                E=easting,N=northing,
                                distMat=as.matrix(distMat),endDayMat=as.matrix(endDayMat)))
  
  if(is.null(mod3Vars)) mod3Vars <- names(oRingMat2) #If no variables for mod3 are provided, uses all in oRingMat2
  
  #Check names of variables
  if(any(!mod3Vars %in% names(oRingMat2))){
    stop(c(paste0(mod3Vars[!mod3Vars %in% names(oRingMat2)],collapse=','),
           paste(' not found in names(oRingMat2)')))
  }
  if(any(!mod4Var %in% names(oRingMat2))){
    stop(c(paste0(mod4Var[!mod4Var %in% names(oRingMat2)],collapse=','),
           paste(' not found in names(oRingMat2)')))
  }
  
  #Assemble data list
  datList <- c(datList,oRingMat2[mod3Vars],oRingMat2[mod4Var])
  
  #"Null model" of just spatiotemporal random effects
  # mod1 <- gam(count~offset(log(trapdays))+
  #               s(endjulian,bs='ts',k=Kvals[[2]][1])+ #Thin plate spline with shrinkage
  #               s(easting,northing,bs='ts',k=Kvals[[2]][2])+
  #               ti(northing,easting,endjulian,bs='ts',k=Kvals[[2]][3]), #Spatiotemporal interaction   
  #             data=datList,family='nb',method='REML')
  mod1 <- gam(count~offset(log(trapdays))+
                s(day,k=Kvals[[2]][1])+ #Thin plate spline with shrinkage
                s(E,N,k=Kvals[[2]][2])+
                ti(N,E,day,k=Kvals[[2]][3]), #Spatiotemporal interaction
              data=datList,family='nb',select=T,method='ML')
  cat('Finished mod1. ')
  
  
  #Model of landscape effect (trap location only)
  # mod2 <- gam(count~offset(log(trapdays))+
  #               s(day,bs='ts',k=Kvals[[1]][1])+ #Thin plate spline with shrinkage
  #               s(E,N,bs='ts',k=Kvals[[1]][2])+
  #               ti(N,E,day,bs='ts',k=Kvals[[1]][3])+ #Spatiotemporal interaction   
  #               trapLoc,
  #             data=datList,family='nb',method='REML',select=T)
  mod2 <- gam(count~offset(log(trapdays))+
                s(day,k=Kvals[[1]][1])+ #Thin plate spline with shrinkage
                s(E,N,k=Kvals[[1]][2])+
                ti(N,E,day,k=Kvals[[1]][3])+ #Spatiotemporal interaction   
                trapLoc,
              data=datList,family='nb',select=T,method='ML')
  cat('Finished mod2. ')

  #Entire set of cover classes - "kitchen sink model"
  #Model with extra shrinkage takes longer to run. Maybe 3-5 mins?
  
  #Create model formula
  # mod3Formula <- "count~offset(log(trapdays))+
  #   s(day,bs='ts',k=Kvals[[3]][1])+ 
  #   s(E,N,bs='ts',k=Kvals[[3]][2])+
  #   ti(N,E,day,bs='ts',k=Kvals[[3]][3])+ 
  #   trapLoc"
  # for(i in 1:length(mod3Vars)){ #Add in specified terms (s and ti)
  #   mod3Formula <- paste0(mod3Formula,"+ s(distMat,by=",mod3Vars[i],",bs='ts')",
  #                        " + ti(distMat,endDayMat,by=",mod3Vars[i],",bs='ts')")
  # }
  mod3Formula <- "count~offset(log(trapdays))+
    s(day,k=Kvals[[3]][1])+
    s(E,N,k=Kvals[[3]][2])+
    ti(N,E,day,k=Kvals[[3]][3])+
    trapLoc"
  for(i in 1:length(mod3Vars)){ #Add in specified terms (s and ti)
    mod3Formula <- paste0(mod3Formula,"+ s(distMat,by=",mod3Vars[i],")",
                         " + ti(distMat,endDayMat,by=",mod3Vars[i],")")
  }
  
  mod3Formula <- as.formula(mod3Formula) #Convert to formula object
  
  mod3 <- gam(formula=mod3Formula,data=datList,family='nb',select=T,method='ML')
  cat('Finished mod3. ')
  
  #Model using "Noncrop" only:
  if(!is.null(mod4Var)){ #If mod4Var is defined
    #Assemble formula
    # mod4Formula <- as.formula(paste0("count~offset(log(trapdays))+s(day,bs='ts',k=Kvals[[4]][1])+ 
    # s(E,N,bs='ts',k=Kvals[[4]][2])+ti(N,E,day,bs='ts',k=Kvals[[4]][3])+ 
    # trapLoc + ","+ s(distMat,by=",mod4Var,",bs='ts')"," + ti(distMat,endDayMat,by=",mod4Var,",bs='ts')"))
    mod4Formula <- as.formula(paste0("count~offset(log(trapdays))+s(day,k=Kvals[[4]][1])+
    s(E,N,k=Kvals[[4]][2])+ti(N,E,day,k=Kvals[[4]][3])+
    trapLoc + ","+ s(distMat,by=",mod4Var,")"," + ti(distMat,endDayMat,by=",mod4Var,")"))
    
    #Fit model
    mod4 <- gam(formula=mod4Formula,data=datList,family='nb',select=T,method='ML')
    cat('Finished mod4.')
  } else { #If variable not defined
    cat('Variable for mod4 not found.')
    mod4 <- NULL
  }
  
  # #What if there are interactions between trapping location and surrounding landscape? (eg. nearby wetlands only affect abundance in canola)
  # #This requires some lengthy coding in mgcv, but can be done. Here's an example for noncrop land
  # 
  # #Set up matrix of distance values (composition values in rows with that trap location, 0 otherwise) 
  # mmat <- model.matrix(~trapLoc+0,data=tempTrap)
  # mmat <- lapply(1:5,function(x) matrix(rep(mmat[,x],ncol(oRingNoncropProp)),ncol=ncol(oRingNoncropProp)))
  # mmat <- lapply(mmat,function(x) x*oRingNoncropProp)
  # names(mmat) <- levels(tempTrap$trapLoc)
  # 
  # mod5 <- gam(count~offset(log(trapdays))+
  #               s(day,bs='ts',k=Kvals[[4]][1])+ #Thin plate spline with shrinkage
  #               s(E,N,bs='ts',k=Kvals[[4]][2])+
  #               ti(N,E,day,bs='ts',k=Kvals[[4]][3])+ #Spatiotemporal interaction   
  #               trapLoc+ 
  #               s(distMat,by=oRingNoncropProp,bs='ts')+ 
  #               s(distMat,by=mmat$canola,bs='ts')+
  #               s(distMat,by=mmat$ditch,bs='ts')+
  #               s(distMat,by=mmat$native,bs='ts')+
  #               s(distMat,by=mmat$pivot,bs='ts')+
  #               s(distMat,by=mmat$wetland,bs='ts')+
  #               ti(distMat,endDayMat,by=oRingNoncropProp,bs='ts') #Noncrop
  #             ,data=datList,family='nb')
  # summary(mod5)
  # plot(mod5,rug=F,scheme=2,pages=1,all.terms=T)
  # #Problem: this requires 5 new terms to be written for each cover category. Leaving this out for now.
  
  #List of objects to return
  retList <- list(datList=datList,tempTrap=tempTrap,mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4,Kvals=Kvals)
  return(retList)
}

# Load everything ---------------------------------------------------------
load('./data/cleanData.Rdata') #Load site, trap, arth data
load('./data/geoData.Rdata') #Load NNdist and oring data
load('./data/geoDataAAFC.Rdata') #Load oring data extracted from AAFC data

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
oRingMat2Prop$Shrubland <- NULL; oRingMat2Prop$Forest <- NULL
# #Water + Wetland
# oRingMat2Prop$WetlandWater <- oRingMat2Prop$Wetland + oRingMat2Prop$Water
# oRingMat2Prop$Wetland <- NULL; oRingMat2Prop$Water <- NULL

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

#Get all 4 models for P. melanarius
mod3Vars <- c('Grassland','Cereal','Canola','Pasture','Wetland','TreeShrub','Pulses','Urban')
mod4Var <- 'NonCrop'

debugonce(runMods)
PteMelMod <- runMods(tempArth,trap,nnDistMat,oRingMat2Prop,
                     Kvals=rep(list(c(5,70,4)),4),mod3Vars=mod3Vars,mod4Var=mod4Var)
# save(PteMelMod,file='./data/PteMelMod.Rdata')
load('./data/PteMelMod.Rdata')

#Check models
attach(PteMelMod)

AIC(mod1,mod2,mod3,mod4) #Best model has separate land cover types + SpatioTemporal effect




#Model 3 
summary(mod3); AIC(mod3)

#Check k values
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F)

#Check for multicollinearity
concurvity(mod3)
checkMC <- concurvity(saveMod,full=F)$worst
termNames <- gsub('s\\(distMat(,endDayMat)?\\)','s',gsub('\\:oRingMat2Prop\\$',':',rownames(checkMC)))
termNames <- gsub('ti\\(distMat(,endDayMat)?\\)','ti',termNames)
termNames <- gsub('s\\(endjulian\\)','s:time',gsub('s\\(easting,northing\\)','s:space',termNames))
termNames <- gsub('ti\\(northing,easting,endjulian\\)','ti:spacetime',termNames)

# c(5, 4, 4, 2) + 0.1. #Default mar
par(mar=c(1, 6, 6, 1)); matrixplot(checkMC,termNames); par(mar=c(5, 4, 4, 2) + 0.1)


coef(mod1)





#Plot spatial/temporal effects
png(file = './figures/P_melanarius_raneff.png',width=2000,height=800,pointsize=20)
par(mfrow=c(1,3))
plot(mod3,scale=0,scheme=2,resid=F,select=1,xlab='Day of year',ylab='Effect',main='Temporal random effect',shade=T,shift=coef(mod3)[1])
plot(mod3,scale=0,scheme=2,resid=F,select=2,main='Spatial random effect',cex=0.5,pch=4,shift=coef(mod3)[1])
plot(mod3,scale=0,scheme=2,resid=F,select=3,main='Spatiotemporal interaction',shift=coef(mod3)[1])
dev.off(); par(mfrow=c(1,1))

#Plot landscape effects
png(file = './figures/P_melanarius_fixef.png',width=1200,height=1200,pointsize=20)
par(mfrow=c(3,2))

plot(mod3,scale=0,scheme=2,rug=F,select=5,xlab='Distance',ylab='Day of year',main='Grassland (p=0.002)')
plot(mod3,scale=0,scheme=1,rug=F,select=8,xlab='Distance',ylab='Coefficient',main='Pasture (p=0.05)'); abline(h=0,lty='dashed',col='red')
# plot(mod3,scale=0,scheme=1,rug=F,select=12,xlab='Distance',ylab='Coefficient',main='Tree/Shrub (s) (p=0.06)'); abline(h=0,lty='dashed',col='red')
plot(mod3,scale=0,scheme=2,rug=F,select=13,xlab='Distance',ylab='Day of year',main='Tree/Shrub (ti) (p<0.001)');  
plot(mod3,scale=0,scheme=1,rug=F,select=14,xlab='Distance',ylab='Coefficient',main='Pulses (p=0.004)'); abline(h=0,lty='dashed',col='red')
plot(mod3,scale=0,scheme=1,rug=F,select=16,xlab='Distance',ylab='Coefficient',main='Urban (s) (p<0.001)'); abline(h=0,lty='dashed',col='red')
plot(mod3,scale=0,scheme=2,rug=F,select=17,xlab='Distance',ylab='Day of year',main='Urban (ti) (p<0.001)'); # 
# plot(mod3,scale=0,scheme=1,rug=F,select=4,xlab='Distance',ylab='Coefficient',main='Pasture'); abline(h=0,lty='dashed',col='red')
# plot(mod3,scale=0,scheme=1,rug=F,select=6,xlab='Distance',ylab='Coefficient',main='Pulses'); abline(h=0,lty='dashed',col='red')
dev.off(); par(mfrow=c(1,1))

#Outlier plot over time at each site
tempTrap %>% mutate(resid=resid(mod3)) %>% 
  ggplot(aes(endjulian,resid))+geom_point()+
  facet_wrap(~BLID)+geom_hline(yintercept=0)

#Summed outlier plot (like RSS but in log-space) over space
tempTrap %>% mutate(resid=resid(mod3)) %>% 
  group_by(ID) %>% summarize(sumRes=sum(abs(resid))) %>% ungroup() %>% 
  # filter(sumRes>quantile(sumRes,0.2)) %>% 
  ggplot()+
  geom_sf(aes(size=sumRes,col=sumRes),alpha=0.5,show.legend=F)+
  scale_colour_gradient(low='blue',high='red')

#Variogram of summed residuals - not much of a pattern?
tempTrap %>% mutate(resid=resid(mod3)) %>% 
  group_by(ID) %>% summarize(sumRes=sum(abs(resid))) %>% ungroup() %>% 
  as_Spatial() %>% variogram(sumRes~1,.) %>% plot()

#Look at locations with weird temporal outliers only
weirdPoints <- c(10348,11891,12000,12754,13825,14025,17882,20007,25196) %>% 
  factor(.,levels=levels(tempTrap$BLID))

tempTrap %>% mutate(resid=resid(mod3)) %>% 
  filter(BLID %in% weirdPoints) %>%
  ggplot(aes(endjulian,resid))+geom_point()+
  facet_wrap(~BLID)+geom_hline(yintercept=0)

summary(mod4)
plot(mod4,pages=1,scale=0,scheme=2)

#How do the 4 models compare?
AIC(mod1,mod2,mod3,mod4)

#Looks like the ring model of landscape does better. Important landscape features seem to be:
# Urban (spatial + temporal), Pasture, Pulses, Tree/Shrubs (weak)
detach(PteMelMod)

# Pardosa distincta (wolf spider) -----------------------------------------------------------------

#What spp of spiders are present?
arth %>% filter(grepl('PF',BTID),arthOrder=='Araneae') %>%
  mutate(genSpp=paste(genus,species,sep=' ')) %>% group_by(family,genus,species) %>% 
  summarize(n=n()) %>% arrange(family,genus,desc(n)) %>% data.frame()
#Lots of Pardosa distincta and P. moesta

# #Select only wolf spiders
# tempArth <- arth %>% filter(family=='Lycosidae') %>% group_by(BTID) %>% summarize(n=n())

#Select only Pardosa distincta
tempArth <- arth %>% filter(genus=='Pardosa',species=='distincta') %>% group_by(BTID) %>% summarize(n=n())
tempArth <- arth %>% filter(genus=='Pardosa',species=='moesta') %>% group_by(BTID) %>% summarize(n=n())

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
  mutate(easting=(easting-mean(easting))/1000,northing=(northing-mean(northing))/1000) #Center and scale coordinates to km

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

mod1 <- gam(count~offset(log(trapdays))+
              s(endjulian,bs='ts',k=10)+ #Thin plate spline with shrinkage
              s(easting,northing,bs='gp',k=10)+
              ti(northing,easting,endjulian,bs='ts')+ #Spatiotemporal interaction   
              # s(grass)+s(wetlands)+
              trapLoc,
            data=tempTrap,family='nb',method='REML')
par(mfrow=c(2,2)); gam.check(mod1); par(mfrow=c(1,1))
summary(mod1)  
plot(mod1,scheme=2,pages=1,all.terms=F,residuals=F)

#Model of landscape effect (ring composition from AAFC rasters)

#Proportion area within each ring
oRingMat2Prop <-  lapply(oRingMat2,function(x) x/Reduce('+',oRingMat2))
#Arrange matrices from oRingMat2 to correspond with rows of tempTrap
oRingMat2Prop <- lapply(oRingMat2Prop,function(x){
  x %>% as.data.frame() %>% rownames_to_column('ID') %>% 
    left_join(st_drop_geometry(select(tempTrap,ID)),by='ID') %>% 
    select(-ID) %>% as.matrix() })
#Add proportion "noncrop"
nonCropClasses <- c('Grassland','Pasture','Wetland','Urban','Shrubland','Forest','Water','Barren')
oRingNoncropProp <- Reduce('+',oRingMat2Prop[names(oRingMat2Prop) %in% nonCropClasses])
#Add proportion Trees + Shrubs
oRingMat2Prop$TreeShrub <- oRingMat2Prop$Shrubland + oRingMat2Prop$Forest
#Distance matrix to use in functional regression
distMat <- matrix(rep(as.numeric(gsub('d','',colnames(oRingMat2Prop[[1]]))),each=nrow(oRingMat2Prop[[1]])),
                  ncol=ncol(oRingMat2Prop[[1]]))
#Matrix of end ("julian") days
endDayMat <- matrix(rep(tempTrap$endjulian,times=ncol(oRingMat2Prop[[1]])),
                    ncol=ncol(oRingMat2Prop[[1]]))

#Problem with landscape matrix: mismatches between "trapLoc" and landscape composition at "0m" (cell that trap is located in)
tempTrap %>% st_drop_geometry() %>% select(trapLoc,ID) %>% cbind(.,sapply(oRingMat2Prop,function(x) x[,1])) %>% distinct() %>% 
  pivot_longer(cols=Grassland:TreeShrub) %>% filter(value!=0) %>% group_by(trapLoc,name) %>% summarize(s=sum(value)) %>% data.frame()

#Solution: change Om column to match trapLoc column. ditch -> Urban, native -> Grassland, pivot -> Grassland, canola -> Canola, wetland -> Wetland
oRingMat2Prop <- lapply(oRingMat2Prop,function(x){x[,1] <- rep(0,nrow(x)); return(x)}) #Change 0m values to 0

#Dataframe of 0m values
zeroCol <- tempTrap %>% st_drop_geometry() %>% select(trapLoc) %>% cbind(.,model.matrix(~trapLoc+0,data=tempTrap)) %>% 
  rename(Canola=trapLoccanola,Urban=trapLocditch,Grassland=trapLocnative,Grassland2=trapLocpivot,Wetland=trapLocwetland) %>% 
  mutate(Grassland=Grassland+Grassland2) %>% select(-Grassland2)

for(i in 1:length(names(zeroCol[,-1]))) oRingMat2Prop[[names(zeroCol[,-1])[i]]][,1] <- zeroCol[,-1][,i] #Replace 0m values
remove(zeroCol,i) #Cleanup



#"Null model" - just spatiotemporal component
mod2 <- gam(count~offset(log(trapdays))+
              s(endjulian,bs='ts',k=10)+ #Thin plate spline with shrinkage
              s(easting,northing,bs='gp',k=10)+
              ti(northing,easting,endjulian,bs='ts'), #Spatiotemporal interaction   
            data=tempTrap,family='nb',select=F,method='REML')
summary(mod2); AIC(mod2)
plot(mod2,pages=1,scheme=2,resid=F)
par(mfrow=c(2,2)); gam.check(mod2); par(mfrow=c(1,1))

#Landscape model
mod3b <- gam(count~offset(log(trapdays))+
              s(endjulian,bs='ts',k=10)+ #Thin plate spline with shrinkage
              s(easting,northing,bs='ts',k=10)+
              ti(northing,easting,endjulian,bs='ts')+ #Spatiotemporal interaction          
              s(distMat,by=Grassland,bs='ts')+ ti(distMat,endDayMat,by=Grassland,bs='ts')+
              s(distMat,by=Canola,bs='ts')+ ti(distMat,endDayMat,by=Canola,bs='ts')+
              s(distMat,by=Pasture,bs='ts')+ti(distMat,endDayMat,by=Pasture,bs='ts')+
              s(distMat,by=Wetland,bs='ts')+ ti(distMat,endDayMat,by=Wetland,bs='ts')+ 
              s(distMat,by=TreeShrub,bs='ts')+ ti(distMat,endDayMat,by=TreeShrub,bs='ts')+
              s(distMat,by=Pulses,bs='ts')+ ti(distMat,endDayMat,by=Pulses,bs='ts')+
              s(distMat,by=Urban,bs='ts')+ ti(distMat,endDayMat,by=Urban,bs='ts')
            ,
            data=datList,family='nb',select=F,method='REML')
summary(mod3b); AIC(mod3b)
plot(mod3,pages=1,scheme=2,rug=F,scale=0)
gam.check(mod3)

#Check for multicollinearity
checkMC <- concurvity(mod3,full=F)$worst
#Looks OK. Some correlation between grassland and wetland
termNames <- gsub('(te|s)\\(distMat(,endDayMat)?\\)','f',gsub('\\:oRingMat2Prop\\$',':',rownames(checkMC)))
matrixplot(checkMC,termNames)

#Plot spatial/temporal effects
png(file = './figures/Lycosidae_raneff.png',width=2000,height=800,pointsize=20)
par(mfrow=c(1,3))
plot(mod3,scale=0,scheme=2,resid=F,select=1,xlab='Day of year',ylab='Effect',main='Temporal random effect',shade=T,shift=coef(mod3)[1])
plot(mod3,scale=0,scheme=2,resid=F,select=2,main='Spatial random effect',cex=0.5,pch=4,shift=coef(mod3)[1])
plot(mod3,scale=0,scheme=2,resid=F,select=3,main='Spatiotemporal interaction',shift=coef(mod3)[1])
dev.off()

#Plot landscape effects
png(file = './figures/Lycosidae_fixef.png',width=1200,height=1200,pointsize=20)
par(mfrow=c(2,2))
plot(mod3,scale=0,shade=T,rug=F,select=4,xlab='Distance',ylab='Effect',main='Canola'); abline(h=0,lty='dashed',col='red')
plot(mod3,scale=0,scheme=2,rug=F,select=5,xlab='Distance',ylab='Day of year',main='Wetland')
plot(mod3,scale=0,shade=T,rug=F,select=6,xlab='Distance',ylab='Effect',main='Pulses'); abline(h=0,lty='dashed',col='red')
plot(mod3,scale=0,scheme=2,rug=F,select=7,xlab='Distance',ylab='Day of year',main='Urban')
dev.off()




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

