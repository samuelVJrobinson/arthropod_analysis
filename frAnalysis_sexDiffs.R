#Script to examine how different sexes of arthropods use the landscape differently. Marc et al 1999 note that male/female spiders have different life histories, and possibly have different habitat requirements. How does this translate into landscape patterns?

# Pterostichus melanarius ------------------------------------

#Sex not determined
arth %>% filter(genus=='Pterostichus',species=='melanarius',year==2017) %>% group_by(sex) %>% summarize(n=n())

# Pardosa distincta ------------------------

arth %>% filter(genus=='Pardosa',species=='distincta',year==2017) %>% group_by(sex) %>% summarize(n=n())

#Do males and females have different responses to landscape?
tempArthF <- arth %>% filter(genus=='Pardosa',species=='distincta',year==2017,sex=='F') %>% group_by(BTID) %>% summarize(n=n())
tempArthM <- arth %>% filter(genus=='Pardosa',species=='distincta',year==2017,sex=='M') %>% group_by(BTID) %>% summarize(n=n())

# #Tried running with/without cereal. Very little difference, so dropping cereal terms (also, concurved).
# ParDisModF <- runMods(tempArthF,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
# ParDisModM <- runMods(tempArthM,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
# save(ParDisModF, file='./data/ParDisModF.Rdata')
# save(ParDisModM, file='./data/ParDisModM.Rdata')

load('./data/ParDisModF.Rdata')
load('./data/ParDisModM.Rdata')

# Females 
attach(ParDisModF)

AIC(mod1,mod2,mod3,mod4) #Landscape model performs better
#Model 3 
summary(mod3) #Still relatively low explained deviance. 
matrixplot(abs(cov2cor(sp.vcov(mod3))),c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75) #ti(TreeShrub) has NA entries
matrixplot(sp.vcov(mod3),c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75) #ti(TreeShrub) has NA entries

#Variogram of residuals - no pattern
tempTrap %>% mutate(resid=resid(mod3)) %>% 
  # group_by(ID) %>% summarize(resid=sum(abs(resid))) %>% ungroup() %>% 
  as_Spatial() %>% variogram(resid~1,.) %>% plot()

#Check k values
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F,seWithMean=T)

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

detach(ParDisModF)

# Males
attach(ParDisModM)
AIC(mod1,mod2,mod3,mod4) #Landscape model much better

summary(mod3); AIC(mod3)
matrixplot(abs(cov2cor(sp.vcov(mod3))),c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75) #ti(Flax) has NA entries
matrixplot(sp.vcov(mod3),c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75) 

#Variogram of residuals - no pattern
tempTrap %>% mutate(resid=resid(mod3)) %>% 
  # group_by(ID) %>% summarize(resid=sum(abs(resid))) %>% ungroup() %>% 
  as_Spatial() %>% variogram(resid~1,.) %>% plot()

#Check k values
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F,seWithMean=T)

#Check concurvity
concurvity(mod3) #High concurvity
checkMC <- with(ParMoeModM2,concurvity(mod3,full=F)$estimate)
termNames <- rownames(checkMC)
termNames[1] <- 'Linear terms'
N <- 1:nrow(checkMC) #Plot all smooth terms
# N <- 4:nrow(checkMC) #Plot with only landscape smoothers
lN <- length(N)
matrixplot(checkMC[N,N],mar=c(1, 10, 10, 1),termNames[N])
abline(h=seq(0-1/lN/2,1+1/lN/2,length.out=1+lN/3)) #Lines to separate terms
abline(v=seq(0-1/lN/2,1+1/lN/2,length.out=1+lN/3))

detach(ParDisModM)

compareFM(ParDisModF$mod3,ParDisModM$mod3)

#Plot important landscape effects

#Get order of terms to plot
with(ParDisModF,data.frame(round(anova(mod3)$s.table,3))) %>% rownames_to_column('smoother')
with(ParDisModM,data.frame(round(anova(mod3)$s.table,3))) %>% rownames_to_column('smoother')

#Days to display on plots (early,mid,late)
dispDays <- data.frame(doy=c(173,203,232)) %>% 
  mutate(date=c('June 20','July 20','August 20')) #Actually June/July 22, but close enough...

#Intercepts

#Labels for x-axis
xlabs <- firstUpper(levels(ParDisModF$tempTrap$trapLoc))

#Trap location 
p1 <- with(ParDisModF,data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod3,type='terms',terms='trapLoc',se.fit=T),sex='F')) %>% 
  bind_rows(with(ParDisModM,data.frame(trapLoc=tempTrap$trapLoc,pred=predict(mod3,type='terms',terms='trapLoc',se.fit=T),sex='M'))) %>% 
  rename('pred'='pred.trapLoc','se'='pred.trapLoc.1') %>% distinct() %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96) %>% 
  mutate(trapLoc=as.character(trapLoc)) %>% 
  mutate(trapLoc=factor(trapLoc,labels=xlabs)) %>%
  ggplot(aes(x=trapLoc,y=pred,col=sex))+geom_pointrange(aes(ymax=upr,ymin=lwr),position=position_dodge(width=0.3))+
  labs(x='Trap location',y='Activity Density',col='Sex')+scale_colour_manual(values=c('red','blue'))

#GrassWetland (day of year)
p2 <- smoothPred(data.frame(endDayMat=149:241,GrassWetland=1,y=NA,sex='F'),ParDisModF$mod3,whichSmooth=4) %>% 
  bind_rows(smoothPred(data.frame(endDayMat=149:241,GrassWetland=1,y=NA,sex='M'),ParDisModM$mod3,whichSmooth=4)) %>% 
  rename(x=endDayMat) %>% 
  mutate(x=as.Date(paste0(x,'-2017'),format='%j-%Y')) %>% 
  ggplot(aes(x=x,y=pred,fill=sex))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(col=sex))+geom_hline(yintercept=0,linetype='dashed')+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  labs(y='Grass/Wetland Effect',fill='Sex',col='Sex',x='Day of Year')

#Canola (distance, day, interaction)
p3 <- smoothPred(expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Canola=1,sex='F'),ParDisModF$mod3,whichSmooth=6:8) %>% 
  bind_rows(smoothPred(expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Canola=1,sex='M'),ParDisModM$mod3,whichSmooth=6:8)) %>% 
  rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  ggplot(aes(x=x,y=pred,fill=sex,linetype=y))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.2)+
  geom_line(aes(col=sex))+geom_hline(yintercept=0,linetype='dashed')+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  labs(y='Canola Effect',fill='Sex',col='Sex',x='Distance (m)',linetype='Date')

#Pasture (distance, day, interaction)
p4 <- smoothPred(expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Pasture=1,sex='F'),ParDisModF$mod3,whichSmooth=9:11) %>% 
  bind_rows(smoothPred(expand.grid(endDayMat=c(173,203,232),distMat=seq(30,1500,30),Pasture=1,sex='M'),ParDisModM$mod3,whichSmooth=9:11)) %>% 
  rename(x=distMat,y=endDayMat) %>% 
  mutate(y=factor(y,labels=dispDays$date)) %>% 
  ggplot(aes(x=x,y=pred,fill=sex,linetype=y))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.2)+
  geom_line(aes(col=sex))+geom_hline(yintercept=0,linetype='dashed')+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  labs(y='Pasture Effect',fill='Sex',col='Sex',x='Distance (m)',linetype='Date')

#TreeShrub (day of year)
p5 <- smoothPred(data.frame(endDayMat=149:241,TreeShrub=1,y=NA,sex='F'),ParDisModF$mod3,whichSmooth=13) %>% 
  bind_rows(smoothPred(data.frame(endDayMat=149:241,TreeShrub=1,y=NA,sex='M'),ParDisModM$mod3,whichSmooth=13)) %>% 
  rename(x=endDayMat) %>% 
  mutate(x=as.Date(paste0(x,'-2017'),format='%j-%Y')) %>% 
  ggplot(aes(x=x,y=pred,fill=sex))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(col=sex))+geom_hline(yintercept=0,linetype='dashed')+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  labs(y='Tree/Shrub Effect',fill='Sex',col='Sex',x='Day of Year')

#Flax (distance)
p6 <- smoothPred(expand.grid(distMat=seq(30,1500,30),Flax=1,y=NA,sex='F'),ParDisModF$mod3,whichSmooth=18) %>% 
  bind_rows(smoothPred(expand.grid(distMat=seq(30,1500,30),Flax=1,y=NA,sex='M'),ParDisModM$mod3,whichSmooth=18)) %>% 
  rename(x=distMat) %>%
  ggplot(aes(x=x,y=pred,fill=sex))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.2)+
  geom_line(aes(col=sex))+geom_hline(yintercept=0,linetype='dashed')+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  labs(y='Flax Effect',fill='Sex',col='Sex',x='Distance (m)')

#Urban (distance)
p7 <- smoothPred(expand.grid(distMat=seq(30,1500,30),Urban=1,y=NA,sex='F'),ParDisModF$mod3,whichSmooth=21) %>% 
  bind_rows(smoothPred(expand.grid(distMat=seq(30,1500,30),Urban=1,y=NA,sex='M'),ParDisModM$mod3,whichSmooth=21)) %>% 
  rename(x=distMat) %>%
  ggplot(aes(x=x,y=pred,fill=sex))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.2)+
  geom_line(aes(col=sex))+geom_hline(yintercept=0,linetype='dashed')+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  labs(y='Urban Effect',fill='Sex',col='Sex',x='Distance (m)')

ggarrange(p1,p2,p3,p4,p5,p6,p7,
          labels=letters[1:7],nrow=2,ncol=4,
          legend='bottom',common.legend=T,font.label=list(size=landscapeLabel)) 

# Pardosa moesta (wolf spider) --------------------------------------------

arth %>% filter(genus=='Pardosa',species=='moesta',year==2017) %>% group_by(sex) %>% summarize(n=n())

#Do males and females have different responses to landscape?
tempArthF <- arth %>% filter(genus=='Pardosa',species=='moesta',year==2017,sex=='F') %>% group_by(BTID) %>% summarize(n=n())
tempArthM <- arth %>% filter(genus=='Pardosa',species=='moesta',year==2017,sex=='M') %>% group_by(BTID) %>% summarize(n=n())

#Trying with cereal for now. Males have stronger response, but may be effect of concurvity. 
#Does concurvity impact results from other critters? If so, it may be better to leave cereal in and let shrinkage take care of it. 
ParMoeModF <- runMods(tempArthF,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
ParMoeModM <- runMods(tempArthM,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
save(ParMoeModF, file='./data/ParMoeModF.Rdata')
save(ParMoeModM, file='./data/ParMoeModM.Rdata')
load('./data/ParMoeModF.Rdata')
load('./data/ParMoeModM.Rdata')

#Females
attach(ParMoeModF)

AIC(mod1,mod2,mod3,mod4) #Landscape model performs better
#Model 3 
summary(mod3) 
matrixplot(abs(cov2cor(sp.vcov(mod3))),c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75) #ti(TreeShrub) has NA entries
matrixplot(sp.vcov(mod3),c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75) #ti(TreeShrub) has NA entries

#Variogram of residuals - no pattern
tempTrap %>% mutate(resid=resid(mod3)) %>% 
  # group_by(ID) %>% summarize(resid=sum(abs(resid))) %>% ungroup() %>% 
  as_Spatial() %>% variogram(resid~1,.) %>% plot()

#Check k values
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F,seWithMean=T)

#Check concurvity
concurvity(mod3) #High concurvity
checkMC <- with(ParMoeModM2,concurvity(mod3,full=F)$estimate)
termNames <- rownames(checkMC)
termNames[1] <- 'Linear terms'
N <- 1:nrow(checkMC) #Plot all smooth terms
# N <- 4:nrow(checkMC) #Plot with only landscape smoothers
lN <- length(N)
matrixplot(checkMC[N,N],mar=c(1, 10, 10, 1),termNames[N])
abline(h=seq(0-1/lN/2,1+1/lN/2,length.out=1+lN/3)) #Lines to separate terms
abline(v=seq(0-1/lN/2,1+1/lN/2,length.out=1+lN/3))

detach(ParMoeModF)

#Males
attach(ParMoeModM)

AIC(mod1,mod2,mod3,mod4) #Landscape model performs better
#Model 3 
summary(mod3) 
matrixplot(abs(cov2cor(sp.vcov(mod3))),c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75) 
matrixplot(sp.vcov(mod3),c('scale',names(mod3$sp)),mar=c(1,10,10,1),numSize=0.75) 

#Variogram of residuals - no pattern
tempTrap %>% mutate(resid=resid(mod3)) %>% 
  # group_by(ID) %>% summarize(resid=sum(abs(resid))) %>% ungroup() %>% 
  as_Spatial() %>% variogram(resid~1,.) %>% plot()

#Check k values
par(mfrow=c(2,2)); gam.check(mod3); par(mfrow=c(1,1))
plot(mod3,scheme=2,shade=T,pages=1,all.terms=T,rug=F,seWithMean=T)

detach(ParMoeModM)

compareFM(ParMoeModF$mod3,ParMoeModM$mod3,vlines=2.5+rep(0:7)*3)

# Phalangium opilio -----------------------

#Many unsexed P. opilio
arth %>% filter(arthOrder=='Opiliones',year==2017) %>% group_by(sex) %>% summarize(n=n())

#Are they all from the same locations? Looks like it. Would likely have to drop time aspect of models.
arth %>% filter(arthOrder=='Opiliones',year==2017) %>% group_by(BLID,pass,sex,.drop=FALSE) %>% summarize(n=n()) %>% ungroup() %>% 
  mutate(sex=ifelse(sex=='','na',tolower(sex))) %>% 
  pivot_wider(names_from=sex,values_from=n,values_fill=list(n=0)) %>% mutate(propNA=na/(na+f+i+m)) %>%
  inner_join(select(site,BLID,lat,lon),by='BLID') %>% 
  st_as_sf(coords=c('lon','lat'),crs=4326) %>% 
  ggplot()+geom_sf(aes(col=propNA))+facet_wrap(~pass)

tempArthF <- arth %>% filter(arthOrder=='Opiliones',year==2017,sex=='F') %>% group_by(BTID) %>% summarize(n=n())
tempArthM <- arth %>% filter(arthOrder=='Opiliones',year==2017,sex=='M') %>% group_by(BTID) %>% summarize(n=n())


