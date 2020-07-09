#Script to examine how different sexes of arthropods use the landscape differently. Marc et al 1999 note that male/female spiders have different life histories, and possibly have different habitat requirements. How does this translate into landscape patterns?

# Pterostichus melanarius ------------------------------------



# Pardosa distincta ------------------------

#Do males and females have different responses to landscape?
tempArthF <- arth %>% filter(genus=='Pardosa',species=='distincta',year==2017,sex=='F') %>% group_by(BTID) %>% summarize(n=n())
tempArthM <- arth %>% filter(genus=='Pardosa',species=='distincta',year==2017,sex=='M') %>% group_by(BTID) %>% summarize(n=n())

#Tried running with/without cereal. Very little difference, so dropping cereal terms (also, concurved).
ParDisModF <- runMods(tempArthF,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)
ParDisModM <- runMods(tempArthM,trap,nnDistMat,oRingMat2Prop,formulas=modFormulas,basisFun='ts'); beep(1)

save(ParDisModF, file='./data/ParDisModF.Rdata')
save(ParDisModM, file='./data/ParDisModM.Rdata')

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
par(mfrow=c(3,1)); for(i in 12:14) plot(mod3,scheme=2,shade=T,rug=F,seWithMean=T,select=i); par(mfrow=c(1,1))

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

detach(ParDisModM)

#Compare parameters
compareFM <- function(female,male,topTitle=''){
  #F/M intercepts
  fInt <- summary(female)[c('p.coeff','p.t','p.pv')] %>% do.call('cbind',.) %>% 
    cbind(.,p.se=unname(summary(female)[[2]][c(1:5)]),f=1)
  mInt <- summary(male)[c('p.coeff','p.t','p.pv')] %>% do.call('cbind',.) %>% 
    cbind(.,p.se=unname(summary(female)[[2]][c(1:5)]),f=0)
  
  #Intercepts
  p1 <- rbind(fInt,mInt) %>% data.frame(.) %>% rownames_to_column('loc') %>% 
    mutate(loc=gsub('.1','',loc)) %>% mutate(loc=gsub('trapLoc','',loc)) %>% 
    mutate(upr=p.coeff+p.se,lwr=p.coeff-p.se) %>% mutate(f=factor(f,labels=c('M','F'))) %>% 
    ggplot(aes(x=loc))+geom_pointrange(aes(y=p.coeff,ymax=upr,ymin=lwr,col=f))+
    labs(x='Location',y='Intercept',col='Sex')+ scale_colour_manual(values=c('blue','red'))
  
  #F/M smoothing terms
  fS <- summary(female)[c('chi.sq','s.pv','edf')] %>% do.call('cbind',.) %>% cbind(f=1)
  mS <- summary(male)[c('chi.sq','s.pv','edf')] %>% do.call('cbind',.) %>% cbind(f=0)
  
  #Smoothers
  p2 <- rbind(fS,mS) %>% data.frame %>% rownames_to_column('s') %>% 
    mutate(s=gsub('.1','',s),s=gsub('s.','',s),s=gsub('ti.','',s)) %>%
    mutate(s=factor(s,levels=s[f==1],labels=rownames(fS))) %>% 
    mutate(f=factor(f,labels=c('M','F'))) %>% 
    ggplot(aes(x=s,y=chi.sq))+geom_col(aes(fill=f),position=position_dodge())+
    #facet_wrap(~f,ncol=1)+
    theme(axis.text.x=element_text(angle=90))+
    labs(x='Smoother',y='Chi square value',fill='Sex')+
    scale_fill_manual(values=c('blue','red'))
  
  ggarrange(p1,p2,ncol=1,common.legend=TRUE,legend='right') 
}

compareFM(ParDisModF$mod3,ParDisModM$mod3)
