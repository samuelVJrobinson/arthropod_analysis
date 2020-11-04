# EXAMPLES OF FUNCTIONAL REGRESSION PLOTS
# USED FOR BOX 1 IN PAPER, EXPLAINING HOW TO INTERPRET FUNCTIONAL REGRESSION PLOTS
# SR FALL 2020

# Load everything ---------------------------------------------------------
library(tidyverse)
library(ggpubr)
theme_set(theme_classic()) #Classic theme

#Functions to simulate

#Line
lineFun <- function(int,slope,xRange=c(-1,1),n=100){
  x <- seq(min(xRange),max(xRange),length.out=n)
  y <- int + slope*x
  return(y)
}
plot(lineFun(-1,0.001,xRange=c(0,1000))) #Works

#Exponential
expFun <- function(a,b,d,xRange=c(0,1),n=100){
  x <- seq(min(xRange),max(xRange),length.out=n)
  y <- a*(1-exp(-b/x))+d
  return(y)
}
plot(expFun(1,0.05,-0.1)); abline(h=0)

#Set seed
set.seed(1)

#Positively sloped line
p1 <- apply(mapply(lineFun,int=rnorm(1000,1,0.1),slope=rnorm(1000,0.5,0.1),MoreArgs=list(n=1000)),
      1,function(x) quantile(x,c(0.05,0.5,0.95))) %>% t(.) %>% 
  data.frame(x=seq(0,1,length.out=length(.[,1])),.) %>% 
  rename(lwr=`X5.`,med=`X50.`,upr=`X95.`) %>% 
  ggplot(aes(x,med))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+
  labs(y='Activity Density (y)',x='Proportion cover')+#ylim(-1.5,1.5)+
  annotate('segment',x=0.75,xend=0.75,y=0.5+0.75,yend=0.5+1,linetype='dashed')+
  annotate('segment',x=0.75,xend=1,y=0.5+1,yend=0.5+1,linetype='dashed')+
  annotate('segment',x=0.5,xend=0.7,y=0.7+1,yend=0.55+1,arrow=arrow())+
  annotate('text',x=0.4,y=1.7,label='beta == 0.5',parse=TRUE)+
  annotate('text',x=0.25,y=1.2,label='hat(y) == alpha + beta %*% Proportion~Cover',parse=TRUE)
  
#Actual function: y = -0.5 + 1 * x
#Inverse: x + 0.5 = y

#Flat line
p2 <- apply(mapply(lineFun,int=rnorm(1000,0.5,0.1),slope=rnorm(1000,0,0.1),MoreArgs=list(n=1000)),
            1,function(x) quantile(x,c(0.05,0.5,0.95))) %>% t(.) %>% 
  data.frame(x=seq(0,1000,length.out=length(.[,1])),.) %>% 
  rename(lwr=`X5.`,med=`X50.`,upr=`X95.`) %>% 
  ggplot(aes(x,med))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+geom_hline(yintercept=0,linetype='dashed')+
  annotate('text',x=500,y=-0.75,label='Effect of cover is positive at all distances')+
  annotate('segment',x=400,xend=300,y=-0.5,yend=0.2,arrow=arrow())+
  annotate('segment',x=600,xend=700,y=-0.5,yend=0.2,arrow=arrow())+
  labs(x='Distance from sampling location',y=expression(paste('Slope of proportion cover (',beta,')')))+ylim(-1.5,1.5)

#Exponential
p3 <- apply(mapply(expFun,a=rnorm(1000,1,0.2),b=rnorm(1000,0.05,0.002),d=rnorm(1000,-0.02,0.1),MoreArgs=list(n=1000)),
            1,function(x) quantile(x,c(0.05,0.5,0.95))) %>% t(.) %>% 
  data.frame(x=seq(0,1000,length.out=length(.[,1])),.) %>% 
  rename(lwr=`X5.`,med=`X50.`,upr=`X95.`) %>% 
  ggplot(aes(x,med))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+geom_hline(yintercept=0,linetype='dashed')+
  annotate('text',x=700,y=1,label='Positive effect of nearby cover')+
  annotate('segment',x=350,xend=100,y=1,yend=0.6,arrow=arrow())+
  annotate('text',x=500,y=-0.75,label='Neutral effect of far-away cover')+
  annotate('segment',x=600,xend=600,y=-0.5,yend=-0.2,arrow=arrow())+
  labs(x='Distance from sampling location',y=expression(paste('Slope of proportion cover (',beta,')')))+ylim(-1.5,1.5)

#Negatively sloped line
p4 <- apply(mapply(lineFun,int=rnorm(1000,0,0.1),slope=rnorm(1000,-0.5,0.1),MoreArgs=list(n=1000)),
      1,function(x) quantile(x,c(0.05,0.5,0.95))) %>% t(.) %>% 
  data.frame(x=seq(0,1000,length.out=length(.[,1])),.) %>% 
  rename(lwr=`X5.`,med=`X50.`,upr=`X95.`) %>% 
  ggplot(aes(x,med))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+geom_hline(yintercept=0,linetype='dashed')+
  annotate('text',x=700,y=1,label='Positive early effect of cover')+
  annotate('segment',x=450,xend=200,y=1,yend=0.6,arrow=arrow())+
  annotate('text',x=300,y=-1,label='Negative late effect cover')+
  annotate('segment',x=525,xend=800,y=-1,yend=-0.6,arrow=arrow())+
  scale_x_continuous(breaks=seq(0,1000,250),labels=seq(100,200,25))+
  labs(x='Day of year',y=expression(paste('Slope of proportion cover (',beta,')')))+ylim(-1.5,1.5)

pAll <- ggarrange(p1,p2,p3,p4,labels=letters[1:4],nrow=2,ncol=2) #Plot spatial/temporal effects
ggsave('./figures/frExample.png',pAll,height=8,width=10)
