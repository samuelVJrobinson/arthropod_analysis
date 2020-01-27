#Simulation to look at dispersal distances from seminatural features, SR 2020

library(tidyverse)

# Functions for simulation -------------------------------------------------------

#Idea: 2D world of habitable patches surrounded by uninhabitable matrix. 
#Organisms start in a patch, reproduce according to logistic growth equation (K related to size of patch), and then disperse according to some univariate normal distribution
#If they land in another/the same patch, they add to its population
#Basically just the theory of island biogeography + metapopulation stuff

#Function to generate distribution of propagules from some source
#patches - dataframe of patches: xmin, xmax, n = starting population
#worldLims - vector of x-limits
#r - rate of growth at n=0
#dispType - dispersal kernel ('normal','t','cauchy','ft')
#dispPars - vector of parameters for dispersal (SD, c(SD,nu) for t-dist, c(a,b) for fat-tailed)
#tim - time steps to take
#kfac - factor to multiply patch size by (raises K-constant of patch)
#propDisperse - proportion of organisms that disperse, rather than staying in the patch (between 0 and 1)
#figure - what type of figure should be displayed ('hist','dense','none')
#lastStep - should only last step be saved (T/F)
#returnDist - should locations of dispersers be returned as a vector (T/F)
mkDist <- function(patches,worldLims,r,dispType,dispPars,tim,kfac,propDisperse=1,figure='hist',lastStep=F,returnDist=F){
  #Helper functions
  
  #Logistic growth equation (n at t+1)
  ntPlus <- function(nt,r,K){
    nt+(r*nt*(1-(nt/K)))
  }
  
  #Generate random numbers from 2-sided fat tailed distribution
  #1-sided: function(d,a,b) 1/(1+a*d^b)
  rft2 <- function(n,a,b){
    d <- runif(n,0,1) #Random uniform
    posneg <- sample(c(-1,1),n,T) #Positive/negative 1
    return(posneg*((1-d)/(a*d))^(1/b))
  }
  
  #Parameters for patches
  patches$size <- with(patches,xmax-xmin) #Size of patch
  patches$mid <- rowMeans(patches[,c('xmax','xmin')]) #Middle of patch
  npatches <- nrow(patches) #number of patches
  
  # plot(0,type='n',xlim=range(patches[,c('xmax','xmin')]),ylim=c(0,30),xlab='',ylab='')
  # for(i in 1:length(patches)){
  #   with(patches,{
  #     lines(c(xmin[i],xmax[i]),c(0,0),col='forestgreen',lwd=3)
  #     points(rep(mean(c(xmin[i],xmax[i])),n[i]),1:(n[i]),pch=19)
  #   })
  # }
  
  recordDisp <- list() #List to store dispersal
  
  for(t in 1:tim){
    #Reproduce
    for(i in 1:npatches){
      patches$n[i] <- max(0,round(ntPlus(patches$n[i],r,patches$size[i]))) #Population cannot be negative
    }
    
    #Disperse
    
    # dLocs <- with(patches,unlist(mapply(rnorm,n,mid,disp))) #Disperse from middle of patch
    dLocs <- numeric(0) #Empty set to store locations
    for(i in 1:npatches){
      #Dispersal location for migrants - all individuals disperse if propDisperse == 1
      nDispersers <- round(patches$n[i]*propDisperse)
      
      dLocs <- with(patches,c(dLocs,switch(dispType, #Choose dispersal kernel
        normal=rnorm(nDispersers,runif(nDispersers,xmin[i],xmax[i]),dispPars[1]), #Normal distribution
        cauchy=rcauchy(nDispersers,runif(nDispersers,xmin[i],xmax[i]),dispPars[1]), #Cauchy
        t=rt(nDispersers,dispPars[2])*dispPars[1] + runif(nDispersers,xmin[i],xmax[i]), #Student's t (regularized)
        ft=rft2(nDispersers,dispPars[1],dispPars[2]) + runif(nDispersers,xmin[i],xmax[i]) #2-sided fat-tailed
      )))
                     
      #Should try fat-tailed as well, as in Chapman et al 2006
      #Formula f(d|a,b) = 1/(1+a*d^b); authors found a = 0.3-0.4, b = 1.5-1.6
      #Turns out this is impossible to integrate. Better to use t dist or cauchy dist.
      #https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2656.2006.01172.x
      
      patches$n[i] <- patches$n[i] - nDispersers #Remove migrants
      if(patches$n[i]<0) stop('Population negative')
    }
    
    #Assign migrants to new patches
    for(i in 1:npatches){
      #New population = residents + immigrants
      patches$n[i] <- patches$n[i] + sum(dLocs>=patches$xmin[i] & dLocs<=patches$xmax[i]) 
    }
    
    #Record location of all migrants (sucessful/unsuccessful)
    recordDisp[[t]] <- dLocs
  }
  
  if(lastStep) recordDisp <- recordDisp[[length(recordDisp)]] #Only saves the last time step
  
  recordDisp <- unlist(recordDisp) #Converts to vector
  recordDisp <- recordDisp[recordDisp>=min(worldLims)&recordDisp<=max(worldLims)] #Get rid of records outside of world range
  
  if(figure=='hist'|figure=='dens'){ #Make a figure of the dispersal data
    # plot(0,type='n',xlim=range(patches[,c('xmax','xmin')]),ylim=c(0,max(patches$n)),xlab='',ylab='') #Empty plot
    if(figure=='hist') hist(unlist(recordDisp),freq=F,main='',breaks=100,xlim=worldLims,xlab='Distance') #Histogram plot
    
    if(figure=='dens') plot(density(unlist(recordDisp)),zero.line=T,main='',xlim=worldLims,xlab='Distance') #Density plot
    
    
    for(i in 1:npatches){
      with(patches,{
        lines(c(xmin[i],xmax[i]),c(0,0),col='forestgreen',lwd=3)
        # if(n[i]>0) points(rep(mean(c(xmin[i],xmax[i])),n[i]),1:(n[i]),pch=19)
      })
    }
  }
  
  if(returnDist) return(recordDisp)

}

patches <- data.frame(xmin=c(0,30,70),xmax=c(10,40,100),n=c(5,5,15))
worldLims <- c(0,100)
r <- 1 #Rate of growth

# #normal distribution
# dispType <- 'normal'
# dispPars <- c(5)

# #t distribution
# dispType <- 't'
# dispPars <- c(5,3) #dispersal parameters

# #cauchy distribution 
# dispType <- 'cauchy'
# dispPars <- c(1) #dispersal parameters

#2-sided fat-tailed distribution
dispType <- 'ft'
dispPars <- c(0.4,1.5) #dispersal parameters

tim <- 1000 #Time steps
kfac <- 10 #Factor to multiply k by (change carrying capacity for all patches)

#Works
mkDist(patches,worldLims,r,dispType,dispPars,tim,kfac,propDisperse=0.1,figure='hist',lastStep=F,returnDist=F)

#Function to generate random 1-D landscapes
# worldLim: vector of x-limits to world
# nPatches: number of patches to generate
# nWorlds: number of random landscapes to generate (Go then. There are other worlds than these.)
# leaveClear: location(s) to check if patches intersect (reject intersecting worlds). NA indicates no patches needed
mkLandscape <- function(worldLims,nPatches,nWorlds=1,leaveClear=NA){
  #Width of world
  worldWidth <- max(worldLims)-min(worldLims)
  #Number of non-overlapping points
  nPoints <- ifelse(any(is.na(leaveClear)),0,length(leaveClear))
  
  require(gtools)
  nBreaks <- (nPatches*2)+1 #Dimensions for dirichlet process
  outputList <- list()
  listInd <- 1 #Index for list
  
  while(length(outputList)<nWorlds){
    #Generate "cut up string" of locations using dirichlet process
    breakLocs <- cumsum(c(0,rdirichlet(1,rep(1,nBreaks))))*worldWidth + min(worldLims) 
    
    #Assemble into data frame
    patches <- data.frame(xmin=breakLocs[1:length(breakLocs)-1],xmax=breakLocs[2:length(breakLocs)])
    patches <- patches[seq(2,nrow(patches)-1,2),] #Choose only "middle" patches
    rownames(patches) <- NULL #Reset row names
    
    #Check intersection criteria
    pass <- nPoints==0 || !any(rep(patches$xmax,nPoints)>=rep(leaveClear,each=nPatches) & 
                   rep(leaveClear,each=nPatches)>=rep(patches$xmin,nPoints))
    
    if(pass) { #If it passes, append to list
      outputList[[listInd]] <- patches 
      # patches$size <- patches$xmax-patches$xmin #Calculate patch size
      listInd <- listInd + 1 #Increment counter
    }
  }
  return(outputList)
}

#Function to extract "onion skin" proportions of patch cover at distances from cent
# patches: dataframe with xlims of patches
# cent: vector of centre locations
# rad: length of radius
# nSlice: number of slices to take from 0 to rad
onionSkin <- function(patches,cent,rad,nSlice,scale=F){
  #Helper function
  overlap <- function(t1,t2,b1,b2){ #How much of a target line (t1,t2) overlaps the boundary line (b1,b2)
    if(t1>=t2|b1>=b2) stop('Second argument (t2/b2) must be larger than first (t1/b1)')
    max(0,min(t2,b2)-max(t1,b1))
  }
  
  #Setup 
  nPatch <- nrow(patches) #Number of patches specified
  nCent <- length(cent) #Number of point centres to use
  wSlice <- rad/nSlice #Width of a slice
  
  #Matrix for storing values
  results <- matrix(rep(NA,nCent*nSlice),nrow=nCent,ncol=nSlice)
  rownames(results) <- cent
  colnames(results) <- paste0('d',format(round(seq(wSlice,rad,wSlice),2),nsmall=2,trim=T))
  
  for(i in 1:nCent){ #For each point centre
    for(j in 1:nSlice){ #For each slice
      #Boundaries for each slice
      left <- c(cent-wSlice*j,cent-wSlice*(j-1))
      right <- c(cent+wSlice*(j-1),cent+wSlice*j)
      results[i,j] <- 0 #Initialize to zero
      for(k in 1:nPatch){
        results[i,j] <- results[i,j] + #Add overlap 
          overlap(patches$xmin[k],patches$xmax[k],left[1],left[2]) + 
          overlap(patches$xmin[k],patches$xmax[k],right[1],right[2])
      }
    }
  }
  if(scale){ #Scales results to proportion of area in slices
    results <- results/(rad*2/nSlice)
  }
  return(results) #Output
}

#Negative exponential decay function
#d: distance
#eta: intercept at d=0
#rho: rate of decay with distance (estimate on log scale)
#lambda: power term to use with distance
negExp <- function(d,eta,rho,lambda) eta*exp(-(rho^2)*(d^lambda))
par(mfrow=c(4,1)); curve(negExp(x,3,exp(-2),2),0,50,main='Squared (lambda=2)',ylim=c(0,3),xlab='Distance',ylab='Value') #squared exponential
curve(negExp(x,3,exp(-3),2),0,50,col='red',add=T) 
curve(negExp(x,3,exp(-1),2),0,50,col='blue',add=T) 
curve(negExp(x,3,exp(-1),1),0,50,main='Linear (lambda=1)',ylim=c(0,3),xlab='Distance',ylab='Value')
curve(negExp(x,3,exp(-1.5),1),0,50,add=T,col='red')
curve(negExp(x,3,exp(0),1),0,50,add=T,col='blue') 
curve(negExp(x,3,exp(-0.5),0.5),0,50,main='Square root (lambda=0.5)',ylim=c(0,3),xlab='Distance',ylab='Value')
curve(negExp(x,3,exp(-1),0.5),0,50,add=T,col='red')
curve(negExp(x,3,exp(0.5),0.5),0,50,add=T,col='blue') 
curve(negExp(x,3,exp(-0.5),-0.5),0,50,main='Inverse square root (lambda=-0.5)',ylim=c(0,3),xlab='Distance',ylab='Value')
curve(negExp(x,3,exp(-1),-0.5),0,50,add=T,col='red')
curve(negExp(x,3,exp(0.5),-0.5),0,50,add=T,col='blue') 
par(mfrow=c(1,1))


# Simple situation --------------------------------------------------------

theme_set(theme_classic())

#1 large patch on left side
patches <- data.frame(xmin=c(0),xmax=c(20),n=c(5))
worldLims <- c(0,100)
r <- 1 #Rate of growth
tim <- 2000 #Time steps
kfac <- 10 #Factor to multiply k by

#Simulate dispersal at 3 distances (25,50,75)
dispDist <- lapply(1:100,function(x){
  data.frame(d=mkDist(patches,worldLims,r,dispType='ft',dispPars=c(0.4,1.5),tim,kfac,propDisperse=0.1,figure='none',lastStep=F,returnDist=T)) %>% 
    filter(d>20) %>% mutate(d=d-20) %>% #Removes dispersals within patch
    # ggplot(aes(x=d))+geom_histogram()
    mutate(d=case_when(d>20&d<30 ~ 'near', d>45&d<55 ~ 'mid', d>70&d<80 ~ 'far')) %>% filter(!is.na(d)) %>% 
    mutate(d=factor(d,levels=c('near','mid','far'))) %>% group_by(d,.drop=F) %>% summarize(n=n()) %>% ungroup() %>% 
    mutate(d=case_when(d=='near' ~ 25, d=='mid' ~ 50, d=='far' ~ 75)) %>% 
    mutate(rep=x)
})
dispDist <- do.call('rbind',dispDist) %>% 
  mutate(rep=factor(rep))

mod1 <- dispDist %>% mutate(d=scale(d)) %>% 
  glm(n~d,data=.,family='poisson')

mod1 <- dispDist %>% mutate(logn=log(n+0.1),d=d) %>% 
  lm(n~sqrt(d),data=.)

# mod1 <- dispDist %>% mutate(d=scale(d)) %>% 
#   glmer(n~d+(d|rep),data=.,family='poisson')
summary(mod1)
par(mfrow=c(2,1)); plot(mod1,which=c(1,2))


# Radial functional regression --------------------------------------------

#Step 1: simulate data

#Create some fake landscapes
patches <- mkLandscape(c(0,100),nPatches=1,nWorlds=30,leaveClear=50)
patches <- c(patches,mkLandscape(c(0,100),nPatches=2,nWorlds=30,leaveClear=50))
patches <- c(patches,mkLandscape(c(0,100),nPatches=3,nWorlds=30,leaveClear=50))
#Get composition in rings around site
rings <- lapply(patches,onionSkin,cent=50,rad=50,nSlice=10,scale=T)
rings <- as.data.frame(do.call('rbind',rings))
row.names(rings) <- NULL

#Parameters for generation
int <- 1.5 #intercept
dists <- seq(5,50,5)-2.5 #Distances from centre
eta <- 3
rho <- 0.1
coefs <- negExp(dists,eta,rho,2)
mu <- exp(int+(as.matrix(rings) %*% coefs))
counts <- rpois(length(mu),mu)

#Step 2: estimate coefficients

# #log-likelihood function
# ll <- function(pars,dat,type='poisson'){
#   #unpack coefficients
#   int <- pars[1] #intercept
#   eta <- pars[2] #height of increase
#   rho <- exp(pars[3]) #distance decay
#   lambda <- pars[4] #lambda
#   if(type=='negbin') theta <- pars[5]
#   #unpack dependent vars
#   y <- dat[[1]] #counts
#   d <- dat[[2]] #distances
#   mat <- dat[[3]] #matrix of composition at distances
#   
#   yhat <- exp(int + (mat %*% negExp(d,eta,rho,lambda))) 
#   nll <- switch(type,poisson=-sum(dpois(y,yhat,log=T)),
#                 negbin=-sum(dnbinom(y,mu=yhat,size=theta,log=T)))
#   return()
#   return(nll)
# }

#negative log-likelihood function for onion skin method
nllOnion <- function(pars,dat,pdf='poisson',intercept=F){
  #unpack dependent vars
  y <- dat[[1]] #counts
  d <- dat[[2]] #distances
  mat <- dat[[3]] #matrix of composition at distances
  scal <- dat[[4]] #vector of weights (represent size of radii)
  
  #unpack coefficients
  if(intercept){ #If using intercept
    int <- pars[1] #First parameter is intercept
    getPars <- 2 #Offsets rank of other parameters
  } else {
    int <- 0 #No intercept
    getPars <- 1
  }
  eta <- pars[getPars] #height of increase
  rho <- exp(pars[getPars+1]) #distance decay
  lambda <- pars[getPars+2] #lambda parameter
  
  if(pdf=='negbin'){ #If using a negative binomial distribution
    theta <- exp(pars[getPars+3]) # theta
  }
  
  #Multiply matrix by scal
  mat <- mat*matrix(rep(scal,nrow(mat)),ncol=length(scal),byrow=T)
  
  #Expected value
  yhat <- exp(int + mat %*% negExp(d,eta,rho,lambda))
  
  #Log likelihood
  nll <- switch(pdf, 
                poisson=-sum(dpois(y,yhat,log=T)),
                negbin=-sum(dnbinom(y,mu=yhat,size=theta,log=T)))
  return(nll)
}

#Test
datList <- list(counts,dists,as.matrix(rings),scal=rep(1,length(dists)))
nllOnion(c(1.5,3,log(0.1),2),datList,intercept=T) #Looks ok
# debugonce(nllOnion)

#Estimate
est <- optim(c(1.5,3,0.1,2),nllOnion,dat=datList,intercept=T,hessian=T) 
if(any(eigen(est$hessian)$values<=0)) print('Non-positive definite Hessian')

data.frame(parName=c('int','eta','rho','lambda'),est=est$par, #Estimates look OK
           se=sqrt(diag(solve(est$hessian))),actual=c(1.5,3,log(0.1),2)) %>% 
  ggplot(aes(x=parName,y=est))+geom_pointrange(aes(ymax=est+se*1.96,ymin=est-se*1.96),col='red')+
  geom_point(aes(y=actual),col='black',size=5)+
  facet_wrap(~parName,scales='free')

#Bootstrap CIs for curve
bootRep <- replicate(200,expr={
  N <- length(datList[[1]])
  s <- sample(c(1:N),N,TRUE)
  dat <- list(datList[[1]][s],datList[[2]],datList[[3]][s,])
  bootEst <- optim(c(1.5,3,0.1,2),ll,dat=dat,hessian=T)$par
  negExp(0:100,bootEst[2],exp(bootEst[3]),bootEst[4])
  })

#Coefficient curve looks OK
data.frame(dist=0:100,t(apply(bootRep,1,function(x) quantile(x,c(0.5,0.95,0.05)))),
           actual=negExp(0:100,3,0.1,2)) %>% 
  rename('med'=X50.,'upr'=X95.,'lwr'=X5.) %>% 
  ggplot(aes(x=dist))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='red')+
  geom_line(aes(y=med),col='red')+
  geom_line(aes(y=actual))+xlim(0,50)+
  labs(x='Distance',y='Coefficient')

#Predicted vs actual values - looks OK
data.frame(counts=counts,pred=exp(est$par[1]+(as.matrix(rings) %*% negExp(dists,est$par[2],exp(est$par[3]),est$par[4])))) %>% 
  ggplot(aes(x=counts+1,y=pred))+geom_point()+
  geom_abline(intercept=0,slope=1,col='red')+
  scale_x_log10()+scale_y_log10()

#Negbin version:
theta <- 2.5 #Dispersion for negbin
counts <- rnbinom(length(mu),mu=mu,size=theta) #Same generating parameters, but use negative binomial instead of poisson

#Test
datList <- list(counts,dists,as.matrix(rings),scal=rep(1,length(dists)))
nllOnion(c(1.5,3,log(0.1),2,2.5),datList,pdf='negbin',intercept=T) #Looks ok

#Estimate
est <- optim(c(1.5,3,log(0.1),2,2.5),nllOnion,pdf='negbin',dat=datList,intercept=T,hessian=T) 
if(any(eigen(est$hessian)$values<=0)) print('Non-positive definite Hessian')

data.frame(parName=c('int','eta','rho','lambda','theta'),est=est$par, #Estimates look OK
           se=sqrt(diag(solve(est$hessian))),actual=c(1.5,3,log(0.1),2,log(2.5))) %>% 
  ggplot(aes(x=parName,y=est))+geom_pointrange(aes(ymax=est+se*1.96,ymin=est-se*1.96),col='red')+
  geom_point(aes(y=actual),col='black',size=5)+
  facet_wrap(~parName,scales='free')

#Bootstrap CIs for curve
bootRep <- replicate(200,expr={
  N <- length(datList[[1]])
  s <- sample(c(1:N),N,TRUE)
  dat <- list(datList[[1]][s],datList[[2]],datList[[3]][s,],datList[[4]])
  bootEst <- optim(est$par,nllOnion,pdf='negbin',dat=dat,hessian=T)$par
  negExp(0:100,bootEst[2],exp(bootEst[3]),bootEst[4])
})

#Coefficient curve looks OK
data.frame(dist=0:100,t(apply(bootRep,1,function(x) quantile(x,c(0.5,0.95,0.05)))),
           actual=negExp(0:100,3,0.1,2)) %>% 
  rename('med'=X50.,'upr'=X95.,'lwr'=X5.) %>% 
  ggplot(aes(x=dist))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='red')+
  geom_line(aes(y=med),col='red')+
  geom_line(aes(y=actual))+xlim(0,50)+
  labs(x='Distance',y='Coefficient')

#Predicted vs actual values - looks OK
data.frame(counts=counts,pred=exp(est$par[1]+(as.matrix(rings) %*% negExp(dists,est$par[2],exp(est$par[3]),est$par[4])))) %>% 
  ggplot(aes(x=counts+1,y=pred))+geom_point()+
  geom_abline(intercept=0,slope=1,col='red')+
  scale_x_log10()+scale_y_log10()

#Having trouble with this. Seems to be right on the border of what's estimable. Perhaps intercept 
  

# RRG with simulated data -------------------------------------------------

#Step 1: generate worlds

sampleLims <- c(47.5,52.5)

worlds <- list() 
for(n in 1:5) worlds <- c(worlds,mkLandscape(worldLims=c(0,100),nPatches=n,nWorlds=100,leaveClear=sampleLims[1]:sampleLims[2]))

#Step 2: simulate dispersion and extract counts of dispersers

#Add 1 population to each world
worlds <- lapply(worlds,function(x) data.frame(x,n=1))

#Simulates dispersion, and extracts counts between sampleLims
counts <- lapply(worlds,function(x,sL) {
  temp <- mkDist(x,worldLims=c(0,100),r=1,dispType='ft',dispPars=c(0.4,1.5),tim=400,kfac=10,propDisperse=0.1,figure='none',lastStep=F,returnDist=T)
  temp <- temp[temp>=min(sL)&temp<=max(sL)]
  return(temp)},sL=sampleLims)
counts <- sapply(counts,length)

#Step 3: extract landscape data - nearest neighbour, 5m rings, and various radii

#Nearest neighbour distance
nnDist <- sapply(worlds,function(x,sL){
  x <- unlist(x[,-3])
  return(min(abs(c(sL[1]-x,sL[2]-x))))
},sL=sampleLims)
hist(nnDist)

#Ring composition in 5m slices
ringComp <- t(sapply(worlds,function(x) onionSkin(x,cent=50,rad=50,nSlice=10,scale=T)))
ringDist <- seq(5,50,length.out=10)-2.5
colnames(ringComp) <- paste0('d',ringDist)
par(mfrow=c(5,2)); for(i in 1:10) hist(ringComp[,i],main=i); par(mfrow=c(1,1)) #Looks OK

#Composition within a 10m-50m radius, by 5m increments
totalComp <- t(apply(ringComp,1,function(x) cummean(x)[-1]))
colnames(totalComp) <- paste0('d',seq(10,50,5))

#Step 4: compare types of regression to see which one works best



par(mfrow=c(1,1)); curve(negExp(x,3,exp(-2),1),0,50,ylim=c(0,3)); curve(negExp(x,3,exp(-1),1),0,50,add=T,col='red')

datList <- list(counts,ringDist,ringComp,scal=rep(1,length(ringDist)))

#Works - eta,rho,lambda,theta
nllOnion(c(3,1,1),dat=datList,pdf='negbin',intercept=F,lambda=1)
debugonce(nllOnion)

#NegBin - no intercept
(omod1 <- optim(c(3,1,1),nllOnion,
                dat=datList,pdf='negbin',intercept=F,lambda=1,hessian=T))
if(any(eigen(omod1$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod1$par,se=sqrt(diag(solve(omod1$hessian)))) #Looks OK
curve(negExp(x,omod1$par[1],exp(omod1$par[2]),omod1$par[3]),0,50,xlab='Dist',ylab='Coef') #Looks OK

#NegBin - no intercept - squared decay
(omod2 <- optim(c(3,-1,0),nllOnion,
                # method='L-BFGS-B',lower=c(-10,-10,-10),upper=c(10,10,10),
                dat=datList,type='sq',pdf='negbin',intercept=F,hessian=T))
if(any(eigen(omod2$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod2$par,se=sqrt(diag(solve(omod2$hessian)))) #Looks OK
curve(sqExp(x,omod2$par[1],exp(omod2$par[2])),0,sqrt(50),xlab='sqrt Dist',ylab='Coef',add=T,col='red') #Looks OK

#NegBin - intercept - linear decay
(omod3 <- optim(c(0,5.6,-1,0.8),nllOnion,
                # method='L-BFGS-B',lower=c(-10,-10,-10,-10),upper=c(10,10,10,10),
                dat=datList,type='ln',pdf='negbin',intercept=T,hessian=T))
if(any(eigen(omod3$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod3$par,se=sqrt(diag(solve(omod3$hessian)))) #Looks OK
curve(lnExp(x,omod3$par[2],exp(omod3$par[3]))+omod3$par[1],0,sqrt(50),xlab='sqrt Dist',ylab='Coef') #Looks OK

#NegBin - intercept - squared decay
(omod4 <- optim(c(0,5.6,-1,0.8),nllOnion,
                # method='L-BFGS-B',lower=c(-10,-10,-10,-10),upper=c(10,10,10,10),
                dat=datList,type='sq',pdf='negbin',intercept=T,hessian=T))
if(any(eigen(omod4$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod4$par,se=sqrt(diag(solve(omod4$hessian)))) #Looks OK
curve(sqExp(x,omod4$par[2],exp(omod4$par[3]))+omod3$par[1],0,sqrt(50),xlab='sqrt Dist',ylab='Coef',add=T,col='red') #Looks OK

#Standard nnDist GLMs - highly overdispersed, so using negbin
mod1 <-MASS::glm.nb(counts~dist,data=data.frame(count=counts,dist=sqrt(nnDist))) 
par(mfrow=c(2,1)); plot(mod1,which=c(1,2)); par(mfrow=c(1,1)) #Better, but still not great

#Compare different fixed radii as predictors
mod2 <- list()
mod2[[1]] <- MASS::glm.nb(counts~ringComp[,1])
par(mfrow=c(2,1)); plot(mod2[[1]],which=c(1,2))
for(i in 1:ncol(totalComp)){
  mod2[[i+1]] <- MASS::glm.nb(counts~totalComp[,i])
  plot(mod2[[i+1]],which=c(1,2))
}
par(mfrow=c(1,1))

#Assemble model results
modResults <- data.frame(model=c('Onion - No Intercept - Linear','Onion - No Intercept - Square','Onion - Intercept - Linear','Onion - Intercept - Square'),
                         logLik=sapply(list(omod1,omod2,omod3,omod4),function(x) -x$value),
                         nPars=c(3,3,4,4)) %>% bind_rows(
                           data.frame(model=c('LM - Int - Nearest neighbour',paste0('LM - Int - radius',ringDist+2.5)),
                                      logLik=c(logLik(mod1),sapply(mod2,logLik)),nPars=rep(3,1+length(mod2)))
                         ) %>% mutate(AIC=(-2*logLik)+(2*nPars),model=factor(model,levels=model)) %>% 
  mutate(deltaAIC=AIC-min(AIC),rank=rank(deltaAIC)) %>% separate(model,c('modType','intercept','predictor'),sep=' - ') %>% 
  mutate(predictor=ifelse(predictor=='Linear'|predictor=='Square',paste0(predictor,' decay'),predictor))
  
#Likelihood plot
plot(ringDist+2.5,sapply(mod2,logLik),ylab='log likelihood',xlab='Distance',pch=19,type='b',ylim=range(modResults$logLik))
abline(h=logLik(mod1),lty='dashed')
text(40,logLik(mod1),'Nearest patch')
text(40,logLik(mod2[[8]])+25,'Varying\nradii')
abline(h=modResults$logLik[grepl('Onion',modResults$modType)],col='red',lty='dashed')
# text(40,-1465,'"Onion-ring" models',col='red')
with(modResults,
     text(40,logLik[grepl('Onion',modType)],paste(predictor[grepl('Onion',modType)],intercept[grepl('Onion',modType)],sep='-'),col='red')
     )
legend('bottom',legend=c('Standard GLM','"Onion-skin" GLM'),fill=c('black','red'), y.intersp=0.8,x.intersp=0.5)

#Check results of best model (Onion Intercept Linear Decay)
omod3$par

#Bootstrap CIs for slope coefficients
bootCoef <- replicate(100,{
  samps <- sample(1:length(datList[[1]]),length(datList[[1]]),replace=T)
  tempDat <- list(datList[[1]][samps],datList[[2]],datList[[3]][samps,],datList[[4]])
  tempPars <- optim(omod3$par,nllOnion,dat=tempDat,type='ln',pdf='negbin',intercept=T,hessian=T)$par
  #Predictions for coefs at different distances
  tempRet <- lnExp(ringDist,tempPars[2],exp(tempPars[3]))
  names(tempRet) <- colnames(ringComp)
  return(tempRet)
})

data.frame(dist=ringDist+2.5,
           med=apply(bootCoef,1,median),upr=apply(bootCoef,1,function(x) quantile(x,0.95)),lwr=apply(bootCoef,1,function(x) quantile(x,0.05))) %>% 
  ggplot(aes(dist,med))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+
  labs(x='Distance',y='Slope coefficient')

#Bootstrap CIs for predicted values
bootPred <- replicate(100,{
  samps <- sample(1:length(datList[[1]]),length(datList[[1]]),replace=T) #Samples to take
  tempDat <- list(datList[[1]][samps],datList[[2]],datList[[3]][samps,],datList[[4]]) #Re-sampled data
  tempPars <- optim(omod3$par,nllOnion,dat=tempDat,type='ln',pdf='negbin',intercept=T,hessian=T)$par
  # #Predictions for coefs at different distances
  # tempRet <- lnExp(ringDist,tempPars[2],exp(tempPars[3]))
  # names(tempRet) <- colnames(ringComp)
  #Generate expected values
  tempRet <- as.vector(exp(tempPars[1] + datList[[3]] %*% lnExp(ringDist,tempPars[2],exp(tempPars[3]))))
  # plot(as.vector(exp(tempPars[1] + datList[[3]] %*% lnExp(ringDist,tempPars[2],exp(tempPars[3])))),counts,xlab='pred',ylab='actual')
  # abline(0,1,col='red')
  return(tempRet)
  rm(samps,tempDat,tempPars,tempRet)
})

#This looks really overdispersed
data.frame(counts=counts,
           patches=rep(1:5,each=100),
           pred=as.vector(exp(omod3$par[1] + ringComp %*% lnExp(ringDist,omod3$par[2],exp(omod3$par[3])))),
           med=apply(bootPred,1,median),
           upr=apply(bootPred,1,function(x) quantile(x,0.95)), lwr=apply(bootPred,1,function(x) quantile(x,0.05))) %>% 
  ggplot(aes(x=counts+1,y=pred))+ # geom_pointrange(aes(x=counts+0.1,y=med,ymax=upr,ymin=lwr))+
  geom_point()+
  facet_wrap(~patches)+
  labs(x='Actual',y='Predicted')
  # geom_smooth(method=lm,formula=y~exp(x),se=F)+
  # scale_y_log10()+scale_x_log10()
  

#Try this again, but with sqrt-transformed distances
datList <- list(counts,sqrt(ringDist),ringComp,scal=rep(1,length(ringDist)))

#NegBin - intercept - linear decay
(omod5 <- optim(c(0,5.6,-1,0.8),nllOnion,
                # method='L-BFGS-B',lower=c(-10,-10,-10,-10),upper=c(10,10,10,10),
                dat=datList,type='ln',pdf='negbin',intercept=T,hessian=T))
if(any(eigen(omod5$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod5$par,se=sqrt(diag(solve(omod5$hessian)))) #Looks OK
curve(lnExp(x,omod5$par[2],exp(omod5$par[3]))+omod5$par[1],log(1),log(50),xlab='log(Dist)',ylab='Coef') #Looks OK

as.vector(exp(omod5$par[1] + ringComp %*% lnExp(ringDist,omod5$par[2],exp(omod5$par[3]))))
