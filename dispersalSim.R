#Simulation to look at dispersal distances from seminatural features, SR 2020

library(tidyverse)
library(mgcv)
library(TMB)
# library(devtools)
# install_github('https://github.com/jdyen/FREE')
# library(FREE) #Functional regression

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

#negative log-likelihood function for onion skin method
# intercept: estimate intercept? (T/F)
# lambda: fix lambda at set value (estimate if NA)
#Order of pars, () indicates optional parameters:
#(intercept),eta,rho,(lambda),(theta),(otherCoefs)
#Order of dat:
#counts, ringDist,ringMatrix,scale,(otherVars)

nllOnion <- function(pars,dat,pdf='poisson',intercept=F,lambda=NA){
  #unpack variables
  y <- dat[[1]] #counts - dependent variable
  ringDist <- dat[[2]] #distances of rings
  ringMat <- dat[[3]] #matrix of composition at distances
  scal <- dat[[4]] #vector of weights (represent size of radii)
  if(length(dat)==5) { #If another matrix provided
    otherVars <- dat[[5]] #Matrix of variables for linear regression
  } else otherVars <- NA
  
  #unpack coefficients
  if(intercept){ #If using intercept
    int <- pars[1] #First parameter is intercept
    getPars <- 2 #Offsets rank of other parameters
  } else {
    int <- 0 #No intercept
    getPars <- 1
  }
  eta <- pars[getPars] #height of increase
  getPars <- getPars + 1
  rho <- exp(pars[getPars]) #distance decay
  getPars <- getPars + 1
  if(is.na(lambda)){
    lambda <- pars[getPars] #lambda parameter
    getPars <- getPars + 1
  } 
  if(pdf=='negbin'){ #If using a negative binomial distribution
    theta <- exp(pars[getPars]) # theta
    getPars <- getPars + 1
  }
  if(is.matrix(otherVars)){ #If other variables are being used
    otherCoefs <- pars[seq(getPars,getPars+ncol(otherVars)-1,1)] #Get other coefficients
  }
  
  #Multiply matrix by scal
  ringMat <- ringMat*matrix(rep(scal,nrow(ringMat)),ncol=length(scal),byrow=T)
  
  #Expected value
  yhat <- int + (ringMat %*% negExp(ringDist,eta,rho,lambda)) + ifelse(is.matrix(otherVars),otherVars %*% otherCoefs,0) 
  
  #Link function (log)
  yhat <- exp(yhat)
  
  #Log likelihood
  nll <- switch(pdf, 
                poisson=-sum(dpois(y,yhat,log=T)),
                negbin=-sum(dnbinom(y,mu=yhat,size=theta,log=T)))
  return(nll)
}


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

#Parameters for generation - Poisson 
int <- 1.5 #intercept
dists <- seq(5,50,5)-2.5 #Distances from centre
eta <- 3
rho <- 0.1
coefs <- negExp(dists,eta,rho,2)
mu <- exp(int+(as.matrix(rings) %*% coefs))
counts <- rpois(length(mu),mu)

#Step 2: estimate coefficients

#Test nll function
datList <- list(counts,dists,as.matrix(rings),scal=rep(1,length(dists)))
nllOnion(c(1.5,3,log(0.1)),dat=datList,pdf='poisson',lambda=2,intercept=T) #Looks ok
# debugonce(nllOnion)

#Estimate
#NOTE: lambda can be estimated separately, but it's close to being non-separable. Probably better to fix term or estimate WRT other terms.
est <- optim(c(1.5,3,log(0.1)),nllOnion,dat=datList,lambda=2,intercept=T,hessian=T) 
if(any(eigen(est$hessian)$values<=0)) print('Non-positive definite Hessian')

data.frame(parName=c('int','eta','rho'),est=est$par, #Estimates look OK
           se=sqrt(diag(solve(est$hessian))),actual=c(1.5,3,log(0.1))) %>% 
  ggplot(aes(x=parName,y=est))+geom_pointrange(aes(ymax=est+se*1.96,ymin=est-se*1.96),col='red')+
  geom_point(aes(y=actual),col='black',size=5)+
  facet_wrap(~parName,scales='free')

#Bootstrap CIs for curve
bootRep <- replicate(200,expr={
  N <- length(datList[[1]])
  s <- sample(c(1:N),N,TRUE)
  dat <- list(datList[[1]][s],datList[[2]],datList[[3]][s,],datList[[4]])
  bootEst <- optim(c(1.5,3,log(0.1)),nllOnion,dat=dat,lambda=2,intercept=T,hessian=T) $par
  negExp(0:100,bootEst[2],exp(bootEst[3]),2)
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
data.frame(counts=counts,pred=exp(est$par[1]+(as.matrix(rings) %*% negExp(dists,est$par[2],exp(est$par[3]),2)))) %>% 
  ggplot(aes(x=counts+1,y=pred))+geom_point()+
  geom_abline(intercept=0,slope=1,col='red')+
  scale_x_log10()+scale_y_log10()

#Negbin version
theta <- 2.5 #Dispersion for negbin
counts <- rnbinom(length(mu),mu=mu,size=theta) #Same generating parameters, but use negative binomial instead of poisson

#Test
datList <- list(counts,dists,as.matrix(rings),scal=rep(1,length(dists)))
#(intercept),eta,rho,(lambda),(theta)
nllOnion(c(1.5,3,log(0.1),2,2.5),datList,pdf='negbin',intercept=T) #Looks ok

#Estimate
est <- optim(c(1.5,3,log(0.1),2.5),nllOnion,pdf='negbin',lambda=2,dat=datList,intercept=T,hessian=T) 
if(any(eigen(est$hessian)$values<=0)) print('Non-positive definite Hessian')

data.frame(parName=c('int','eta','rho','theta'),est=est$par, #Estimates look OK
           se=sqrt(diag(solve(est$hessian))),actual=c(1.5,3,log(0.1),log(2.5))) %>% 
  ggplot(aes(x=parName,y=est))+geom_pointrange(aes(ymax=est+se*1.96,ymin=est-se*1.96),col='red')+
  geom_point(aes(y=actual),col='black',size=5)+
  facet_wrap(~parName,scales='free')

#Bootstrap CIs for curve
bootRep <- replicate(200,expr={
  N <- length(datList[[1]])
  s <- sample(c(1:N),N,TRUE)
  dat <- list(datList[[1]][s],datList[[2]],datList[[3]][s,],datList[[4]])
  bootEst <- optim(est$par,nllOnion,pdf='negbin',lambda=2,dat=dat,intercept=T,hessian=T)$par
  negExp(0:100,bootEst[2],exp(bootEst[3]),2)
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
data.frame(counts=counts,pred=exp(est$par[1]+(as.matrix(rings) %*% negExp(dists,est$par[2],exp(est$par[3]),2)))) %>% 
  ggplot(aes(x=counts,y=pred))+geom_point()+
  geom_abline(intercept=0,slope=1,col='red')+
  scale_x_log10()+scale_y_log10()

#NegBin version using TMB
setwd('~/Documents/arthropod_analysis/TMBscripts')
compile("ringNegBin.cpp",flags="-O0 -g")
dyn.load(dynlib("ringNegBin"))

#List of data
datList <- list(y=counts, coefMat= matrix(model.matrix(~ 1,data=data.frame(counts))),
                ringDist=dists,ringMat=as.matrix(rings))

#Starting parameters
# startPars <- list(coefVec=as.vector(1.5),eta=3,logRho=0.1,lambda=2,logTheta=2.5)
startPars <- list(coefVec=1.5,eta=3,logRho=0.1,logTheta=2.5)

nllOnionTMB <- MakeADFun(data=datList, parameters=startPars, DLL="ringNegBin")

#Check function
nllOnionTMB$fn(startPars) #Works
nllOnionTMB$gr(startPars)
nllOnionTMB$he(startPars)

est <- nlminb(startPars,nllOnionTMB$fn,nllOnionTMB$gr) #Estimate parameters
estSD <- sdreport(nllOnionTMB) #SD report

#Check convergence
final_gradient <- nllOnionTMB$gr(est$par) #Gradient at convergence
if(any(abs(final_gradient)>0.0001) | estSD$pdHess==FALSE ) stop("Not converged")

est$par

str(est)
str(estSD)
estSD$value

#This isn't estimating theta correctly

summary(estSD) %>% as.data.frame() %>% 
  rename('est'='Estimate','SE'='Std. Error') %>% 
  rownames_to_column() %>% 
  mutate(upr=est+SE,lwr=est-SE) %>% 
  mutate(actual=c(1.5,3,NA,NA,0.1,2.5,coefs)) %>% 
  filter(!is.na(actual)) %>% 
  ggplot(aes(x=rowname))+geom_pointrange(aes(y=est,ymax=upr,ymin=lwr),col='red')+
  geom_point(aes(y=actual))+
  facet_wrap(~rowname,scales='free')

int <- 1.5 #intercept
dists <- seq(5,50,5)-2.5 #Distances from centre
eta <- 3
rho <- 0.1




#Extract the intercept and SE
ParHat = as.list(nllOnionTMB$SD, "Estimate")
SEHat  = as.list(nllOnionTMB$SD, "Std. Error")

opt$objective #nll
logLik(glm.fit) #ll from glm 

opt



# RFR (radial functional regression) with simulated data -------------------------------------------------

#Step 1: generate worlds

# sampleLims <- c(45,55)
sampleLims <- c(47.5,52.5)

worlds <- list() 
for(n in 1:5) worlds <- c(worlds,mkLandscape(worldLims=c(0,100),nPatches=n,nWorlds=100,leaveClear=seq(min(sampleLims),max(sampleLims),0.1)))

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
par(mfrow=c(1,1)); hist(nnDist); summary(nnDist)

#Ring composition in 5m slices
ringComp <- t(sapply(worlds,function(x) onionSkin(x,cent=50,rad=50,nSlice=10,scale=T)))
ringDist <- seq(5,50,length.out=10)-2.5
colnames(ringComp) <- paste0('d',ringDist)
par(mfrow=c(5,2)); for(i in 1:10) hist(ringComp[,i],main=i); par(mfrow=c(1,1)) #Looks OK

#Composition within a 10m-50m radius, by 5m increments
totalComp <- t(apply(ringComp,1,function(x) cummean(x)[-1]))
colnames(totalComp) <- paste0('d',seq(10,50,5))

#Step 4: compare types of regression to see which one works best
par(mfrow=c(1,1)); curve(negExp(x,3,exp(-1.5),1),0,50,ylim=c(0,3),ylab='Coef'); curve(negExp(x,3,exp(-1.5),2),0,50,add=T,col='red');
legend('topright',c('Linear','Square'),fill=c('black','red'))

datList <- list(counts=counts,ringDist=ringDist,ringComp=ringComp,scal=rep(1,length(ringDist)))

#Works
#(intercept),eta,rho,(lambda),(theta)
nllOnion(c(3,log(0.1),1),dat=datList,pdf='negbin',intercept=F,lambda=1)
# debugonce(nllOnion)

#Initial set of RFR, using only onion-ring composition
omod1 <- list()

#NegBin - no intercept - linear decay
(omod1[[1]] <- optim(c(3,log(0.1),1),nllOnion,dat=datList,pdf='negbin',intercept=F,lambda=1,hessian=T))
if(any(eigen(omod1[[1]]$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod1[[1]]$par,se=sqrt(diag(solve(omod1[[1]]$hessian)))) #Looks OK

#NegBin - no intercept - squared decay
(omod1[[2]] <- optim(c(3,log(0.1),1),nllOnion,dat=datList,pdf='negbin',intercept=F,lambda=2,hessian=T))
if(any(eigen(omod1[[2]]$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod1[[2]]$par,se=sqrt(diag(solve(omod1[[2]]$hessian)))) #Looks OK

#NegBin - intercept - linear decay
(omod1[[3]] <- optim(c(0,5.6,-1,0.8),nllOnion,dat=datList,pdf='negbin',intercept=T,lambda=1,hessian=T))
if(any(eigen(omod1[[3]]$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod1[[3]]$par,se=sqrt(diag(solve(omod1[[3]]$hessian)))) #Looks OK

#NegBin - intercept - squared decay
(omod1[[4]] <- optim(c(0,5.6,-1,0.8),nllOnion,dat=datList,pdf='negbin',intercept=T,lambda=2,hessian=T))
if(any(eigen(omod1[[4]]$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod1[[4]]$par,se=sqrt(diag(solve(omod1[[4]]$hessian)))) #Looks OK

#Second set of RFR using onion-ring composition + NNdist
omod2 <- list()
datList2 <- datList
datList2$nnDist <- matrix(scale(sqrt(nnDist)),ncol=1)

#NegBin - no intercept - linear decay
(omod2[[1]] <- optim(c(3,log(0.1),1,-1),nllOnion,dat=datList2,pdf='negbin',intercept=F,lambda=1,hessian=T))
if(any(eigen(omod2[[1]]$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod2[[1]]$par,se=sqrt(diag(solve(omod2[[1]]$hessian)))) #Looks OK

#NegBin - no intercept - squared decay
(omod2[[2]] <- optim(c(3,log(0.1),1,-1),nllOnion,dat=datList2,pdf='negbin',intercept=F,lambda=2,hessian=T))
if(any(eigen(omod2[[2]]$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod2[[2]]$par,se=sqrt(diag(solve(omod2[[2]]$hessian)))) #Looks OK

#NegBin - intercept - linear decay - NONESTIMABLE
(omod2[[3]] <- optim(c(omod1[[3]]$par,0),nllOnion,dat=datList2,pdf='negbin',intercept=T,lambda=1,hessian=T))
if(any(eigen(omod2[[3]]$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod2[[3]]$par,se=sqrt(diag(solve(omod2[[3]]$hessian)))) #Not good

#NegBin - intercept - squared decay - NONESTIMABLE
(omod2[[4]] <- optim(c(omod1[[4]]$par,0),nllOnion,dat=datList2,pdf='negbin',intercept=T,lambda=2,hessian=T))
if(any(eigen(omod2[[4]]$hessian)$values<=0)) print('Non-positive definite Hessian')
data.frame(est=omod2[[4]]$par,se=sqrt(diag(solve(omod2[[4]]$hessian)))) #Not good

#Functional regression using mgcv
datList$ringMat <- matrix(rep(ringDist,each=nrow(ringComp)),ncol=ncol(ringComp)) #Matrix of distance values
datList$nnDist <- nnDist
omod3 <- gam(counts~s(ringMat,by=ringComp),data=datList,family='nb',method='ML') #Just ring composition
omod4 <- gam(counts~nnDist+s(ringMat,by=ringComp),data=datList,family='nb',method='ML') #Ring composition + NN distance
AIC(omod3,omod4)
summary(omod3)
summary(omod4)
plot(omod3); abline(h=0,lty=2)
plot(omod4); abline(h=0,lty=2)


#Plot coefficient curves
par(mfrow=c(3,1))
plot(omod3,ylab='Coef',xlab='Dist',main='Radial Functional Regression - Radii only',se=F,col='blue',rug=F)
abline(h=0,lty=2)
curve(negExp(x,omod1[[1]]$par[1],exp(omod1[[1]]$par[2]),1),0,50,add=T) #Looks OK
curve(negExp(x,omod1[[2]]$par[1],exp(omod1[[2]]$par[2]),2),0,50,add=T,col='red') #Looks OK
curve(negExp(x,omod1[[3]]$par[2],exp(omod1[[3]]$par[3]),1)+omod1[[3]]$par[1],0,50,add=T,lty=2) #Looks OK
curve(negExp(x,omod1[[4]]$par[2],exp(omod1[[4]]$par[3]),2)+omod1[[4]]$par[1],0,50,add=T,lty=2,col='red') #Looks OK


legend('topright',c('Linear','Squared','Linear + Int','Squared + Int','Spline'),col=c('black','red','black','red','blue'),lty=c(1,1,2,2,1))
#Radial component
curve(negExp(x,omod2[[1]]$par[1],exp(omod2[[1]]$par[2]),1),0,50,xlab='Dist',ylab='Coef',main='Radial Functional Regression - Radii component') #Looks OK
curve(negExp(x,omod2[[2]]$par[1],exp(omod2[[2]]$par[2]),2),0,50,add=T,col='red')
#Linear component
yhat <- as.vector(apply(ringComp,2,mean) %*% negExp(ringDist,omod2[[1]]$par[1],exp(omod2[[1]]$par[2]),2)) + datList2$nnDist*omod2[[1]]$par[4]
plot(sort(nnDist),sort(exp(yhat)),xlab='Nearest patch',ylab='Effect',type='l',ylim=c(0.5,12),main='Radial Functional Regression - Linear component') 
yhat <- as.vector(apply(ringComp,2,mean) %*% negExp(ringDist,omod2[[2]]$par[1],exp(omod2[[2]]$par[2]),2)) + datList2$nnDist*omod2[[2]]$par[4]
lines(sort(nnDist),sort(exp(yhat)),col='red') 

#Standard nnDist GLMs - highly overdispersed, so using negbin
mod1 <-MASS::glm.nb(counts~dist,data=data.frame(count=counts,dist=nnDist)) 
par(mfrow=c(2,1)); plot(mod1,which=c(1,2)); par(mfrow=c(1,1)) #Better, but still not great
plot(predict(mod1,type='response'),counts,ylab='Actual',xlab='Predicted'); abline(0,1,col='red')

#Compare different fixed radii as predictors
mod2 <- list()
mod2[[1]] <- MASS::glm.nb(counts~ringComp[,1])
par(mfrow=c(3,2)); plot(mod2[[1]],which=c(1,2),main=paste0('Dist:',ringDist[1]+2.5))
for(i in 1:ncol(totalComp)){
  mod2[[i+1]] <- MASS::glm.nb(counts~totalComp[,i])
  plot(mod2[[i+1]],which=c(1,2),main=paste0('Dist:',ringDist[i+1]+2.5))
}
par(mfrow=c(1,1))

#Fixed radii + nnDist
mod3 <- list()
mod3[[1]] <- MASS::glm.nb(counts~ringComp[,1]+nnDist)
par(mfrow=c(3,2)); plot(mod3[[1]],which=c(1,2),main=paste0('Dist:',ringDist[1]+2.5))
for(i in 1:ncol(totalComp)){
  mod3[[i+1]] <- MASS::glm.nb(counts~totalComp[,i]+nnDist)
  plot(mod3[[i+1]],which=c(1,2),main=paste0('Dist:',ringDist[i+1]+2.5))
}
par(mfrow=c(1,1))

sapply(mod2,car::vif)

#Assemble model results
modResults <- data.frame(model=c('Onion - No Intercept - Linear','Onion - No Intercept - Square','Onion - Intercept - Linear','Onion - Intercept - Square'),
                         logLik=-sapply(omod1,function(x) x$value),nPars=c(3,3,4,4)) %>% 
  bind_rows(data.frame(model=c('Onion - Intercept - Spline','Onion - Intercept - Spline + NN'),
                       logLik=c(logLik(omod3),logLik(omod4)),nPars=c(length(coef(omod3)),length(coef(omod4))))) %>% 
  # #Linear + NN models have same logLik as models above - no improvement
  # bind_rows(data.frame(model=c('Onion - No Intercept - Linear+NN','Onion - No Intercept - Square+NN'),
  #                      logLik=-sapply(omod2,function(x) x$value)[1:2],nPars=c(4,4))) %>%
  bind_rows(data.frame(model=c('LM - Int - Nearest neighbour',paste0('LM - Int - radius',ringDist+2.5),paste0('LM - Int - radius+NN',ringDist+2.5)),
                       logLik=c(logLik(mod1),sapply(mod2,logLik),sapply(mod3,logLik)),
                       nPars=c(rep(3,1+length(mod2)),rep(4,length(mod3))))) %>%
  mutate(AIC=(-2*logLik)+(2*nPars),model=factor(model,levels=model)) %>% 
  mutate(deltaAIC=AIC-min(AIC),rank=rank(deltaAIC)) %>% separate(model,c('modType','intercept','predictor'),sep=' - ') %>% 
  mutate(predictor=ifelse(predictor=='Linear'|predictor=='Square',paste0(predictor,' decay'),predictor))
  
#Likelihood plot
plot(ringDist+2.5,sapply(mod2,logLik),ylab='log likelihood',xlab='Distance',pch=19,type='b',ylim=range(modResults$logLik))
points(ringDist+2.5,sapply(mod3,logLik),pch=19,type='b')
abline(h=logLik(mod1),lty='dashed')
text(40,logLik(mod1)-5,'Nearest patch')
text(40,logLik(mod2[[8]])+25,'Varying\nradii')
text(40,logLik(mod3[[8]])-15,'Nearest patch +\nVarying radii')
abline(h=modResults$logLik[grepl('decay',modResults$predictor)],col='red',lty='dashed')
with(filter(modResults,grepl('decay',predictor)),
     text(40,logLik,paste(predictor,intercept,sep='-'),col='red')
     )
abline(h=modResults$logLik[grepl('Spline',modResults$predictor)],col='blue',lty='dashed')
with(filter(modResults,grepl('Spline',predictor)),text(40,logLik,paste(predictor,intercept,sep='-'),col='blue'))
legend('bottom',legend=c('Standard GLM','"Onion-skin" GLM'),fill=c('black','red'), y.intersp=0.8,x.intersp=0.5)

#AIC plot
plot(ringDist+2.5,sapply(mod2,AIC),ylab='AIC',xlab='Distance',pch=19,type='b',ylim=range(modResults$AIC),main='Model Performance')
points(ringDist+2.5,sapply(mod3,AIC),pch=19,type='b')
abline(h=AIC(mod1),lty='dashed')
text(40,AIC(mod1)+5,'Nearest patch')
text(40,AIC(mod2[[8]])-45,'Varying\nradii')
text(40,AIC(mod3[[8]])+25,'Nearest patch +\nVarying radii')
abline(h=modResults$AIC[grepl('decay',modResults$predictor)],col='red',lty='dashed')
with(filter(modResults,grepl('decay',predictor)),
     text(40,AIC,paste(predictor,intercept,sep='-'),col='red')
)
abline(h=modResults$AIC[grepl('Spline',modResults$predictor)],col='blue',lty='dashed')
with(filter(modResults,grepl('Spline',predictor)),
     text(40,AIC,paste(predictor,intercept,sep='-'),col='blue')
)
legend('top',legend=c('Standard GLM','"Onion-skin" GLM'),fill=c('black','red'), y.intersp=0.8,x.intersp=0.5)

#Check results of best model (Onion Intercept Linear Decay)
omod3$par

#Bootstrap CIs for slope coefficients
bootCoef <- replicate(100,{
  samps <- sample(1:length(datList[[1]]),length(datList[[1]]),replace=T)
  tempDat <- list(datList[[1]][samps],datList[[2]],datList[[3]][samps,],datList[[4]])
  tempPars <- optim(omod3$par,nllOnion,dat=tempDat,pdf='negbin',lambda=1,intercept=T,hessian=T)$par
  #Predictions for coefs at different distances
  tempRet <- negExp(ringDist,tempPars[2],exp(tempPars[3]),1)
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
  tempPars <- optim(omod3$par,nllOnion,dat=tempDat,pdf='negbin',lambda=1,intercept=T,hessian=T)$par
  # #Predictions for coefs at different distances
  # tempRet <- lnExp(ringDist,tempPars[2],exp(tempPars[3]))
  # names(tempRet) <- colnames(ringComp)
  #Generate expected values
  tempRet <- as.vector(exp(tempPars[1] + datList[[3]] %*% negExp(ringDist,tempPars[2],exp(tempPars[3]),1)))
  # plot(as.vector(exp(tempPars[1] + datList[[3]] %*% lnExp(ringDist,tempPars[2],exp(tempPars[3])))),counts,xlab='pred',ylab='actual')
  # abline(0,1,col='red')
  return(tempRet)
  rm(samps,tempDat,tempPars,tempRet)
})

#This looks really overdispersed
data.frame(counts=counts,
           patches=rep(1:5,each=100),
           pred=as.vector(exp(omod3$par[1] + ringComp %*% negExp(ringDist,omod3$par[2],exp(omod3$par[3]),1))),
           med=apply(bootPred,1,median),
           upr=apply(bootPred,1,function(x) quantile(x,0.95)), lwr=apply(bootPred,1,function(x) quantile(x,0.05))) %>% 
  ggplot(aes(x=counts+1,y=pred))+ # geom_pointrange(aes(x=counts+0.1,y=med,ymax=upr,ymin=lwr))+
  geom_point()+
  # facet_wrap(~patches)+
  labs(x='Actual',y='Predicted')
  # geom_smooth(method=lm,formula=y~exp(x),se=F)+
  # scale_y_log10()+scale_x_log10()

#Next step: try functional regression to see if splines do any better  

# #Bootstrap CIs for slope coefficients
# bootCoef <- replicate(100,{
#   samps <- sample(1:length(datList[[1]]),length(datList[[1]]),replace=T)
#   tempDat <- list(datList[[1]][samps],datList[[2]],datList[[3]][samps,],datList[[4]])
#   tempPars <- optim(omod3$par,nllOnion,dat=tempDat,pdf='negbin',lambda=1,intercept=T,hessian=T)$par
#   #Predictions for coefs at different distances
#   tempRet <- negExp(ringDist,tempPars[2],exp(tempPars[3]),1)
#   names(tempRet) <- colnames(ringComp)
#   return(tempRet)
# })
# 
# data.frame(dist=ringDist+2.5,
#            med=apply(bootCoef,1,median),upr=apply(bootCoef,1,function(x) quantile(x,0.95)),lwr=apply(bootCoef,1,function(x) quantile(x,0.05))) %>% 
#   ggplot(aes(dist,med))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
#   geom_line()+
#   labs(x='Distance',y='Slope coefficient')
# 
# #Bootstrap CIs for predicted values
# bootPred <- replicate(100,{
#   samps <- sample(1:length(datList[[1]]),length(datList[[1]]),replace=T) #Samples to take
#   tempDat <- list(datList[[1]][samps],datList[[2]],datList[[3]][samps,],datList[[4]]) #Re-sampled data
#   tempPars <- optim(omod3$par,nllOnion,dat=tempDat,pdf='negbin',lambda=1,intercept=T,hessian=T)$par
#   # #Predictions for coefs at different distances
#   # tempRet <- lnExp(ringDist,tempPars[2],exp(tempPars[3]))
#   # names(tempRet) <- colnames(ringComp)
#   #Generate expected values
#   tempRet <- as.vector(exp(tempPars[1] + datList[[3]] %*% negExp(ringDist,tempPars[2],exp(tempPars[3]),1)))
#   # plot(as.vector(exp(tempPars[1] + datList[[3]] %*% lnExp(ringDist,tempPars[2],exp(tempPars[3])))),counts,xlab='pred',ylab='actual')
#   # abline(0,1,col='red')
#   return(tempRet)
#   rm(samps,tempDat,tempPars,tempRet)
# })
# 
# #This looks really overdispersed
# data.frame(counts=counts,
#            patches=rep(1:5,each=100),
#            pred=as.vector(exp(omod3$par[1] + ringComp %*% negExp(ringDist,omod3$par[2],exp(omod3$par[3]),1))),
#            med=apply(bootPred,1,median),
#            upr=apply(bootPred,1,function(x) quantile(x,0.95)), lwr=apply(bootPred,1,function(x) quantile(x,0.05))) %>% 
#   ggplot(aes(x=counts+1,y=pred))+ # geom_pointrange(aes(x=counts+0.1,y=med,ymax=upr,ymin=lwr))+
#   geom_point()+
#   # facet_wrap(~patches)+
#   labs(x='Actual',y='Predicted')
#   # geom_smooth(method=lm,formula=y~exp(x),se=F)+
#   # scale_y_log10()+scale_x_log10()
#   
# 
# #Try this again, but with sqrt-transformed distances
# datList <- list(counts,sqrt(ringDist),ringComp,scal=rep(1,length(ringDist)))
# 
# #NegBin - intercept - linear decay
# (omod5 <- optim(c(0,5.6,-1,0.8),nllOnion,
#                 # method='L-BFGS-B',lower=c(-10,-10,-10,-10),upper=c(10,10,10,10),
#                 dat=datList,type='ln',pdf='negbin',intercept=T,hessian=T))
# if(any(eigen(omod5$hessian)$values<=0)) print('Non-positive definite Hessian')
# data.frame(est=omod5$par,se=sqrt(diag(solve(omod5$hessian)))) #Looks OK
# curve(lnExp(x,omod5$par[2],exp(omod5$par[3]))+omod5$par[1],log(1),log(50),xlab='log(Dist)',ylab='Coef') #Looks OK
# 
# as.vector(exp(omod5$par[1] + ringComp %*% lnExp(ringDist,omod5$par[2],exp(omod5$par[3]))))
