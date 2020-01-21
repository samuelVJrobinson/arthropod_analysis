#Simulation to look at dispersal distances from seminatural features, SR 2020

library(tidyverse)

# Set up simulation -------------------------------------------------------

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
      patches$n[i] <- round(ntPlus(patches$n[i],r,patches$size[i]))
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

mkDist(patches,worldLims,r,dispType,dispPars,tim,kfac,propDisperse=0.1,figure='hist',lastStep=F,returnDist=F)

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

#Function to generate random 1-D landscapes
# worldLim: vector of x-limits to world
# nPatches: number of patches to generate
# nWorlds: number of random landscapes to generate (Go then. There are other worlds than these.)
# leaveClear: location to check if patches intersect (reject intersecting worlds)
mkLandscape <- function(worldLims,nPatches,nWorlds=1,leaveClear=NA){
  #Width of world
  worldWidth <- max(worldLims)-min(worldLims)
  
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
    
    #Check intersection criteria (if present)
    if(is.na(leaveClear) | !any(patches$xmax>=leaveClear & leaveClear>=patches$xmin)) {
      patches$size <- patches$xmax-patches$xmin #Calculate patch size
      outputList[[listInd]] <- patches #If it passes, append to list
      listInd <- listInd + 1 #Increment counter
    }
  }
  return(outputList)
}

#Works
mkLandscape(worldLims=c(0,100),nPatches=2,nWorlds=3,leaveClear=50)










