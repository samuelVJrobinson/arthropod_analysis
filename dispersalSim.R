#Simulation to look at dispersal distances from seminatural features, SR 2020

#Idea: 2D world of habitable patches surrounded by uninhabitable matrix. 
#Organisms start in a patch, reproduce according to logistic growth equation (K related to size of patch), and then disperse according to some univariate normal distribution
#If they land in another/the same patch, they add to its population
#Basically just the theory of island biogeography + metapopulation stuff

#Logistic growth equation (n at t+1)
ntPlus <- function(nt,r,K){
  nt+(r*nt*(1-(nt/K)))
}

#Function to generate distribution of propagules from some source


mkDist <- function(patches,r,disp,T,kfac,propDisperse=1,figure=T,lastStep=F,returnDist=F){
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
  
  for(t in 1:T){
    #Reproduce
    for(i in 1:npatches){
      patches$n[i] <- round(ntPlus(patches$n[i],r,patches$size[i]))
    }
    
    #Disperse
    
    # dLocs <- with(patches,unlist(mapply(rnorm,n,mid,disp))) #Disperse from middle of patch
    dLocs <- numeric(0) #Disperse from random point within patch
    for(i in 1:npatches){
      #Dispersal location for migrants - all individuals disperse if propDisperse == 1
      
      #Normal distribution
      dLocs <- c(dLocs,with(patches,rnorm(round(n[i]*propDisperse),runif(n[i],xmin[i],xmax[i]),disp)))
      
      #Should try fat-tailed as well, as in Chapman et al 2006
      #https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2656.2006.01172.x
      
      patches$n[i] <- round(n[i]*(1-propDisperse)) #Remove migrants
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
  
  if(figure){ #Make a figure of the dispersal data
    # plot(0,type='n',xlim=range(patches[,c('xmax','xmin')]),ylim=c(0,max(patches$n)),xlab='',ylab='') #Empty plot
    # plot(density(unlist(recordDisp)),zero.line=T,main='',xlim=range(patches[,c('xmax','xmin')])) #Density plot
    hist(unlist(recordDisp),freq=F,main='',breaks=100,xlim=range(patches[,c('xmax','xmin')]))
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
r <- 1 #Rate of growth
disp <- 10 #SD of dispersal
T <- 100 #Time steps
kfac <- 30 #Factor to multiply k by (change carrying capacity for all patches)

mkDist(patches,r,disp,T,kfac,lastStep=T,propDisperse=0.2)


