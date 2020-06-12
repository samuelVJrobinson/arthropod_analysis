# Helper functions for project

#Function to display matrix (m) with row/column names (nam)
# Useful for looking at correlation matrices
matrixplot <- function(m,nam=NULL,mar=NULL,numSize=1,numCol='red'){
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
  image(t(apply(m,2,rev)),axes=F,col=gray.colors(12,start=0.2,end=1,rev=T))
  
  mtext(text=nam[n:1], side=2, line=0.3, at=seq(0,1,length.out=nrow(m)), las=1, cex=0.8)
  mtext(text=nam, side=3, line=0.3, at=seq(0,1,length.out=nrow(m)), las=2, cex=0.8)
  
  # #shadowtext trick from:
  # #https://stackoverflow.com/questions/25631216/r-plots-is-there-any-way-to-draw-border-shadow-or-buffer-around-text-labels/39002911
  # theta <-  seq(0, 2*pi, length.out=60) #Vector of angles to use
  # r <- 0.15 #Proportion of character size to extend
  # x0 <- r*strwidth('0')
  # y0 <- r*strheight('0')
  
  if(numSize>0){
    #Add row text
    for(i in 1:n){ #Column of matrix to access
      cpos <- (i-1)/(n-1) #Col position in figure
      j <- 1:n #Rows of matrix to access
      rpos <- ((j-1)/(n-1))[n:1] #Column position in figure
      # #White background text
      # for(k in 1:length(theta)){
      #   text(cpos+cos(k)*x0,rpos+sin(k)*y0,round(m[j,i],2),col='white')
      # }
      text(cpos,rpos,round(m[j,i],2),col=numCol,cex=numSize) #Text
      # }
    }
  }
  par(mfrow=savepar$mfrow,mar=savepar$mar)
}

#Takes a df of arthropod counts at each site (unique to each spp), then runs 4 landscape-level models of abundance
#fitMethod = 'ML' or 'REML'
#fitFam = family of distribution to use
#basisFun = name of basis function to be used ('ts' = thin plate spline with extra shrinkage)
#doublePenalize = should gam use double-penalization? (don't use if using 'ts' or 'cr' as basis functions)
runMods <- function(tempArth,trap,nnDistMat,oRingMat2,formulas=NULL,
                    fitMethod='REML',fitFam='nb',basisFun='ts',doublePenalize=F){
  if(is.null(formulas)) stop('No gam formulas provided')
  
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
      right_join(st_drop_geometry(select(tempTrap,ID)),by='ID') %>% 
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
  
  #Assemble data list
  datList <- c(datList,oRingMat2)
  
  #Fit models
  mod1 <- gam(formula=as.formula(formulas[1]), 
              data=datList,family=fitFam,method=fitMethod,select=doublePenalize)
  cat('Finished mod1. ')
  
  mod2 <- gam(formula=as.formula(formulas[2]),
              data=datList,family=fitFam,method=fitMethod,select=doublePenalize)
  cat('Finished mod2. ')
  
  #Entire set of cover classes - "kitchen sink model"
  #Model with extra shrinkage takes longer to run. Maybe 3-5 mins?
  mod3 <- gam(formula=as.formula(formulas[3]),
              data=datList,family=fitFam,method=fitMethod,select=doublePenalize)
  cat('Finished mod3. ')
  
  #Model using "Noncrop" only:
  mod4 <- gam(formula=as.formula(formulas[4]),
              data=datList,family=fitFam,method=fitMethod,select=doublePenalize)
  cat('Finished mod4.')
  
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
  #               s(day,bs=basisFun,k=Kvals[[4]][1])+ #Thin plate spline with shrinkage
  #               s(E,N,bs=basisFun,k=Kvals[[4]][2])+
  #               ti(N,E,day,bs=basisFun,k=Kvals[[4]][3])+ #Spatiotemporal interaction   
  #               trapLoc+ 
  #               s(distMat,by=oRingNoncropProp,bs=basisFun)+ 
  #               s(distMat,by=mmat$canola,bs=basisFun)+
  #               s(distMat,by=mmat$ditch,bs=basisFun)+
  #               s(distMat,by=mmat$native,bs=basisFun)+
  #               s(distMat,by=mmat$pivot,bs=basisFun)+
  #               s(distMat,by=mmat$wetland,bs=basisFun)+
  #               ti(distMat,endDayMat,by=oRingNoncropProp,bs=basisFun) #Noncrop
  #             ,data=datList,family='nb')
  # summary(mod5)
  # plot(mod5,rug=F,scheme=2,pages=1,all.terms=T)
  # #Problem: requires 5 new terms to be written for each cover category, and is likely inestimable. Leaving this out for now.
  
  #List of objects to return
  retList <- list(datList=datList,tempTrap=tempTrap,mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4)
  return(retList)
}

#Capitalize first letter in a string (vector of strings)
firstUpper <- function(x){
  paste0(toupper(substring(x,1,1)),substring(x,2,nchar(x)))
}

#Extract data for partial effects plots of smoothing terms
# dat = dataframe of predictor data (think "newdata" from predict.gam) + column of ones with title of "by" matrix (if doing functional regression)
#   * must have the same name as predictors
# m = GAM model
# whichSmooth = smoothing terms (numeric) to use; all terms are aggregated (useful for interaction plots)
# ci = multiplicitive factor for SE bounds (default = 1.96)
smoothPred <- function(dat,m,whichSmooth,ci=1.96){ 
  #Predictor matrices
  predMat <- lapply(whichSmooth,function(i) PredictMat(m$smooth[[i]],data=dat)) #Get predictor matrices from each smoother, using dat
  predMat <- do.call('cbind',predMat) #Amalgamate into single matrix
  #Coefficients
  coefRange <- do.call('c',lapply(m$smooth[whichSmooth],function(x) x$first.para:x$last.para)) #Get coefficients to use
  coefs <- coef(m)[coefRange] #Extract coefficient values
  #Predicted value - predictor matrix X coefficients
  dat$pred <- predMat %*% coefs
  #SE - swiped from plot.gam code
  dat$se <- sqrt(pmax(0,rowSums(predMat %*% m$Vp[coefRange,coefRange] * predMat)))
  #Confidence intervals
  dat$upr <- dat$pred + dat$se*ci; dat$lwr <- dat$pred - dat$se*ci 
  return(dat) #Return entire dataframe
}

#Helper ggplot function to plot functional regression effects
#Optional legend, and different colours
effectPlot <- function(dat,leg=T,cols=NULL){
  if(sum(is.na(dat$y))==nrow(dat)){ #B/W version
    dat$y <- 'a'
    cols <- 'black'
  }
  p <- ggplot(data=dat,aes(x=x,y=pred))+
    geom_ribbon(aes(ymax=upr,ymin=lwr,fill=y),alpha=0.3,show.legend=leg)+
    geom_line(aes(col=y),show.legend=leg)+
    geom_hline(yintercept=0,linetype='dashed')+
    labs(y='Effect',fill='Date',col='Date')
  if(!is.null(cols)) { #Change colours, if necessary
    p <- p+scale_colour_manual(values=cols)+scale_fill_manual(values=cols)
  }
  return(p)
}