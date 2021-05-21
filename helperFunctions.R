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
runMods <- function(tempTrap,nnDistMat,oRingMat2,formulas=NULL,
                    fitMethod='REML',fitFam='nb',basisFun='ts',doublePenalize=FALSE){
  if(is.null(formulas)) stop('No gam formulas provided')
  
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
  retList <- list(datList=datList,#tempTrap=tempTrap,
                  mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4)
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

#Compare parameters in male/female models
compareFM <- function(female,male,topTitle='',vlines=NA){
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
    geom_vline(xintercept=vlines)+
    theme(axis.text.x=element_text(angle=90))+
    labs(x='Smoother',y='Chi square value',fill='Sex')+
    scale_fill_manual(values=c('blue','red'))
  
  ggarrange(p1,p2,ncol=1,common.legend=TRUE,legend='right') 
}

#Function to extract marginal variance from terms from a GAM (a la Nakagawa and Schielzeth 2013)
getR2Terms <- function(mod,intVal,thetaVal,individual=FALSE){ #Model object, plus intercept value from null model.
  #Note from Appendix: intVal (exp(beta0)) should be obtained either from a model with centred or scaled variables (sense Schielzeth 2010), or an intercept‐only model while including all random effects. Not clear whether this should be done with theta as well, but I assume so.
  #Assumes models follow same structure as used in frAnalysis.R
  
  #Distribution-specific variance:
  
  #For a poisson distribution (see Appendix 1 in https://doi.org/10.1111/j.2041-210x.2012.00261.x):
  # E(x) = lambda
  # var(x) = lambda
  # var(ln(x)) = ln(1 + var(x)/lambda^2)
  #            = ln(1 + lambda/lambda^2)
  #            = ln(1 + 1/lambda)
  
  #For a negative binomial distribution (see Appendix 1 in https://doi.org/10.1101/095851):
  # E(x) = mu
  # var(x) = mu(1+mu/theta)
  # var(ln(x)) = ln(1 + var(x)/mu^2)
  #            = ln(1 + mu(1+mu/theta)/mu^2)
  #            = ln(1 + 1/mu + 1/theta)
  #            = ln(1 + 1/exp(Intercept) + 1/theta)
  
  c1 <- mod$coefficients #Get coefs
  m1 <- model.matrix(mod) #Get model matrix
  
  #Linear terms
  lt <- attributes(mod$pterms)$term.labels
  ltVar <- sapply(lt,function(x) var(m1[,grepl(x,colnames(m1))] %*% c1[grepl(x,names(c1))])) 
  
  #Smooth terms
  st <- sapply(mod$smooth,function(x) x$label) #Get smoothing terms
  st <- unique(gsub('(s|ti)\\((distMat|endDayMat|distMat,endDayMat)\\)\\:','',st)) #Reduce to unique terms
  sVar <- sapply(st,function(x) var(m1[,grepl(x,colnames(m1),fixed=TRUE)] %*% c1[grepl(x,names(c1),fixed=TRUE)]))
  
  #Variance for each model component (sigma^2_f)
  sigmaF <- data.frame(terms=c(lt,st),var=c(ltVar,sVar),row.names=NULL) #"Fixed" terms
  
  #Distribution-specific variance (sigma^2_d) - assumes intercept is 0
  sigmaD <- log(1+(1/intVal)+(1/thetaVal))
  
  if(individual){ #Return individual variance components, and let the user figure it out
    return(rbind(sigmaF,data.frame(terms='sigmaD',var=sigmaD)))
  } else {
    r2 <- c(sum(sigmaF$var[!grepl('(day|E,N)',sigmaF$terms)])/(sum(sigmaF$var)+sigmaD), #"Marginal" R2 (no spatial/temporal smoothers)
            sum(sigmaF$var)/(sum(sigmaF$var)+sigmaD)) #Conditional R2
    names(r2) <- c('Marginal','Conditional')
    return(r2)
  }
}

#Capitalize first letter in a string
firstUpper <- function(x){
  paste0(toupper(substring(x,1,1)),substring(x,2,nchar(x)))
}

#Rounds to dig digits, converts to character, anything less than digits becomes "<0.xx1" (useful for tables)
roundLess <- function(x,dig){
  x <- as.character(round(x,dig))
  x[x=='0'] <- paste0('<0.',strrep('0',dig-1),'1',collapse='')
  return(x)
}

#Makes text bold in rMarkdown
boldRMD <- function(x) paste0('**',x,'**') 

#Makes text bold in LaTeX
boldLaTeX <- function(x) paste0('\\textbf{',x,'}') 

# Convenience function to make FR plots from model structure
#mod = model
#term = which landscape term? (eg "Pulses")
#type = both, space, time
#days = range of days to use
#dists = range of distances to use
#ylab = custom y label
#showLegend = show legend for "both" plots?
makeFRplot <- function(mod,term=NULL,type=NULL,days=NULL,dists=seq(30,1500,30),ylab=NULL,dateLabs=NULL,showLegend=F){
  if(xor(is.null(days),is.null(dateLabs))) stop('Must specify days and dateLabs, or neither')
  if(is.null(days)&is.null(dateLabs)){
    if(type=='time'){
      days <- 149:241
    } else {
      days <- c(173,232)
      dateLabs <- c('Early','Late')
    }
  }
  smLabs <- sapply(mod$smooth,function(x) x$label) #Smoother labels
  whichSmooths <- which(grepl(term,smLabs)) #Select terms
  if(is.null(ylab)) ylab <- paste0(term,' effect')
  if(type=='time') {
    whichSmooths <- whichSmooths[2]
    d <- data.frame(endDayMat=days,Temp=1,y=NA) 
    names(d)[2] <- term  
    p <- d %>% smoothPred(mod,whichSmooth=whichSmooths) %>% rename(x=endDayMat) %>%
      mutate(x=as.Date(paste0(x,'-2017'),format='%j-%Y')) %>%
      effectPlot(leg=showLegend)+labs(x='Time of year',y=ylab)
  } else if(type=='space') {
    whichSmooths <- whichSmooths[1]
    d <- data.frame(distMat=dists,Temp=1,y=NA) 
    names(d)[2] <- term
    p <- d %>% smoothPred(mod,whichSmooth=whichSmooths) %>% rename(x=distMat) %>% 
      effectPlot(leg=showLegend)+labs(x='Distance from trap location (m)',y=ylab)
  } else {
    d <- expand.grid(endDayMat=days,distMat=dists,Temp=1) 
    names(d)[3] <- term
    p <- d %>% smoothPred(mod,whichSmooth=whichSmooths) %>% rename(x=distMat,y=endDayMat) %>%
      mutate(y=factor(y,labels=dateLabs)) %>%
      effectPlot(leg=showLegend,cols=smoothCols)+labs(x='Distance from trap location (m)',y=ylab)
  }
  return(p)
}

#Get compact letter display from GAM linear terms
#Idea from: 
# https://stats.stackexchange.com/questions/376237/correcting-for-multiple-pairwise-comparisons-with-gam-objects-mgcv-in-r
cldGam <- function(mod){ 
  require(multcomp)
  coefNames <- names(coef(mod)[1:5])
  comparisons <- rep('',sum(seq(length(coefNames)-1)))
  n <- 1
  for(i in 1:(length(coefNames)-1)){
    for(j in (i+1):length(coefNames)){
      comparisons[n] <- paste0(coefNames[i],' - ',coefNames[j],' = 0') 
      n <- n+1
    }
  }
  modcomp <- glht(mod,linfct=comparisons)
  # summary(mod3comp)
  modcomp$focus <- 'trapLoc' #Trick to get CLD from GAM
  letDisp <- cld(modcomp)$mcletters$Letters
  data.frame(trapLoc=names(letDisp),labs=unname(letDisp))
}

