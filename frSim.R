#Simple example of functional regression
#Sam Robinson, Spring 2020

# Functions ---------------------------------------------------------------

library(mgcv)
logit <- function(x) log(x/(1-x))
invLogit <- function(x) exp(x)/(1+exp(x))

# Simulate data -----------------------------------------------------------

#Goal: fit yield data to time series of NDVI (or other indices)

set.seed(1)

Ndays <- 100 #Number of days of data (columns of matrix)
Nsamp <- 1000 #Number of samples (rows of matrix)

#Say that the "true" effect of NDVI on yield is something like this:
coefFun <- function(t) 3*exp(-0.03*t)
curve(coefFun,0,100,xlab='Day',ylab='Coefficient')
#i.e the effect of NDVI is positive early on, and zero later on

#Generate some NDVI data - double logistic function from doi:10.1016/j.rse.2005.10.021
ndviFun <- function(t,min,max,mS,mA,S,A){
  #t = day, min/max = min/max NDVI,
  #mS = spring slope, mA = autumn slope, S = spring start, A = autumn start
  a <- 1/(1+exp(-mS*(t-S))) #Part 1
  b <- 1/(1+exp(mA*(t-A))) #Part 2
  min + (max-min)*(a+b-1)
}
#Example
curve(ndviFun(x,min=0.1,max=1,mS=0.2,mA=0.4,S=20,A=80),0,100,xlab='Day',y='NDVI')

#Parameters for NDVI generation
ndviPars <- c(0.1,1,0.2,0.4,20,80)

#"True" NDVI values
ndviTrue <- replicate(Nsamp,{
  sapply(1:Ndays,ndviFun,min=ndviPars[1],
         max=ndviPars[2]-1/(1+exp(rnorm(1,2,2))), #Add some variation in max NDVI
         mS=ndviPars[3], mA=ndviPars[4], 
         S=ndviPars[5] + rnorm(1,0,5), #Add variation in start time
         A=ndviPars[6])
},simplify=F)

#Received NDVI values (white noise added)
ndvi <- lapply(ndviTrue,function(x){
  invLogit(logit(x)+rnorm(length(x),0,0.5))
})

#Convert both into matrices
ndviTrue <- do.call('rbind',ndviTrue)
ndvi <- do.call('rbind',ndvi)

#Example plots
par(mfrow=c(3,1))
for(i in c(1,50,100)){
  plot(1:Ndays,ndviTrue[i,],xlab='Day',ylab='NDVI',col='red',
       type='l',ylim=c(0,1),main=paste('Pixel',i))
  points(1:Ndays,ndvi[i,],pch=19)
}
par(mfrow=c(1,1))

#Generate some yield data
yield <- ndviTrue %*% coefFun(1:100)


# Fit model ---------------------------------------------------------------

#Create matrix of day values to use in model
dayMat <- matrix(rep(1:Ndays,each=Nsamp),ncol=Ndays)

#Fit model
par(mfrow=c(2,2))

mod1 <- gam(yield ~ s(dayMat,by=ndvi,bs='tp'))
summary(mod1)
plot(mod1,shade=T,rug=F,ylab='Coefficient',xlab='Day',main='Thin-plate spline')
lines(1:100,coefFun(1:100),col='red',lwd=2) #Pretty close fit
legend('topright',c('Fit','Actual'),fill=c('black','red'))
plot(predict(mod1),yield,xlab='Predicted',ylab='Actual'); abline(0,1,col='red')


mod2 <- gam(yield ~ s(dayMat,by=ndvi,bs='cr'))
summary(mod2)
plot(mod2,shade=T,rug=F,ylab='Coefficient',xlab='Day',main='Cubic regression spline')
lines(1:100,coefFun(1:100),col='red',lwd=2) #Similar fit from a cubic regression spline
legend('topright',c('Fit','Actual'),fill=c('black','red'))
plot(predict(mod1),yield,xlab='Predicted',ylab='Actual'); abline(0,1,col='red')

#Other things to try:
# Fiddle around with variances
# Fit only part of data (say, every 10 days)
