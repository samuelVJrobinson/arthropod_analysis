setwd("~/Projects/UofC postdoc/arthropod_analysis/TMBscripts")

library(TMB)
library(tidyverse)
theme_set(theme_classic())

#Generate data
set.seed(1234)
Nsamp <- 100
modMat <- model.matrix(~runif(Nsamp,0,1))
coefs <- c(2,1.5)
theta <- 2
yhat <- exp(modMat %*% coefs)
counts <- rnbinom(Nsamp,mu=yhat,size=theta)
plot(modMat[,2],counts)

#Create model
compile('negBinTest.cpp')
dyn.load(dynlib("negBinTest"))

datList <- list(y=counts,coefMat=modMat) #Data for model
parList <- list(coefVec=coefs,logTheta=1) #Starting parameters

obj <- MakeADFun(data=datList,parameters=parList,DLL="negBinTest")

obj$fn(obj$par)
obj$gr(obj$par)
obj$he(obj$par)

opt <- nlminb(obj$par,obj$fn,obj$gr)
optSD <- sdreport(obj)
str(opt)
str(optSD)

finalGradient <- obj$gr(opt$par)
if(any(abs(finalGradient)>0.0001) | optSD$pdHess==FALSE ) stop("Not converged")

obj$fn(opt$par) #LogLik
opt$par #Parameter estimates


data.frame(par=c('int','slope','logTheta'),actual=c(2,1.5,log(2)),est=opt$par,
           se=sqrt(diag(solve(obj$he(obj$par))))) %>% 
  mutate(upr=est+1.96*se,lwr=est-1.96*se) %>% 
  ggplot(aes(x=par))+geom_pointrange(aes(y=est,ymax=upr,ymin=lwr),col='red')+
  geom_point(aes(y=actual))





