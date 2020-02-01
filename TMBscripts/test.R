library(TMB)
#NegBin version using TMB
compile("./TMBscripts/ringNegBin.cpp")
dyn.load(dynlib("./TMBscripts/ringNegBin"))

#List of data
datList <- list(y=counts, coefMat= model.matrix(~ 1,data=data.frame(counts)),
                ringDist=dists,ringMat=as.matrix(rings))

#Starting parameters
startPars <- list(coefVec=1.5,eta=3,logRho=log(0.1),lambda=2,logTheta=2.5)

nllOnionTMB <- MakeADFun(data=datList, parameters=startPars, DLL="ringNegBin")

#Check function
nllOnionTMB$fn(nllOnionTMB$par)
nllOnionTMB$gr(nllOnionTMB$par)