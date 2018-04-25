#####################provide the main simulation results, including MSE, bias and boxplots.
library("parallel")


source("data.generate.R")
source("RC_IPTW.R")
source("RC_subclassification.R")
source("RC_matching.R")

dir.create(paste(getwd(),"/Simulation_result", sep=""))
setwd(paste(getwd(),"/Simulation_result", sep=""))


############save results for MSE and bias of each methods, along with boxplots shown in paper
subclass.MSE.bias<-mclapply(1:7,subclass.10level.fun,mc.cores=7)

IPTW.MSE.bias<-mclapply(1:7,IPTW.trun1.fun,mc.cores=7)

matching.MSE.bias<-mclapply(1:7,matching.fun,mc.cores=7)

save(subclass.MSE.bias,matching.MSE.bias,matching.MSE.bias,file="RC_GPS.MSE.bias.RData")


