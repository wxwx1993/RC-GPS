library("parallel")

source（"data.generate.R")
source（"RC_IPTW_trun.R")
source（"RC.adjust.10level.R")
source（"RC_matching.R")

setwd("~/Dropbox/PM2.5 Research/Simulation/More_simulation/Simulation_plot")


subclass.MSE.bias<-mclapply(1:7,subclass.10level.fun,mc.cores=7)

IPTW.MSE.bias<-mclapply(1:7,IPTW.trun1.fun,mc.cores=7)

matching.MSE.bias<-mclapply(1:7,matching.fun,mc.cores=7)


