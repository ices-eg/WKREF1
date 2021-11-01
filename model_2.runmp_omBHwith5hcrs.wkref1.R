
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# R script to run short-mse for WKREF1
# @author Henning Winker
# @email henning.winker@ec.europa.eu
# Licence: EUPL
# mse version: https://github.com/flr/mse/tree/FLombf
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


library(FLCore) # Please update
library(ggplotFL)
library(mse)
library(gridExtra)
library(FLBRP)

source("utilities.R")

# Runs bevHolt OM with 3 hrcs
load(file="model/om5mps_wkref1.Rdata",verbose=T)

stks@names==names(om)
names(sc.ctrl$om) == names(om)

library(doParallel)

# WARNING: Choose an appropriate number of cores
registerDoParallel(35)

# 
runs <- foreach(s = 1:length(stks)) %dopar% {
  # Only run if not exist already
  if(file.exists(paste0("model/mps.",om[[s]]@stock@name,".Rdata"))==F){
  cat("Running",om[[s]]@stock@name)
  mps = FLStocks(lapply(sc.ctrl,function(x){ # for each here?
    out = mp(om[[s]], ctrl=x[[s]], args=om[[s]]@stock@args,parallel=FALSE)
    stock(out@om)
  }))
  mps@names = paste0(om[[s]]@stock@name,"_",hcr.names) 
  save(mps,file=paste0("model/stocks/mps.",om[[s]]@stock@name,".Rdata"))
  return(mps)
  }
} # end of loop 

save(runs, file="model/runs.Rdata", compress="xz")
