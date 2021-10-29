# R script Summary Table FLStocks
# @author Henning Winker (JRC) & Massimiliano Cardinale (SLU)
# @email henning.winker@ec.europa.eu
# Distributed under the terms of the EUPL-1.2

## Before: data/*.RData
## After: 


mkdir("data")

library(FLCore)
library(FLBRP)

load("bootstrap/data/ices.stks.n78.Rdata", verbose=T)

# --- CREATE stock table

wkref1.stock.table = do.call(rbind,lapply(stks,function(x){
data.frame(stock=x@name,
           species=x@species,
           model = strsplit(x@desc,",")[[1]][2],
           assessed = strsplit(x@desc,",")[[1]][3],
           styr = range(x)[["minyear"]],
           endyr = range(x)[["maxyear"]],
           amin = range(x)[["min"]],
           amax = range(x)[["max"]],
           minfbar = range(x)[["minfbar"]],
           maxfbar = range(x)[["maxfbar"]],
           catch.endr = round(an(tail(catch(x))),1),  
           ssb.endyr = round(an(tail(ssb(x))),1),  
           fbar.endyr = round(an(tail(fbar(x))),3),  
           Blim = x@benchmark[["Blim"]],
           Bpa = x@benchmark[["Bpa"]],
           Btrigger = x@benchmark[["Btrigger"]],
           Fmsy =  x@benchmark[["Fmsy"]],
           Bmsy.eqsim = x@eqsim[["BMSY"]],
           B0.eqsim = x@eqsim[["B0"]],
           MSY.eqsim = x@eqsim[["Catchequi"]]
          )}
))

rownames(wkref1.stock.table) = seq(nrow(wkref1.stock.table))

# SAVE

write.taf(wkref1.stock.table, file="data/wkref1.stock.table.csv")
save(wkref1.stock.table,file="data/wkref1.stock.table.Rdata")


# --- CREATE eqsim table

# Exclude stocks with no Fmsy benchmark
excl = do.call(rbind,lapply(stks,function(x) is.na(x@benchmark[["Fmsy"]])))
stks = stks[-which(excl)]
n = length(stks)


# 1. run segmented regression on all stocks using Blim benchmark
sr <- FLSRs(lapply(stks, function(x) { 
  return(fmle(as.FLSR(x,model=segreg),fixed=list(b=as.numeric(x@benchmark["Blim"]))))
}))   

# 2. create ices ref point objects with FLBRP
brp.ices = Map(function(x, y){
  brp = brp(FLBRP(x, y))
  brp = brp + as(x@benchmark[c("Fmsy","Blim","Bpa","Btrigger")],"FLPar") # Add ICES Refpoints
  brp@name = y@name
  brp
},  stks, sr)

# EQSIM
ices.ref = do.call(rbind,lapply(stks,function(x){
  data.frame(
    stock = x@name,
    Blim = x@benchmark[["Blim"]],
    Bpa = x@benchmark[["Bpa"]],
    Btrigger = x@benchmark[["Btrigger"]],
    Bmsy = x@eqsim[["BMSY"]],
    B0 = x@eqsim[["B0"]],
    Fmsy = x@benchmark[["Fmsy"]],
    MSY = x@eqsim[["Catchequi"]],
    type = "ices")
}    
))


# SAVE

write.taf(ices.ref, file="data/ices.ref.csv")
save(ices.ref,file="data/ices.ref.Rdata")

