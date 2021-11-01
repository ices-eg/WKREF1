
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# R script to Performance Evaluation WKREF1 short cut MSE runs for selftest + 5 mp
# @author Henning Winker
# @email henning.winker@ec.europa.eu
# Licence: EUPL
# mse version: https://github.com/flr/mse/tree/FLombf
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
library(FLCore)
library(ggplotFL)
library(mse)
library(mseviz)
library(gridExtra)
library(ggplot2)


load("model/om5mps_wkref1.Rdata")

info = read.csv("data/ices.stock.info.csv")
# Subset
info = info[info$stock%in%names(om),]
# check
info$stock==names(om)

runs.dir = file.path("model", "stocks") 
hcrs = c("selftest","ices","bevholt","sb40","f01","spr40")
n = length(om)
# Derive om reference points
# List of "true" om refs
refs.eval= do.call(rbind,lapply(brp.om,function(x){
  msy <- refpts(x)["msy"]
  z = x
  refpts(z) <- FLPar(NA, dimnames=list(
    refpt=c("0.5MSY","msy","f0.1","virgin","crash"),
    quant=c("harvest", "yield", "rec", "ssb", "biomass", "revenue",
            "cost", "profit"),iter=1))
  
  refpts(z)["0.5MSY", c("harvest", "yield")] <- msy[1,c("harvest","yield")] * c(1.2, 0.5)
  blim = computeRefpts(z)
  s = x@SV["s"]
  r0.7 = an(0.2*0.7*(1-s)/(0.8*s-0.7*(s-0.2)))
  all = model.frame(refpts(x))
  data.frame(stock = x@name,
             s = s,
             msy0.5 = an(blim["0.5MSY","ssb"]),  
             bmsy0.3 = an(blim["msy","ssb"]*0.3), 
             bmsy0.5 = an(blim["msy","ssb"]*0.5), 
             b0.1 = an(blim["virgin","ssb"]*0.1), 
             r0.7 = an(r0.7*blim["virgin","ssb"]), 
             NZ = max(an(blim["msy","ssb"]*0.3),an(blim["virgin","ssb"]*0.1)),
             Bmsy = all[all$quant=="ssb","msy"],
             B0 = all[all$quant=="ssb","virgin"],
             Fmsy = all[all$quant=="harvest","msy"],
             F01 = all[all$quant=="harvest","f0.1"],
             MSY = all[all$quant=="yield","msy"],
             BBmsy = an(x@ssb.obs[,length(x@ssb.obs)])/all[all$quant=="ssb","msy"],
             FFmsy = an(x@fbar.obs[,length(x@fbar.obs)])/all[all$quant=="harvest","msy"],
             BB0 = an(x@ssb.obs[,length(x@ssb.obs)])/all[all$quant=="ssb","virgin"],
             type = "om")
}    
))

# Define Metrics for Evaluation
metrics <- list(SB = ssb, F = fbar, C = landings)


# Performance Statistics
inds=list(
  a = list(~yearMeans(F/FMSY),name="F/F[MSY]",desc="Average annual F/Fmsy"),
  b = list(~apply(iterMeans((SB/SBMSY) > 1),1,mean),name="P1(B>Bmsy)",desc="Probability that SSB > SSBmsy"),
  c = list(~yearMeans(C/MSY),name="Catch/MSY",desc="Mean Catch/MSY over years"),
  d = list(~yearMeans(iav(C)),name="AAV",desc="Average annual variation in catches"),
  
  e = list(~apply(iterMeans((SB/SBmsy0.3) < 1),1,mean),name="P1(B<0.3Bmsy)",desc="Probability that SSB < 0.3xSSBmsy"),
  f = list(~apply(iterMeans((SB/SB0.1) < 1),1,mean),name="P1(B<0.1B0)",desc="Probability that SSB < 0.1xSSB0"),
  g = list(~apply(iterMeans((SB/SBNZ) < 1),1,mean),name="P1(B<max[0.1B0,0.3Bmsy])",desc="Probability that SSB < max(0.1B0,0.3Bmsy)"),
  h = list(~apply(iterMeans((SB/SBmsy0.5) < 1),1,mean),name="P1(B<0.5Bmsy)",desc="Probability that SSB < 0.5xSSBmsy"),
  
  i = list(~apply(iterMeans((SB/SBlim) < 1),1,mean),name="P1(B<Blim)",desc="Probability that SSB < ICES Blim"),
  j = list(~apply(iterMeans(SB/SB0.5msy < 1),1,mean),name="P1(B<B[0.5MSY])",desc="Probability that SSB < SSB(0.5MSY)"),
  k = list(~apply(iterMeans((SB/SBr0.7) < 1),1,mean),name="P1(B<B[0.7rec])",desc="Probability that SSB < B at 0.7rec"),
  
  x = list(~yearMeans(SB/SBMSY),name="SSB/SSB[MSY]",desc="Average annual SSB/SSBmsy"),
  y = list(~yearMeans(F/FMSY),name="F/F[MSY]",desc="Average annual F/Fmsy")
  
  
)

# List of "true" om refs
refs = refs.eval
lref <- split(refs, seq(nrow(refs)))
names(lref) = refs$stock
ices.ref = mp.refs[mp.refs$type=="ices",]

# translate refpts to FLPars
rpts <- FLPars(lapply(lref,function(x){ 
  FLPar(
  SBMSY = x$Bmsy, 
  FMSY = x$Fmsy, 
  SBlim = ices.ref[ices.ref$stock==x$stock,"Blim"],
  SB0.5msy = x$msy0.5,
  SBmsy0.3 = x$bmsy0.3,
  SBmsy0.5 = x$bmsy0.5,
  SB0.1 = x$b0.1,
  SBr0.7 = x$r0.7,
  SBNZ = x$NZ,
  MSY = x$MSY,
  SB0 = x$B0)}))

rpts@names 


# file list from MSE run folder defined above
fileL = split(list.files(runs.dir,pattern ="mps."),seq(length(stks)))

# Load runs (Big File)
runs = lapply(fileL,function(x){load(file.path(runs.dir,x))
  mps})
names(runs) = stks@names

# List with additional data table factors
info.ls = split(info,seq(nrow(info)))

x = runs[[1]]
x@names= paste0(1:6,".",hcrs)
pb = performance(x, refpts = rpts[[x[[1]]@name]], statistics=inds, # Now statistics
                 metrics = metrics, years=list(2059:2078))

x = runs[[47]]
x[[1]]@name



# List of performance 
perfstk = Map(function(x,y){
 x@names= paste0(1:6,".",hcrs)
 pb = performance(x, refpts = rpts[[x[[1]]@name]], statistics=inds, # Now statistics
             metrics = metrics, years=list(2059:2078))
 mp = pb$mp
 pb = pb[,!"mp"]
 # Add things
 pb$common.name = y$common.name
 pb$order = y$order
 pb$resilience = y$resilience
 pb$mp = mp
 return(pb)
 },runs,info.ls)

# Save
save(perfstk,file="model/perf.69stks.ices_2710.rdata")

# Combine in data table
perf = rbindlist(c(perfstk),idcol = "stock")

#=========================================
# Compile performance statistic plots
#=========================================

# Choose orders
orders = unique(perf$order)[c(3,5,4,9)]
# resilience (r)
res = c("All stock",unique(perf$resilience)[c(3,1,2)])
kbcex =function(){theme(plot.title = element_text(size=10),
                legend.key.size = unit(0.3, 'cm'), #change legend key size
               legend.key.height = unit(0.4, 'cm'), #change legend key height
                legend.key.width = unit(0.4, 'cm'), #change legend key width
                legend.text = element_text(size=6)) #change legend text font size
}

# All + Resilience 
# Compile plots
pkbr = pfr =  list()
for(i in 1:4){
  if(i==1) df=perf
  if(i>1) df= perf[perf$resilience%in%res[i],]
  
  pkbr[[i]] = kobeMPs(df,x="x", y="y",xlim=0.5,ylim=1.5)+
    ylab(expression(F/F[MSY]))+xlab(expression(B/B[MSY]))+
    ggtitle(paste0("Kobe: ",res[i]," n = ",length(unique(df$stock))))+
    ylim(0,2.5)+kbcex()+theme()
  
  pfr[[i]] = plotBPs(df,letters[1:8],
                     target = c( a=1,b=0.5, c=1,d=0.1),
                     limit= c( e=0.05,f=0.05, g=0.05,h=0.2),
                     yminmax = c(0.05, 0.95))+theme_bw()+
    facet_wrap(~name,scales = "free_y",ncol=2)+
    ggtitle(paste0("Perfomance: ",res[i]," n = ",length(unique(df$stock))))+
    ylab("Performance statistics")
}  

## Taxomic order
i = 1
pkbo = pfo=  list()
for(i in 1:length(orders)){
  df= perf[perf$order%in%orders[i],]
  pkbo[[i]] = kobeMPs(df,x="x", y="y",xlim=0.5,ylim=1.5)+
    ylab(expression(F/F[MSY]))+xlab(expression(B/B[MSY]))+
    ggtitle(paste0("Kobe: ",orders[i]," n = ",length(unique(df$stock))))+
    ylim(0,2.5)+kbcex()+theme()
  
  pfo[[i]] = plotBPs(df,letters[1:8],
                     target = c( a=1,b=0.5, c=1,d=0.1),
                     limit= c( e=0.05,f=0.05, g=0.05,h=0.2),
                     yminmax = c(0.05, 0.95))+theme_bw()+
    facet_wrap(~name,scales = "free_y",ncol=2)+
    ggtitle(paste0("Perfomance: ",orders[i]," n = ",length(unique(df$stock))))+
    ylab("Performance statistics")
}  

#=========================================
# Plot in PDF
#=========================================
pdf(paste0("report/PerfEval_n = ",length(unique(perf$stock)),".pdf"))

# Plot All
pfr[[1]]

# Overview Resilience
grid.arrange(grobs=pkbr)

# Resilience Perf
pfr[[2]]
pfr[[3]]
pfr[[4]]

# Show by selected orders
grid.arrange(grobs=pkbo)
pfo[[1]]
pfo[[2]]
pfo[[3]]
pfo[[4]] ### More precautionary required

dev.off()


#><>><>><>><>><>><>><>><>><>><>><>><>
# Evals by stock
#><>><>><>><>><>><>><>><>><>><>><>><>

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols = c("black",gg_color_hue(6))



pdf(paste0("report/PerfByStock_n = ",length(unique(perf$stock)),".pdf"))

s = 12
for(s in 1:length(stks)){
  load(file.path(runs.dir,fileL[[s]]))
  run = mps
  run@names= paste0(1:6,".",hcrs)
  
  # Trajectories
  stk = FLStocks(c(FLStocks(om=stks[[s]]),run))
  stk@names = c("stock","selftest","ices","bevholt","sb40","f01","spr40")
  ptrj = plot(stk,metrics=list(Landings=landings,F=fbar,SSB=ssb,Recruitment=rec))+theme_bw()+ 
    geom_flpar(data=FLPars(SSB  = FLPar(Bmsy=rpts[[s]]["SBMSY"],B0=rpts[[s]]["SB0"]),
                           F = FLPar(Fmsy=rpts[[s]]["FMSY"]),Landings=FLPar(MSY=rpts[[s]]["MSY"])),x=c(2002))+
    facet_wrap(~qname, scales="free")+ggtitle(paste0(s,". ",brps[[s]][[1]]@name,": Simulations"))+
    scale_colour_manual(values=cols) # Make sure colors stay consistent by using stock = "black"
  
  print(ptrj)
  
  # Kobe
  pkb = kobeMPs(perfstk[[s]],x="x", y="y",xlim=0.5,ylim=1.5)+
    ylab(expression(F/F[MSY]))+xlab(expression(B/B[MSY]))+
    ggtitle(paste0(s,". Kobe: ",stks[[s]]@name))+
    ylim(0,2.5)
  
  
  ppf = plotBPs(perfstk[[s]],letters[1:8],
                     target = c( a=1,b=0.5, c=1,d=0.1),
                     limit= c( e=0.05,f=0.05, g=0.05,h=0.2),
                     yminmax = c(0.05, 0.95))+theme_bw()+
    facet_wrap(~name,scales = "free_y",ncol=2)+
    ggtitle(paste0(s,".Perfomance: ",stks[[s]]@name))+
    ylab("Performance statistics")
  
  print(ppf)
  print(pkb)
  
} 
dev.off()

