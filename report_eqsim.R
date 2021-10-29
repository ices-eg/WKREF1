# report_eqsim.R - DESC
# /report_eqsim.R

# Copyright Iago MOSQUEIRA (WMR), 2021
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(ggplotFL)
library(FLBRP)
library(gridExtra)

load("data/ices.ref.Rdata", verbose=TRUE)

# Plot by Stock

pdf("report/wkref1.eqsim.pdf")
for(i in 1:length(stks)){
stk = stks[[i]]

# BRPs
rcol=rainbow(5)[c(3,5,4,1,2)]
brp.ices[[i]]
ref = plot(brp.ices[[i]],obs=T,refpts=c("virgin", "Blim","Bpa","Fmsy","f0.1"),colours=rcol,panels=c(2,4),ncol=1)
print(ref+ggtitle(paste0(i,".",stk@name,": Recruitment and Production Function"))+theme_bw())
# Status
end = dims(stk)[["maxyear"]]
st = dims(stk)[["minyear"]]+2
lab1 = st+0.03*length(st:end)
lab2 = st+0.25*length(st:end)
lab3 = st+0.5*length(st:end)

status = plot(stk,metrics=list(Landings=landings,F=fbar,SSB=ssb,Recruitment=rec))+ 
  geom_flpar(data=FLPars(SSB  = FLPar(Blim=ices.ref[i,"Blim"],Bpa=ices.ref[i,"Bpa"],Btrigger=ices.ref[i,"Btrigger"],Bmsy=ices.ref[i,"Bmsy"],B0=ices.ref[i,"B0"]),
                         F = FLPar(Fmsy=ices.ref[i,"Fmsy"]),Landings=FLPar(MSY=ices.ref[i,"MSY"])),
  x=c(lab1,lab2,lab3,lab1,lab1,lab3,lab3),colour = c("red","darkorange","orange","darkgreen","blue","darkgreen","darkgreen"))+
  facet_wrap(~qname, scales="free")
print(status+ggtitle(paste0(i,".",stk@name,": EQSIM"))+theme_bw())
}

dev.off()


pdf("report/wkref1.stock.checks.pdf")
for(i in 1:length(stks)){
  stk = stks[[i]]
  
   advice = plot(stk,metrics=list(Catch=catch,Fbar=fbar,SSB=ssb,Recruitment=rec))+ 
    geom_flpar(data=FLPars(Fbar = FLPar(Fmsy=stk@benchmark[["Fmsy"]]),
                           SSB = FLPar(Blim=stk@benchmark[["Blim"]],Bpa=stk@benchmark[["Bpa"]])),
               x=c(2007,2003,2011),colour=c("blue","red","orange"))+
    facet_wrap(~qname, scales="free")+xlab("Year")+ggtitle(paste0(i,".",stk@name,": ICES Advice, ",strsplit(stk@desc,",")[[1]][2]," (",stk@species,")"))+theme_bw()+ylab("Quantity")
  print(advice)
  
  # Show some Stock Dynamics at age
  dat=as.data.frame(FLQuants(stk,Numbers=stock.n, Weight=catch.wt,M=m,Maturity=mat),drop=T)
  dat$Age = factor(dat$age)
  # Plotting dynamics at age
  bio =ggplot(dat)+
    geom_line(aes(year,data,group=age,col=Age))+
    facet_wrap(~qname,scale="free",ncol=2)+
    theme_bw()+
    xlab("Year")+xlab("Quantities")+ggtitle(paste0(i,".",stk@name,": Biology (",stk@species,")"))+
    theme(legend.position="right",legend.text = element_text(size=8),
          legend.title = element_text(size=9),legend.key.height = unit(.5, 'cm'))

  print(bio)
 
  # By age
  dat = as.data.frame(
    FLQuants(M=m(stk),Weight=stock.wt(stk),
             F=harvest(stk),Selectivity=catch.sel(stk))
  )
  dat$Year = factor(dat$year)
  byage =ggplot(dat)+
    geom_line(aes(age,data,group=Year,col=Year))+
    facet_wrap(~qname,scale="free",ncol=2)+
    theme_bw()+theme(legend.position="none")+
    xlab("Age")+ylab("Value-at-age")+ggtitle(paste0(i,".",stk@name,": Stock-at-age (",stk@species,")"))
  print(byage)
  
  
  # Fishery 
  dat = as.data.frame(
    FLQuants(Biomass=stock.n(stk)*stock.wt(stk), Vuln.Bio=stock.n(stk)*stock.wt(stk)*catch.sel(stk),SSB=stock.n(stk)*stock.wt(stk)*mat(stk),
             Catch = catch.n(stk),F=harvest(stk),
             "Selectivity"=catch.sel(stk)))
  
  dat$Age = factor(dat$age)
  # Plotting dynamics at age
  fish =ggplot(dat)+
    geom_line(aes(year,data,group=age,col=Age))+
    facet_wrap(~qname,scale="free",ncol=2)+
    theme_bw()+theme(legend.position="right")+
    xlab("Year")+xlab("Quantities")+ ggtitle(paste0(i,".",stk@name,": Stock dynamics (",stk@species,")"))+
    theme(legend.position="right",legend.text = element_text(size=8),
          legend.title = element_text(size=9),legend.key.height = unit(.5, 'cm'))
  
  print(fish)
}
dev.off()
 

