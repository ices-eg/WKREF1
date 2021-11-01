
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# R script to build short-mse for WKREF1
# @author Henning Winker
# @email henning.winker@ec.europa.eu
# Licence: EUPL
# mse version: https://github.com/flr/mse/tree/FLombf
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


library(FLCore)
library(ggplotFL)
library(mse)
library(FLBRP)
library(FLSRTMB)
library(gridExtra)

load(file="data/ices.stks.n78.Rdata",verbose=T)

do.call(rbind,lapply(stks,function(x) x@name)) ==stks@names 


names(stks)
# sort
stks = stks[sort(stks@names)]


x= stks[[52]]
x@fishlife
x@benchmark
x@species
x@eqsim
x@desc

# Exclude stocks with no Fmsy benchmark
excl = do.call(rbind,lapply(stks,function(x) is.na(x@benchmark[["Fmsy"]])))
stks = stks[-which(excl)]
n = length(stks)

#-----------------------------------------------------------------
# Plot spr0 vs B/B0 for ICES system
#-----------------------------------------------------------------

d. = do.call(rbind,lapply(stks,function(x){ 
data.frame(species=x@species,stock=x@name, as.data.frame(spr0y(x)/mean(spr0y(x))),bb0=c(ssb(x)/x@eqsim[["B0"]]))
}))

d.$spr0 = d.$data

d.= d.[d.$bb0<=1,]

p1= ggplot(data=d.,aes(bb0,spr0,color=stock))+geom_point()+
  geom_hline(yintercept = 1,size=0.5,linetype="dotted")+
  ylab("spr0")+xlab("B/B0")+theme(legend.position = "none")+
  facet_wrap(~stock, scales="free", ncol=8)
p2 <- p1+ geom_smooth(data =d., aes(bb0,spr0),col=1,se=F, method="lm",show.legend = FALSE)


#-----------------------------------------------------------------
# 3:  Get the SRR
#-----------------------------------------------------------------

# beverton holt (mean spr0y for a,b)
bh.om <- FLSRs(lapply(stks, function(x) { 
  srr = as.FLSR(x,model=bevholtSV)
  
  return(srrTMB(srr,s.est=T,s=x@fishlife[["s"]],
                s.logitsd =x@fishlife[["sd.logit.s"]],spr0=spr0y(x),nyears=NULL,SDreport = F))
}))   


# segmented regression
sr.ices <- FLSRs(lapply(stks, function(x) { 
  return(fmle(as.FLSR(x,model=segreg),fixed=list(b=x@benchmark[["Blim"]])))
}))   

# beverton holt (spr0 for last 3 years)
bh.mp <- FLSRs(lapply(stks, function(x) { 
  srr = as.FLSR(x,model=bevholtSV)
  
  return(srrTMB(srr,s.est=T,spr0=spr0y(x),SDreport = F,nyears=3))
}))   





# for E40 and f01 (spr0 for last 3 years)
sr.b01 <- FLSRs(lapply(stks, function(x) { 
  srr = as.FLSR(x,model=segreg)
  
  return(srrTMB(srr,s.est=T,plim=0.1,spr0=spr0y(x),nyears = 3,SDreport = F))
}))   



srs = list(bh.om=bh.om,bh.mp=bh.mp,sr.ices=sr.ices,sr.b01=sr.b01)

p = list()
bh.scale = NULL
for(i in 1:length(bh.om)){
  p[[i]]=plot(FLSRs(bh.om=bh.om[[i]],bh.mp=bh.mp[[i]],hs.ices=sr.ices[[i]],sr.f01=sr.b01[[i]]))+
    theme(plot.title = element_text(color="blue", size=8),
          legend.position = "right",legend.key.size = unit(0.1, 'cm'),legend.text = element_text(size=7))+xlab("SSB")+ylab("Recruits")+  
    ggtitle(paste0(bh.om[[i]]@name,", s.om=",round(bh.om[[i]]@SV$s,2),", s.mp=",round(bh.mp[[i]]@SV$s,2)))
}
names(p) = stks@names

# plot frames
sel1 = 1:10
sel2 = sel1+10
sel3 = sel2+10
sel4 = sel3+10
sel5 = sel4+10
sel6 = sel5+10
sel7 = (sel6+10)[1:3]


# Check SRR
pdf("report/buildOM_srr.pdf")
grid.arrange(grobs =  p[sel1], ncol = 2)
grid.arrange(grobs =  p[sel2], ncol = 2)
grid.arrange(grobs =  p[sel3], ncol = 2)
grid.arrange(grobs =  p[sel4], ncol = 2)
grid.arrange(grobs =  p[sel5], ncol = 2)
grid.arrange(grobs =  p[sel6], ncol = 2)
grid.arrange(grobs =  p[sel7], ncol = 2)

dev.off()

#><>><>><>><>><>><>><>><>><>><>><>><>><>
# Set up BRPs
#><>><>><>><>><>><>><>><>><>><>><>><>><>

# OM
brp.om = Map(function(x,y){
  brp = brp(FLBRP(x,y))
  brp@name = y@name
  attr(brp,"SV") = y@SV 
  brp
},  stks,bh.om)


do.call(rbind,lapply(stks,function(x){
  range(x)
}))

# BevHolt MP
brp.bh = Map(function(x,y){
  brp = brp(FLBRP(x,y))
  brp = brp
  brp@name = y@name
  attr(brp,"SV") = y@SV 
  brp
},  stks,bh.mp)

# ICES segreg MP, using ICES refpts
brp.ices = Map(function(x,y){
  brp = brp(FLBRP(x,y))
  brp = brp+FLPar(Fmsy =x@benchmark[["Fmsy"]],Blim=x@benchmark[["Blim"]],Bpa=x@benchmark[["Bpa"]],Btrigger=x@benchmark[["Btrigger"]],Bmsy=x@eqsim[["BMSY"]],B0=x@eqsim[["B0"]],CMSY=x@eqsim[["Catchequi"]])
  brp@name = y@name
  brp
},  stks,sr.ices)

# SB40 MP
brp.b40 = Map(function(x,y){
  brp = brp(FLBRP(x,y))
  B0 = an(refpts(brp)["virgin","ssb"])
  Bmsy = 0.4*B0
  brp = brp+ as(data.frame(Bmsy=0.4*B0,Blim=0.2*B0,B0=B0),"FLPar")
  brp@name = y@name
  attr(brp,"SV") = y@SV 
  brp
},  stks,bh.mp)

# MP F0.1
brp.f01 = Map(function(x,y){
  brp = brp(FLBRP(x,y))
  B0 = an(refpts(brp)["virgin","ssb"])
  Fmsy = an(refpts(brp)["f0.1","harvest"])
  brp = brp+ as(data.frame(Fmsy=Fmsy),"FLPar")
  brp@name = y@name
  attr(brp,"SV") = y@SV 
  brp
},  stks,sr.b01)


# spr40 MP
brp.spr40 = lapply(stks,function(x){
  brp = brp(FLBRP(x))
  Bcur = brp@ssb.obs[,length(ssb.obs(brp))]/mean(brp@rec.obs[,(length(brp@rec.obs)-2):length(brp@rec.obs)])
  B0 = an(refpts(brp)["virgin","ssb"])
  brp = brp+FLPar(Bmsy=an(0.4*B0),Bcur=Bcur)
  brp@name = x@name
  brp
})


brps = list()
for(i in 1:length(brp.om)){
brps[[i]] = FLBRPs(list(om=brp.om[[i]],
                        ices=brp.ices[[i]],
                        bh=brp.bh[[i]],
                        b40=brp.b40[[i]],
                        f01=brp.f01[[i]],
                        spr40=brp.spr40[[i]]))
}
names(brps) = stks@names


# Check SRR + SP by stock
pdf("report/buildOM_srr_sp_stock.pdf")

for(s in 1:length(stks)){
  
  sr = (FLSRs(lapply(srs,function(x){
    x[[s]]
  })))
  sr = sr[c(1,3,2,4)]
  
  sr@names = c("1.om","2.ices","3.bevholt","4.hs.f01")
  
  psr = plot(sr)+ylab("Recruits")+xlab("SSB")+
    theme(legend.position = "right")+ggtitle(paste0(s,". ",brps[[s]][[1]]@name,": SRRs, OM s = ",round(sr[[1]]@SV[["s"]],2)))
  
  eqref=do.call(rbind,Map(function(x,y){
    data.frame(model.frame(metrics(x,list(ssb=ssb, harvest=fbar, rec=rec, yield=landings)),drop=FALSE),type=y)
  },brps[[s]][c(1,2,3,5)],as.list(sr@names)))
  
  sp=ggplot(data=eqref,aes(x=ssb,y=yield,group=type))+
    geom_line(aes(color=type))+ggtitle(paste0(s,". ",brps[[s]][[1]]@name,": Production Functions"))
  
  grid.arrange(psr,sp)
}
dev.off()

#--------------------------------------
# Reference Points
#--------------------------------------

br.om= do.call(rbind,lapply(brp.om,function(x){
  all = model.frame(refpts(x))
  data.frame(stock = x@name,
             Blim = all[all$quant=="ssb","msy"]/2,  # Blim = Bmsy/2
             Bpa = all[all$quant=="ssb","msy"]*0.7, # place holder
             Btrigger = all[all$quant=="ssb","msy"]*0.8, 
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


###><> Changed to direct input 
br.ices= do.call(rbind,lapply(stks,function(x){
  data.frame(stock = x@name,
             Blim = an(x@benchmark[["Blim"]]),  # Blim = Bmsy/2
             Bpa = x@benchmark[["Bpa"]],
             Btrigger = x@benchmark[["Btrigger"]], 
             Bmsy = x@eqsim[["BMSY"]],
             B0 = x@eqsim[["B0"]],
             Fmsy = x@benchmark[["Fmsy"]],
             F01 = NA,
             MSY = x@eqsim[["Catchequi"]],
             BBmsy = an(tail(ssb(x)))/x@eqsim[["BMSY"]],
             FFmsy = an(tail(fbar(x)))/x@benchmark[["Fmsy"]],
             BB0 = an(tail(ssb(x)))/x@eqsim[["B0"]],
             type = "ices")
}    
))

br.bh= do.call(rbind,lapply(brp.bh,function(x){
  all = model.frame(refpts(x))
  data.frame(stock = x@name,
             Blim = all[all$quant=="ssb","msy"]/2,  # Blim = Bmsy/2
             Bpa = all[all$quant=="ssb","msy"]*0.7, # place holder
             Btrigger = all[all$quant=="ssb","msy"]*0.8, 
             Bmsy = all[all$quant=="ssb","msy"],
             B0 = all[all$quant=="ssb","virgin"],
             Fmsy = all[all$quant=="harvest","msy"],
             F01 = all[all$quant=="harvest","f0.1"],
             MSY = all[all$quant=="yield","msy"],
             BBmsy = an(x@ssb.obs[,length(x@ssb.obs)])/all[all$quant=="ssb","msy"],
             FFmsy = an(x@fbar.obs[,length(x@fbar.obs)])/all[all$quant=="harvest","msy"],
             BB0 = an(x@ssb.obs[,length(x@ssb.obs)])/all[all$quant=="ssb","virgin"],
             type = "bh")
}    
))

br.b40= do.call(rbind,lapply(brp.b40,function(x){
  all = model.frame(refpts(x))
  data.frame(stock = x@name,
             Blim = all[all$quant=="ssb","Bmsy"]/2,  # Blim = Bmsy/2
             Bpa = all[all$quant=="ssb","Bmsy"]*0.7, # place holder
             Btrigger = all[all$quant=="ssb","Bmsy"]*0.8, 
             Bmsy = all[all$quant=="ssb","Bmsy"],
             B0 = all[all$quant=="ssb","virgin"],
             Fmsy = all[all$quant=="harvest","Bmsy"],
             F01 = all[all$quant=="harvest","f0.1"],
             MSY = all[all$quant=="yield","Bmsy"],
             BBmsy = an(x@ssb.obs[,length(x@ssb.obs)])/all[all$quant=="ssb","Bmsy"],
             FFmsy = an(x@fbar.obs[,length(x@fbar.obs)])/all[all$quant=="harvest","Bmsy"],
             BB0 = an(x@ssb.obs[,length(x@ssb.obs)])/all[all$quant=="ssb","virgin"],
             type = "b40")
}    
))

# Revised with Btrigger based on Blim and sigma = 0.3
br.f01= do.call(rbind,lapply(brp.f01,function(x){
  all = model.frame(refpts(x))
  data.frame(stock = x@name,
             Blim =sr.b01[[x@name]]@params[[2]],  # Blim = b of sr.b01
             Bpa = sr.b01[[x@name]]@params[[2]]*exp(1.645*0.3), # Bpa 
             Btrigger = sr.b01[[x@name]]@params[[2]]*exp(1.645*0.3), # Trigger with sigma = 0.3 
             Bmsy = all[all$quant=="ssb","f0.1"],
             B0 = all[all$quant=="ssb","virgin"],
             Fmsy = all[all$quant=="harvest","f0.1"],
             F01 = all[all$quant=="harvest","f0.1"],
             MSY = all[all$quant=="yield","f0.1"],
             BBmsy = an(x@ssb.obs[,length(x@ssb.obs)])/all[all$quant=="ssb","f0.1"],
             FFmsy = an(x@fbar.obs[,length(x@fbar.obs)])/all[all$quant=="harvest","f0.1"],
             BB0 = an(x@ssb.obs[,length(x@ssb.obs)])/all[all$quant=="ssb","virgin"],
             type = "f01")
}    
))

br.spr40= do.call(rbind,lapply(brp.spr40,function(x){
  all = model.frame(refpts(x))
  data.frame(stock = x@name,
             Blim = all[all$quant=="ssb","Bmsy"]/2,  # Blim = Bmsy/2
             Bpa = all[all$quant=="ssb","Bmsy"]*0.7, # place holder
             Btrigger = all[all$quant=="ssb","Bmsy"]*0.8, 
             Bmsy = all[all$quant=="ssb","Bmsy"],
             B0 = all[all$quant=="ssb","virgin"],
             Fmsy = all[all$quant=="harvest","Bmsy"],
             F01 = all[all$quant=="harvest","f0.1"],
             MSY = all[all$quant=="yield","Bmsy"],
             BBmsy = all[all$quant=="ssb","Bcur"]/all[all$quant=="ssb","Bmsy"],
             FFmsy = an(x@fbar.obs[,length(x@fbar.obs)])/all[all$quant=="harvest","Bmsy"],
             BB0 = all[all$quant=="ssb","Bcur"]/all[all$quant=="ssb","virgin"],
             type = "spr40")
}    
))


mp.refs = rbind(br.om,br.ices,br.bh,br.b40,br.f01,br.spr40)
refs = mp.refs[mp.refs$type=="om",]

# Compare
# Plot
status= mp.refs
for(i in 1:length(unique(status$type))){
  status[status$type%in%unique(status$type)[i],"type"] = paste0(i,".",unique(status$type)[i])
}

pdf("report/comp_status.pdf")
xylim = c(3.5,8.5)
gplot(status, aes(x = BBmsy, y = FFmsy,color=stock))+
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = 0, ymax = 1), 
            colour = "green", fill = "green") + 
  geom_rect(aes(xmin = 0,xmax = 1, ymin = 0, ymax = 1), colour = "yellow", fill = "yellow") +
  geom_rect(aes(xmin = 1, xmax = Inf,ymin = 1, ymax = Inf), colour = "orange", fill = "orange") + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 1, ymax = Inf), colour = "red", fill = "red") + 
  geom_hline(aes(yintercept = 1))+geom_vline(aes(xintercept = 1)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, xylim[2]))+ 
  scale_x_continuous(expand = c(0,0), limits = c(0, xylim[1]))+
  geom_point(aes(fill = stock),size=2,col=1,pch=21)+ylab(expression(F/F[MSY]))+
  xlab(expression(B/B[MSY]))+ggtitle("Stock Status: 6 reference point systems")+facet_wrap(~type,scale="free",ncol=2)+
  theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=6), #change legend title font size
        legend.text = element_text(size=6))+ #change legend text font size
 guides(fill=guide_legend(ncol=2,byrow=TRUE))

dat = mp.refs[mp.refs$type=="ices",]
xylim = c(0.7)
ggplot(dat, aes(x = Blim/B0, y = Bpa/B0,color=stock))+
  geom_abline(linetype="dashed")+
  # crit
  geom_segment(aes(x = 0, xend = 0.1, y = 0.1, yend = 0.1), colour = "red",size=1)+
  geom_segment(aes(x = 0.1, xend = 0.1, y = 0, yend = 0.1), colour = "red",size=1)+
  geom_segment(aes(x = 0.1, xend = 0.1, y = 0, yend = 0.4), colour = "red",size=0.4,linetype="dashed")+
  # 20%
  geom_segment(aes(x = 0, xend = 0.2, y = 0.2, yend = 0.2), colour = "orange",size=1)+
  geom_segment(aes(x = 0.2, xend = 0.2, y = 0, yend = 0.2), colour = "orange",size=1)+
  geom_segment(aes(x = 0.2, xend = 0.2, y = 0, yend = 0.4), colour = "orange",size=0.4,linetype="dashed")+
  
  # 30%
  geom_segment(aes(x = 0, xend = 0.3, y = 0.3, yend = 0.3), colour = "yellow",size=1)+
  geom_segment(aes(x = 0.3, xend = 0.3, y = 0, yend = 0.3), colour = "yellow",size=1)+
  # 40%
  geom_segment(aes(x = 0, xend = 0.4, y = 0.4, yend = 0.4), colour = "green",size=1)+
  geom_segment(aes(x = 0.4, xend = 0.4, y = 0, yend = 0.4), colour = "green",size=1)+
  geom_point(aes(fill = stock),size=2,col=1,pch=21)+
  scale_y_continuous(expand = c(0,0), limits = c(0, xylim[1]))+
  scale_x_continuous(expand = c(0,0), limits = c(0, xylim[1]))+
  ggtitle("ICES Benchmarks: relative Blim and Bpa")+
  ylab(expression(B[pa]/B[0]))+xlab(expression(B[lim]/B[0]))+
  theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=6), #change legend title font size
        legend.text = element_text(size=6),
        legend.position = "bottom")+ #change legend text font size
  guides(fill=guide_legend(nrow=9,byrow=TRUE))
  


dat = mp.refs[mp.refs$type=="ices",]
xylim = c(0.7)
ggplot(dat, aes(x = Btrigger/B0, y = Bmsy/B0,color=stock))+
  geom_abline(linetype="dashed")+
  # crit
  geom_segment(aes(x = 0, xend = 0.1, y = 0.1, yend = 0.1), colour = "red",size=1)+
  geom_segment(aes(x = 0.1, xend = 0.1, y = 0, yend = 0.1), colour = "red",size=1)+
  geom_segment(aes(x = 0.1, xend = 0.1, y = 0, yend = 0.4), colour = "red",size=0.4,linetype="dashed")+
  # 20%
  geom_segment(aes(x = 0, xend = 0.2, y = 0.2, yend = 0.2), colour = "orange",size=1)+
  geom_segment(aes(x = 0.2, xend = 0.2, y = 0, yend = 0.2), colour = "orange",size=1)+
  geom_segment(aes(x = 0.2, xend = 0.2, y = 0, yend = 0.4), colour = "orange",size=0.4,linetype="dashed")+
  
  # 30%
  geom_segment(aes(x = 0, xend = 0.3, y = 0.3, yend = 0.3), colour = "yellow",size=1)+
  geom_segment(aes(x = 0.3, xend = 0.3, y = 0, yend = 0.3), colour = "yellow",size=1)+
  # 40%
  geom_segment(aes(x = 0, xend = 0.4, y = 0.4, yend = 0.4), colour = "green",size=1)+
  geom_segment(aes(x = 0.4, xend = 0.4, y = 0, yend = 0.4), colour = "green",size=1)+
  
  geom_point(aes(fill = stock),size=2,col=1,pch=21)+
  scale_y_continuous(expand = c(0,0), limits = c(0, xylim[1]))+
  scale_x_continuous(expand = c(0,0), limits = c(0, xylim[1]))+
  ggtitle("ICES Benchmarks: relative Btrigger and Bmsy")+
  ylab(expression(B[msy]/B[0]))+xlab(expression(MSYB[trigger]/B[0]))+
  theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=6), #change legend title font size
        legend.text = element_text(size=6),
        legend.position = "bottom")+ #change legend text font size
  guides(fill=guide_legend(nrow=9,byrow=TRUE))


dev.off()

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# SET OM structure and HCR
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


#-------------------------------------------------------------------------------
# 1 : setting dimensions
#-------------------------------------------------------------------------------

it <- 250   # number of iterations                       
ny <- 60    # number of years to project 
nsqy <- 3 # 

# Assign attributes of mpargs to all OMs
stks = FLStocks(lapply(stks,function(x){
  mpargs = list(
    y0=an(range(x)["minyear"]),
    dy=an(range(x)["maxyear"])-1,
    ay=an(range(x)["maxyear"])-1,
    iy=an(range(x)["maxyear"]),
    fy=an(range(x)["maxyear"])+ny,
    it=it,ny=ny,nsqy=nsqy)
  attr(x,"args") = mpargs
  x
}    
))

# Project...takes a second
foms <- FLStocks(lapply(stks,function(x){
  stf(x, ny, nsqy, nsqy)
})) 

# propagate
foms <- FLStocks(lapply(foms,function(x){
  propagate(x,it)
})) 

#-----------------------------------------------------------
#  2. Build OM rec functions
#----------------------------------------------------------

# Bias corrected ar1 after Thorson 2020
ar1rlnorm <- function(rho, years, iters=1, mean=0, margSD=0.6,bias.correct=TRUE) {

  logbias=0
  if(bias.correct)logbias = 0.5*margSD^2
  n <- length(years)
  rhosq <- rho ^ 2
  #
  res <- matrix(rnorm(n*iters, mean=mean, sd=margSD), nrow=n, ncol=iters)
  res <- apply(res, 2, function(x) {
    for(i in 2:n)
      x[i] <- sqrt(rhosq) * x[i-1] + sqrt(1-rhosq) * x[i]
    return(exp(x-logbias))
  })
    return(FLQuant(array(res, dim=c(1,n,1,1,1,iters)),
                 dimnames=list(year=years, iter=seq(1, iters))))
}


# Beverton Holt
set.seed(123)
om = Map(function(x,y){
  res = window(rec(x),end=x@args$fy+1)   
  res = ar1rlnorm(x@fishlife[["rho"]],years=ac(range(x)["minyear"]:range(x)["maxyear"]),it=x@args$it,margSD = x@fishlife[["sigmaR"]])
  res_om = propagate(y,x@args$it)
  residuals(res_om) <- res
  res_om
},foms,bh.om)



#-----------------------------------------------------------
# 3. Build OMs functions
#----------------------------------------------------------
# set projection method for OM
proj <- mseCtrl(method=fwd.om, args=list(maxF=3.))

om =  Map(function(x,y,z){
  FLom(stock=x, sr=y, refpts=refpts(z), projection=proj)  }
  ,foms,om,brp.om)


# Check that args are still there
om[[1]]@stock@args
om[[1]]@stock
om[[1]]@sr
om[[1]]@refpts
om[[1]]@stock@args # mse agrs


#-------------------------------------------------------
# 4. harvest control rule function 
#------------------------------------------------------



# F-HCR will be converted to TACy+1 by mse::tac.is
ICEShcr = function (stk, ftrg, blim, btrigger, fmin = 0.001,as.spr=FALSE, args, tracking) 
{
  ay <- args$ay
  data_lag <- args$data_lag
  man_lag <- args$management_lag
  ssb <- ssb(stk)[, ac(ay - data_lag)]
  if(as.spr){ # Translate biomass into SB/R, R recruitment is geomean of last 3 years
    ssb <-  ssb(stk)[, ac(ay - data_lag)]/exp(mean(log(stock.n(stk)[1,ac((ay - data_lag-2):(ay - data_lag))])))
  }
  fout <- FLQuant(fmin, dimnames = list(iter = dimnames(ssb)$iter))
  fout[ssb >= btrigger] <- ftrg # Fish at Ftrg if > Btrigger
  inbetween <- (ssb < btrigger) & (ssb > blim) # here Blim at origin
  gradient <- (ftrg - fmin)/(btrigger - blim)
  fout[inbetween] <- (ssb[inbetween] - blim) * gradient + fmin
  ctrl <- fwdControl(year = ay + man_lag, quant = "fbar", 
                     value = c(fout))
  list(ctrl = ctrl, tracking = tracking)
}


#-------------------------------------------------------
# 5. Parameterize Harvest Control Rule for MPs
#------------------------------------------------------

# OM hcr (Self-Test)
ref = mp.refs[mp.refs$type=="om",]
y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
hcr.om = lapply(y,function(x){
  mseCtrl(method=ICEShcr,args=list(fmin=0.001,ftrg=x$Fmsy,blim=0.001,btrigger=x$Btrigger))
}) 

# ICES hcr 
ref = mp.refs[mp.refs$type=="ices",]
y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
hcr.ices = lapply(y,function(x){
  mseCtrl(method=ICEShcr,args=list(fmin=0.001,ftrg=x$Fmsy,blim=0.001,btrigger=x$Btrigger))
}) 

# BevHolt hcr
ref = mp.refs[mp.refs$type=="bh",]
y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
hcr.bh = lapply(y,function(x){
  mseCtrl(method=ICEShcr,args=list(fmin=0.001,ftrg=x$Fmsy,blim=0.001,btrigger=x$Btrigger))
}) 


# BevHolt SB40 hcr
ref = mp.refs[mp.refs$type=="b40",]
y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
hcr.b40 = lapply(y,function(x){
  mseCtrl(method=ICEShcr,args=list(fmin=0.001,ftrg=x$Fmsy,blim=0.001,btrigger=x$Btrigger))
}) 


# hs-f0.1 hcr 
ref = mp.refs[mp.refs$type=="f01",]
y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
hcr.f01 = lapply(y,function(x){
  mseCtrl(method=ICEShcr,args=list(fmin=0.001,ftrg=x$Fmsy,blim=0.001,btrigger=x$Btrigger))
}) 


# spr40 hcr
ref = mp.refs[mp.refs$type=="spr40",]
y= setNames(split(ref, seq(nrow(ref))), rownames(ref))
hcr.spr40 = lapply(y,function(x){
  mseCtrl(method=ICEShcr,args=list(fmin=0.001,ftrg=x$Fmsy,blim=0.001,btrigger=x$Btrigger,as.spr=TRUE))
}) 


#-------------------------------------------------------
# 6. set of MP short-cut ctrl
#------------------------------------------------------
# om mp (self-test)
ctrl.om = lapply(hcr.om,function(x){
  mpCtrl(list(
    est = mseCtrl(method=perfect.sa),
    hcr  = x,
    iysy = mseCtrl(method=tac.is)
  ))})

# ices mp
ctrl.ices = lapply(hcr.ices,function(x){
  mpCtrl(list(
    est = mseCtrl(method=perfect.sa),
    hcr  = x,
    iysy = mseCtrl(method=tac.is)
  ))})

ctrl.bh = lapply(hcr.bh,function(x){
  mpCtrl(list(
    est = mseCtrl(method=perfect.sa),
    hcr  = x,
    iysy = mseCtrl(method=tac.is)
  ))})

ctrl.b40 = lapply(hcr.b40,function(x){
  mpCtrl(list(
    est = mseCtrl(method=perfect.sa),
    hcr  = x,
    iysy = mseCtrl(method=tac.is)
  ))})

ctrl.f01 = lapply(hcr.f01,function(x){
  mpCtrl(list(
    est = mseCtrl(method=perfect.sa),
    hcr  = x,
    iysy = mseCtrl(method=tac.is)
  ))})

ctrl.spr40 = lapply(hcr.spr40,function(x){
  mpCtrl(list(
    est = mseCtrl(method=perfect.sa),
    hcr  = x,
    iysy = mseCtrl(method=tac.is)
  ))})


# Short cut control with 6 HCRs
sc.ctrl = list(om=ctrl.om,ices=ctrl.ices,bh=ctrl.bh,b40=ctrl.b40,f01=ctrl.f01,spr40=ctrl.spr40)

#----------------------------------------------------------------------------------------------
# Ready to gofish()
#---------------------------------------------------------------------------------------------

hcrs = c("om","ices","bevholt","sb40","f01","spr40")
hcr.names = hcrs

# Save
save(om,stks,sc.ctrl,mp.refs,hcr.names,srs,brps,brp.om,file="model/om5mps_wkref1.Rdata")
