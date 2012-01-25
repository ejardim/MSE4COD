library(MASS)
library(FLCore)
library(FLash)
#library(FLBRP)

my.dir<-"C:/current/Meetings/2011/stecf_jun/woscod/"

# Activate the history function for the plots
#windows(record=T)
 
source(paste(my.dir,"/R/mseFuncsb.R",sep=""))
source(paste(my.dir,"/R/plotFuncs.R",sep=""))
brun<-list()
frun<-list()
load(paste(my.dir,"data/biolfleetWOScod_1-6.RData",sep=""))
for (i in c(1:6)) {brun[[i]]<-biolrun[[i]];frun[[i]]<-fleetrun[[i]]}
rm(biolrun,fleetrun)
load(paste(my.dir,"data/biolfleetWOScod_7-14.RData",sep=""))
for (i in c(7:14)) {brun[[i]]<-biolrun[[i]];frun[[i]]<-fleetrun[[i]]}
rm(biolrun,fleetrun)

nits<-dims(brun[[i]])$iter
minTAC<-0.8
maxTAC<-1.2
runs <- rbind(expand.grid(OEM=c("wg","m","catch"),OM=c("m","catch"),SR=c(1,0.5),HCR=c("EU"),fcut=c(0.75),minTAC=0.8,maxTAC=1.2),
              expand.grid(OEM=c("catch"),OM=c("catch"),SR=c(1,0.5),HCR=c("EU"),fcut=c(0.75),minTAC=0.001,maxTAC=1000))

for (i in c(1:14)){
  OM     <-as.character(runs[i,"OM", drop=T])
  OEM    <-as.character(runs[i,"OEM",drop=T])
  HCR    <-as.character(runs[i,"HCR",drop=T])
  SR     <-runs[i,"SR", drop=T]
  mnTAC <-runs[i,"minTAC", drop=T]
  mxTAC <-runs[i,"maxTAC", drop=T]
  fcut <-runs[i,"fcut", drop=T]
  if (i == 1) cod4Smry<-cbind(OEM=OEM,OM=OM,HCR=HCR,SR=SR,fcut=fcut,model.frame(smryStats(frun[[i]],brun[[i]],start=2006,end=2030))[,-c(1,3:5)],minTAC=mnTAC,maxTAC=mxTAC)
  else {
    dfOM<-cbind(OEM=OEM,OM=OM,HCR=HCR,SR=SR,fcut=fcut,model.frame(smryStats(frun[[i]],brun[[i]],start=2006,end=2030))[,-c(1,3:5)],minTAC=mnTAC,maxTAC=mxTAC)
    cod4Smry<-rbind(cod4Smry,dfOM)
  }
}
rm(OM,OEM,HCR,SR,fcut,mnTAC,mxTAC,dfOM)
 
Bpa      <-22000
Blim     <-14000
Fmsylo   <-0.17
Fmsy     <-0.19
Fmsyhi   <-0.33
 
probbio <- function(cod,y,brp){
  t.<-tapply(ifelse(cod[cod[,"year"]==y & cod[,"minTAC"]==minTAC & cod[,"maxTAC"]==maxTAC,"ssb"]<brp,0,1),cod[cod[,"year"]==y & cod[,"minTAC"]==minTAC & cod[,"maxTAC"]==maxTAC,c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC")],mean)
  t..<-tapply(ifelse(cod[cod[,"year"]==y & cod[,"minTAC"]==0.001 & cod[,"maxTAC"]==1000.,"ssb"]<brp,0,1),cod[cod[,"year"]==y & cod[,"minTAC"]==0.001 & cod[,"maxTAC"]==1000.,c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC")],mean)
  t1. <- na.omit(cbind(expand.grid(dimnames(t.)),c(t.)))
  t1.. <- na.omit(cbind(expand.grid(dimnames(t..)),c(t..)))
  colnames(t1..)[8]<-"c(t.)"
  t1.<-rbind(t1.,t1..)
  return(t1.)
}
probf <- function(cod,y,frp){
  t.<-tapply(ifelse(cod[cod[,"year"]==y & cod[,"minTAC"]==minTAC & cod[,"maxTAC"]==maxTAC,"fCbar"]>frp,0,1),cod[cod[,"year"]==y & cod[,"minTAC"]==minTAC & cod[,"maxTAC"]==maxTAC,c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC")],mean)
  t..<-tapply(ifelse(cod[cod[,"year"]==y & cod[,"minTAC"]==0.001 & cod[,"maxTAC"]==1000.,"fCbar"]>frp,0,1),cod[cod[,"year"]==y & cod[,"minTAC"]==0.001 & cod[,"maxTAC"]==1000.,c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC")],mean)
  t1. <- na.omit(cbind(expand.grid(dimnames(t.)),c(t.)))
  t1.. <- na.omit(cbind(expand.grid(dimnames(t..)),c(t..)))
  colnames(t1..)[8]<-"c(t.)"
  t1.<-rbind(t1.,t1..)
  return(t1.)
}
stats <- function(cod,y,stat,typ,div){
  t.<-tapply(cod[cod[,"year"]==y & cod[,"minTAC"]==minTAC & cod[,"maxTAC"]==maxTAC,stat],cod[cod[,"year"]==y & cod[,"minTAC"]==minTAC & cod[,"maxTAC"]==maxTAC,c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC")],typ)/div
  t..<-tapply(cod[cod[,"year"]==y & cod[,"minTAC"]==0.001 & cod[,"maxTAC"]==1000.,stat],cod[cod[,"year"]==y & cod[,"minTAC"]==0.001 & cod[,"maxTAC"]==1000.,c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC")],typ)/div
  t1. <- na.omit(cbind(expand.grid(dimnames(t.)),c(t.)))
  t1.. <- na.omit(cbind(expand.grid(dimnames(t..)),c(t..)))
  colnames(t1..)[8]<-"c(t.)"
  t1.<-rbind(t1.,t1..)
  return(t1.)
}

##Probabaility >Blim in yrange[2:end]
#yrange<-c(2008,2010,2012,2015)
yrange<-c(2011:2020)
t1.<-probbio(cod4Smry,yrange[1],Blim)
t2.<-probbio(cod4Smry,yrange[2],Blim)
t2. <- merge(t1.,t2.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
colnames(t2.)[(length(colnames(t2.))-1):length(colnames(t2.))] <- c(paste("p(>Blim)",yrange[1],sep=""),paste("p(>Blim)",yrange[2],sep=""))
for (i in c(yrange[3:length(yrange)])) {
  t1.<-probbio(cod4Smry,i,Blim)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("p(>Blim)",i,sep=""))
}
##Probabaility >Bpa in yrange[2:end]
for (i in c(yrange[1:length(yrange)])) {
  t1.<-probbio(cod4Smry,i,Bpa)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("p(>Bpa)",i,sep=""))
}
##Probabaility <Fmsylo in yrange[2:end]
for (i in c(yrange[1:length(yrange)])) {
  t1.<-probf(cod4Smry,i,Fmsylo)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("p(<Fmsylo)",i,sep=""))
}
##Probabaility <Fmsy in yrange[2:end]
for (i in c(yrange[1:length(yrange)])) {
  t1.<-probf(cod4Smry,i,Fmsy)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("p(<Fmsy)",i,sep=""))
}
##Probabaility <Fmsyhi in yrange[2:end]
for (i in c(yrange[1:length(yrange)])) {
  t1.<-probf(cod4Smry,i,Fmsyhi)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("p(<Fmsyhi)",i,sep=""))
}
##Landings Yield in yrange (avg)
for (i in yrange) {
  t1.<-stats(cod4Smry,i,"ssb","median",1000.)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("med(SSB)",i,sep=""))
}
##Landings Yield in yrange
for (i in yrange) {
  t1.<-stats(cod4Smry,i,"yield","median",1000.)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("med(L)",i,sep=""))
}
##Discards in yrange
for (i in yrange) {
  t1.<-stats(cod4Smry,i,"discards","median",1000.)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("med(D)",i,sep=""))
}
##Catch in yrange
for (i in yrange) {
  t1.<-stats(cod4Smry,i,"catch","median",1000.)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("med(C)",i,sep=""))
}
##Harvest in yrange
for (i in yrange) {
  t1.<-stats(cod4Smry,i,"harvest","median",1.)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("med(H)",i,sep=""))
}
##Landings F in yrange
for (i in yrange) {
  t1.<-stats(cod4Smry,i,"fLbar","median",1.)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("med(FL)",i,sep=""))
}
##Discards F in yrange
for (i in yrange) {
  t1.<-stats(cod4Smry,i,"fDbar","median",1.)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("med(FD)",i,sep=""))
}
##Catch F in yrange
for (i in yrange) {
  t1.<-stats(cod4Smry,i,"fCbar","median",1.)
  t2. <- merge(t2.,t1.,by=c("OEM","SR","OM","HCR","fcut","minTAC","maxTAC"))
  colnames(t2.)[length(colnames(t2.))] <- c(paste("med(FC)",i,sep=""))
}
##Re-calibrate F's using Joe's method
for (i in yrange) {
  t2.[,paste("med(FL)",i,sep="")]<-t2.[,paste("med(FC)",i,sep="")]*t2.[,paste("med(FL)",i,sep="")]/(t2.[,paste("med(FL)",i,sep="")]+t2.[,paste("med(FD)",i,sep="")])
  t2.[,paste("med(FD)",i,sep="")]<-t2.[,paste("med(FC)",i,sep="")]-t2.[,paste("med(FL)",i,sep="")]
}
#=============================================================================
print(t2.)
write.csv(t2.,file=paste(my.dir,"/tables/smry_WoScod_1-14.csv",sep=""))

