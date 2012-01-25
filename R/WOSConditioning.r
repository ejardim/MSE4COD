#Use R-version 2.7.2
library(FLCore)
library(FLash)
library(FLBRP)
library(FLAssess)

#source("C:/FLR/packages/FLBRP/R/FLBRP-class.R")
#source("C:/FLR/packages/FLBRP/R/FLBRP-methods.R")

my.dir<-"../"
source(paste(my.dir,"R/mseFuncs.R",sep=""))

#### read in WG data
cod4.2010 <- readFLStock(paste(my.dir,"input/cod6a.IDX",sep=""))
catch.n(cod4.2010) <- readVPAFile(paste(my.dir,"input/cod6aCN.DAT",sep=""), quiet=FALSE)
catch.wt(cod4.2010) <- readVPAFile(paste(my.dir,"input/cod6aCW.DAT",sep=""), quiet=FALSE)
discards.n(cod4.2010) <- readVPAFile(paste(my.dir,"input/cod6aDN.DAT",sep=""), quiet=FALSE)
discards.wt(cod4.2010) <- readVPAFile(paste(my.dir,"input/cod6aDW.DAT",sep=""), quiet=FALSE)

range(cod4.2010,c("minfbar","maxfbar"))<-c(2,5)
cod4.2010 <-setPlusGroup(cod4.2010 ,7, na.rm=TRUE)

#### Get bootstrap estimates
f.hat <-read.table(paste(my.dir,"Badapt used for MSE/F.CSV",sep=""),header=FALSE,sep=",") [,-c(1,35)]    ### changed MPH
n.hat <-read.table(paste(my.dir,"Badapt used for MSE/NUMBERS.CSV",sep=""),header=FALSE,sep=",") [,-c(1,35,36)] ### changed MPH

#--------------------------------------------------------------------
# using first iter for original named '0'
# dangerous, at least name iter as '1'
# better to use 2 objects, one original, one bootstrap
#--------------------------------------------------------------------

dmns <- expand.grid(age=1:7,iter=0:1000,year=1978:2010) ### changed MPH
harvest(cod4.2010) <- as.FLQuant(cbind(data=as.numeric(as.matrix(f.hat)),dmns),units="f")
stock.n(cod4.2010) <- as.FLQuant(cbind(data=as.numeric(as.matrix(n.hat)),dmns))                     ### changed MPH
rm(f.hat,n.hat)

#### Compare bootstrap and initial runs
# ToBe fixed !
jpeg(paste(my.dir,"/data/KeyRun.jpeg",sep=""))
plot(FLStocks(window(iter(cod4.2010,2:1001), end=2010), window(iter(cod4.2010,1), end=2010)), ylab="Year")
dev.off()

#### fit SR
cod4.sr <- as.FLSR(window(iter(cod4.2010,1), end=2010), model='ricker')
cod4.sr <- fmle(cod4.sr)

jpeg(file=paste(my.dir,"/data/SR.jpeg",sep=""))
plot(cod4.sr)
dev.off()

#### remove 1st run
cod4.2010 <- iter(cod4.2010,2:1001)

#### projections
## SRR
#srPar <- params(cod4.sr)
#srPar[,"b"] <-srPar[,"b"]/1000
#sr.deviates <- FLQuant(exp(rnorm(20000, sd=0.66), dimnames=list(age=1,year=2011:2030,iter=1:1000))*exp((-0.66^2)/2)
sr.deviates <- FLQuant(rlnorm(20000, sd=0.66), dimnames=list(age=1,year=2011:2030,iter=1:1000))
#save(srPar,file=paste(my.dir,"/data/condSR.RData",sep=""))

## average Recruitment
cod4 <- stf(cod4.2010,nyears=20)
trgt <- fwdControl(data.frame(year=2011:2030, quantity="f", val=0.4))
cod4 <- fwd(cod4, trgt, sr=cod4.sr, sr.residuals=sr.deviates)

## halved Recruitment
cod4.2.sr <- cod4.sr 
params(cod4.2.sr)["a"] <- params(cod4.2.sr)["a"]/2
cod4.2 <- fwd(cod4, trgt, sr=cod4.2.sr, sr.residuals=sr.deviates)

jpeg(file=paste(my.dir,"/data/ProjSRs.jpeg",sep=""))
plot(FLStocks(cod4,cod4.2))
dev.off()

#############################################################################################

#### MisReporting ratios 

#--------------------------------------------------------------------
# computeCatch.n is estimating catch.n based on F and N estimates.
# the ratios are done between this and observed catch.n loaded initially.
#--------------------------------------------------------------------

catchMisrep <- computeCatch.n(window(cod4.2010,end=2010))

jpeg(file=paste(my.dir,"/data/MisRepRatio.jpeg",sep="")) ### changed MPH - NOT 100% sure about this. Neither do we !!!
histogram( ~ data | as.character(year), data = as.data.frame(window(catchMisrep/catch.n(cod4.2010),start=1994)), xlab="Mis-reporting Ratio")
dev.off()
### Note that early years ratios are close to 1 so appears as all being 1.

#### OM scenarios
## Mis-specification of catch
cod4.Catch <- window(cod4.2010,end=2010)
catch.n(cod4.Catch) <- computeCatch.n(cod4.Catch)
landings.n(cod4.Catch) <- catch.n(cod4.Catch)*landings.n(cod4.Catch)/(discards.n(cod4.Catch)+landings.n(cod4.Catch))
discards.n(cod4.Catch) <- catch.n(cod4.Catch)-landings.n(cod4.Catch)
stock(cod4.Catch) <- computeStock(cod4.Catch)
catch(cod4.Catch) <- computeCatch(cod4.Catch, "all")
landings(cod4.Catch) <- computeLandings(cod4.Catch)
discards(cod4.Catch) <- computeDiscards(cod4.Catch)

## Mis-specification of M
#--------------------------------------------------------------------
# Hhhmmmm, you may want to take a look at this. M changes from 0.2 to ~0.7:1.2
#--------------------------------------------------------------------
cod4.M <- window(cod4.2010,end=2010)
m(cod4.M) <- sweep(harvest(cod4.M)*sweep(catchMisrep,c(1:2), catch.n(cod4.M), "-")/catchMisrep , 1:2, m(cod4.M), "+")
harvest(cod4.M) <- sweep(harvest(cod4.M)-m(cod4.M),1:2,m(cod4.Catch),"+")
units(harvest(cod4.M)) <- "f"
stock(cod4.M) <- computeStock(cod4.M)
catch(cod4.M) <- computeCatch(cod4.M,"all")
landings(cod4.M) <- computeLandings(cod4.M)
discards(cod4.M) <- computeDiscards(cod4.M)

#### projection for mean recruitment
## average Recruitment
cod4.Catch  <-stf(cod4.Catch,nyrs=20)
cod4.M      <-stf(cod4.M    ,nyrs=20)
trgt        <-fwdTarget(year=2011:2030,quantity="f",value=0.4)
cod4.Catch.1<-fwd(cod4.Catch,trgt,sr.params=srPar,sr.model="ricker",sr.residuals=sr.deviates)
cod4.M.1    <-fwd(cod4.M,    trgt,sr.params=srPar,sr.model="ricker",sr.residuals=sr.deviates)

plot(window(FLStocks(cod4.Catch.1,cod4.M.1),start=2000),xlab="Year")
savePlot(file=paste(my.dir,"/data/Proj1.jpeg",sep=""))

#### OM Misreporting/M scenarios
yrsMC        <-sample(1994:2010, 21*1000, replace=T)+rep((1:1000)*10000,each=21)             ### changed MPH
yrsReformat  <-rep(1994:2010,1000)+rep(1:1000,each=17)*10000                                 ### changed MPH

## Catch misreporting                                                                       
catchMisrepOM<-FLQuant(1,dimnames=list(age="all",year=1978:2031,iter=1:1000))                 ### changed MPH
tmp          <-(catchMisrep/catch.n(cod4.2010))[,ac(1994:2010)]                              ### changed MPH
tmp          <-FLQuant(c(tmp),dimnames=list(age=1:7,year=yrsReformat))
catchMisrepOM[,ac(2011:2031)]<-FLQuant(c(tmp[,ac(yrsMC)]),dimnames=list(age=1:7,year=2011:2031,iter=1:1000))[1,]
catchMisrepOM[,ac(1978:2010)]<-apply(catchMisrep/catch.n(cod4.2010),c(2,6),mean)                                 ### changed MPH
catchMisrepMP<-propagate(apply(catchMisrepOM,2,median),iter=1000)
save(catchMisrepOM,catchMisrepMP,file=paste(my.dir,"/data/condCatch.RData",sep=""))

## M mis-specification
tmp          <-m(cod4.M)[,ac(1994:2010)]                                             ### changed MPH
tmp          <-FLQuant(c(tmp),            dimnames=list(age=1:7,year=yrsReformat))
M_IncrOM     <-FLQuant(c(tmp[,ac(yrsMC)]),dimnames=list(age=1:7,year=2011:2031,iter=1:1000))
M_WG         <-FLQuant(c(apply(m(cod4.Catch)[,ac(1994:2010)],1,mean)),dimnames=list(age=1:7,year=1994:2031))      ### changed MPH
M_Incr       <-M_WG
M_Incr[,ac(1994:2010)]<-apply(m(cod4.M)[,ac(1994:2010)],1:2,median)                         ### changed MPH
M_Incr[,ac(2011:2031)]<-apply(M_IncrOM,1:2,median)

save(M_WG,M_Incr,catchMisrepMP,file=paste(my.dir,"/data/condM.RData",sep=""))

## intermediate year, full recruitment
trgt         <-fwdTarget(year=2011,quantity="f",value=0.75,rel=2010) #compare EC23/2010 to EC57/2011

## F2011 = F2010*.75
cod4.M       <-stf(window(cod4.M,end=2010),nyr=21)
m(cod4.M)[,ac(2011:2031)]<-M_IncrOM
cod4.M       <-fwd(cod4.M,trgt,sr.params=srPar,sr.model="ricker",sr.residuals=sr.deviates)
stock.n(cod4.M)[,ac(2012:2031)]<-NA                               #### Ugly fudge to solve stock.n being 0 instead of NA.

biol         <-as(cod4.M,"FLBiol")
fleet        <-as(cod4.M,"FLFleet")
price(fleet,1,1)[]<-c(1.05,1.05,1.40,2.01,2.759,2.759,3.39)
save(fleet,biol,file=paste(my.dir,"/data/OMM_fullR.RData",sep=""))

## F2011 = F2010*.75
cod4.Catch   <-fwd(stf(window(cod4.Catch,end=2010),nyr=21),trgt,sr.params=srPar,sr.model="ricker",sr.residuals=sr.deviates)
stock.n(cod4.Catch)[,ac(2012:2031)]<-NA                           #### Ugly fudge to solve stock.n being 0 instead of NA.

biol              <-as(cod4.Catch,"FLBiol")
fleet             <-as(cod4.Catch,"FLFleet")
price(fleet,1,1)[]<-c(1.05,1.05,1.40,2.01,2.759,2.759,3.39)
save(fleet,biol,file=paste(my.dir,"/data/OMCatch_fullR.RData",sep=""))

## intermediate year, low recruitment
trgt         <-fwdTarget(year=2011,quantity="f",value=0.75,rel=2010)

## F2011 = F2010*.75
cod4.M       <-stf(window(cod4.M,end=2010),nyr=21)
m(cod4.M)[,ac(2011:2031)]<-M_IncrOM
cod4.M       <-fwd(cod4.M,trgt,sr.params=srParHalf,sr.model="ricker",sr.residuals=sr.deviates)
stock.n(cod4.M)[,ac(2012:2031)]<-NA                               #### Ugly fudge to solve stock.n being 0 instead of NA.

biol         <-as(cod4.M,"FLBiol")
fleet        <-as(cod4.M,"FLFleet")
price(fleet,1,1)[]<-c(1.05,1.05,1.40,2.01,2.759,2.759,3.39)
save(fleet,biol,file=paste(my.dir,"/data/OMM_lowR.RData",sep=""))

## F2011 = F2010*.75
cod4.Catch   <-fwd(stf(window(cod4.Catch,end=2010),nyr=21),trgt,sr.params=srParHalf,sr.model="ricker",sr.residuals=sr.deviates)
stock.n(cod4.Catch)[,ac(2012:2031)]<-NA                           #### Ugly fudge to solve stock.n being 0 instead of NA.

biol              <-as(cod4.Catch,"FLBiol")
fleet             <-as(cod4.Catch,"FLFleet")
price(fleet,1,1)[]<-c(1.05,1.05,1.40,2.01,2.759,2.759,3.39)
save(fleet,biol,file=paste(my.dir,"/data/OMCatch_lowR.RData",sep=""))
save(cod4.Catch,file=paste(my.dir,"/data/Stocks_lowR.RData",sep=""))

## Deviates ####################################################################
## options
## SR
nits     <-1000
StartYr  <-2011 ## i.e. current year
EndYr    <-2030
srCV      <-0.6
srDeviates<-FLQuant(exp(rnorm(nits*length(StartYr:(EndYr+1)), mean=0, sd=srCV)),
                     dimnames=list(age="rec",year=StartYr:(EndYr+1),iter=1:nits))/exp((srCV^2)/2)
## CPUE
idxCV      <-0.3
idxDms     <-list(age=1:range(cod4.2010)["plusgroup"],year=1978:2031,iter=1:nits)
### The above was changed by MPH on 17_2_09 from:
# idxDms     <-list(age=1:plusgroup,year=1963:2028,iter=1:nits)   

idxDeviates<-FLQuant(exp(idxCV*rnorm(prod(unlist(lapply(idxDms,length))))),dimnames=idxDms)
#idxDeviates.Q1<-sweep(iter(Q1.res,1:nits),c(1,6),sqrt(iter(Q1.var,1:nits)),"*")
#idxDeviates.Q3<-sweep(iter(Q3.res,1:nits),c(1,6),sqrt(iter(Q3.var,1:nits)),"*")

save(idxDeviates,srDeviates,file=paste(my.dir,"/data/Deviates.RData",sep=""))
#rm(list=ls())


plot(rec(cod4.sr)~ssb(cod4.sr),pch=19,xlab="SSB",ylab="Recruits at age 1",xlim=c(0,max(ssb(cod4.sr))),ylim=c(0,max(rec(cod4.sr))))
lines(predict(cod4.sr,ssb=FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))~FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))
savePlot(file=paste(my.dir,"/data/SR1.jpeg",sep=""))

plot(rec(cod4.sr)~ssb(cod4.sr),pch=19,xlab="SSB",ylab="Recruits at age 1",xlim=c(0,max(ssb(cod4.sr))),ylim=c(0,max(rec(cod4.sr))))
lines(0.5*predict(cod4.sr,ssb=FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))~FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))
savePlot(file=paste(my.dir,"/data/SR2.jpeg",sep=""))

plot(rec(cod4.sr)~ssb(cod4.sr),pch=19,xlab="SSB",ylab="Recruits at age 1",xlim=c(0,max(ssb(cod4.sr))),ylim=c(0,max(rec(cod4.sr))))
lines(predict(cod4.sr,ssb=FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))~FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))
lines(0.5*predict(cod4.sr,ssb=FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))~FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))
savePlot(file=paste(my.dir,"/data/SR1&2_WOScod.jpeg",sep=""))
