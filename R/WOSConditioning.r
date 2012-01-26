###############################################################################
# Original LK & JdO revised by EJ & IM (20120126) 
# Conditioning WoS MSE
###############################################################################

library(FLCore)
library(FLash)
library(FLAssess)

#==============================================================================
# Read data and create objects
#==============================================================================

my.dir<-"../"

# read in WG data
cod4.2010 <- readFLStock(paste(my.dir,"input/cod6a.IDX",sep=""))
catch.n(cod4.2010) <- readVPAFile(paste(my.dir,"input/cod6aCN.DAT",sep=""), quiet=FALSE)
catch.wt(cod4.2010) <- readVPAFile(paste(my.dir,"input/cod6aCW.DAT",sep=""), quiet=FALSE)
discards.n(cod4.2010) <- readVPAFile(paste(my.dir,"input/cod6aDN.DAT",sep=""), quiet=FALSE)
discards.wt(cod4.2010) <- readVPAFile(paste(my.dir,"input/cod6aDW.DAT",sep=""), quiet=FALSE)

range(cod4.2010,c("minfbar","maxfbar"))<-c(2,5)
cod4.2010 <- setPlusGroup(cod4.2010 ,7, na.rm=TRUE)

# Get bootstrap estimates
f.hat <-read.table(paste(my.dir,"Badapt used for MSE/F.CSV",sep=""),header=FALSE,sep=",") [,-c(1,35)]
n.hat <-read.table(paste(my.dir,"Badapt used for MSE/NUMBERS.CSV",sep=""),header=FALSE,sep=",") [,-c(1,35,36)]

#--------------------------------------------------------------------
# using first iter for original named '0'
# dangerous, it's bypassing validity
# at least name iter '0' 'as '1'
# better option to use 2 objects, one original, one bootstrap
#--------------------------------------------------------------------

dmns <- expand.grid(age=1:7,iter=0:1000,year=1978:2010)
harvest(cod4.2010) <- as.FLQuant(cbind(data=as.numeric(as.matrix(f.hat)),dmns),units="f")
stock.n(cod4.2010) <- as.FLQuant(cbind(data=as.numeric(as.matrix(n.hat)),dmns))
rm(f.hat,n.hat)

# Compare bootstrap and initial runs
jpeg(paste(my.dir,"/data/KeyRun.jpeg",sep=""))
plot(FLStocks(window(iter(cod4.2010,2:1001), end=2010), window(iter(cod4.2010,1), end=2010)), ylab="Year")
dev.off()

#==============================================================================
# S/R
#==============================================================================

# fit
cod4.sr <- as.FLSR(window(iter(cod4.2010,1), end=2010), model='ricker')
cod4.sr <- fmle(cod4.sr)

jpeg(file=paste(my.dir,"/data/SR.jpeg",sep=""))
plot(cod4.sr)
dev.off()

# remove 1st run
cod4.2010 <- iter(cod4.2010,2:1001)

# generate SR deviates
sr.deviates <- FLQuant(rlnorm(20000, sd=0.66), dimnames=list(age=1,year=2011:2030,iter=1:1000))

## project (EJ: why ?? is it just for the plot ??)
#cod4 <- stf(cod4.2010, nyears=20)
#trgt <- fwdControl(data.frame(year=2011:2030, quantity="f", val=0.4))
#cod4 <- fwd(cod4, trgt, sr=cod4.sr, sr.residuals=sr.deviates)

## repeat for halved Recruitment
#cod4.2.sr <- cod4.sr 
#params(cod4.2.sr)["a"] <- params(cod4.2.sr)["a"]/2
#cod4.2 <- fwd(cod4, trgt, sr=cod4.2.sr, sr.residuals=sr.deviates)

#jpeg(file=paste(my.dir,"/data/ProjSRs.jpeg",sep=""))
#plot(FLStocks(cod4,cod4.2))
#dev.off()

#==============================================================================
# MisReporting ratios
#==============================================================================

catchMisrep <- stock.n(cod4.2010)*harvest(cod4.2010)/(harvest(cod4.2010)+m(cod4.2010))*(1-exp(-(harvest(cod4.2010)+m(cod4.2010))))

jpeg(file=paste(my.dir,"/data/MisRepRatio.jpeg",sep=""))
histogram( ~ data | as.character(year), data = as.data.frame(window(catchMisrep/catch.n(cod4.2010),start=1994)), xlab="Mis-reporting Ratio")
dev.off()
### Note that early years ratios are close to 1 so appears as all being 1.

#==============================================================================
# Mis-specification of catch
#==============================================================================

cod4.Catch <- cod4.2010
catch.n(cod4.Catch) <- stock.n(cod4.Catch)*harvest(cod4.Catch)/(harvest(cod4.Catch)+m(cod4.Catch))*(1-exp(-(harvest(cod4.Catch)+m(cod4.Catch))))
landings.n(cod4.Catch) <- catch.n(cod4.Catch)*landings.n(cod4.Catch)/(discards.n(cod4.Catch)+landings.n(cod4.Catch))
discards.n(cod4.Catch) <- catch.n(cod4.Catch)-landings.n(cod4.Catch)
stock(cod4.Catch) <- computeStock(cod4.Catch)
catch(cod4.Catch) <- computeCatch(cod4.Catch, "all")
landings(cod4.Catch) <- computeLandings(cod4.Catch)
discards(cod4.Catch) <- computeDiscards(cod4.Catch)

#==============================================================================
# Mis-specification of M
#==============================================================================
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

#==============================================================================
# projection for mean recruitment
#==============================================================================
#--------------------------------------------------------------------
# Why is this projection being done ?
#--------------------------------------------------------------------
## average Recruitment
#cod4.Catch <- stf(cod4.Catch, nyears=20)
#cod4.M <- stf(cod4.M, nyears=20)
#trgt <- fwdControl(data.frame(year=2011:2030, quantity="f", val=0.4))
#cod4.Catch.1 <- fwd(cod4.Catch, trgt, sr=cod4.sr, sr.residuals=sr.deviates)
#cod4.M.1 <- fwd(cod4.M, trgt, sr=cod4.sr, sr.residuals=sr.deviates)

#jpeg(file=paste(my.dir,"/data/Proj1.jpeg",sep=""))
#plot(window(FLStocks(cod4.Catch.1,cod4.M.1),start=2000),xlab="Year")
#dev.off()

#==============================================================================
# OM Misreporting/M scenarios
#==============================================================================
#--------------------------------------------------------------------
# This code looks odd ... will take a look later
# It seems to bootstrap the misreporting ratios but not considering.
# It's considering ratios to be independent ...
#--------------------------------------------------------------------
yrsMC <- sample(1994:2010, 21*1000, replace=T) + rep((1:1000)*10000, each=21)
yrsReformat <- rep(1994:2010,1000) + rep(1:1000, each=17)*10000

## Catch misreporting
catchMisrepOM <- FLQuant(1, dimnames=list(age="all",year=1978:2031,iter=1:1000))
tmp <- window(catchMisrep/catch.n(cod4.2010), start=1994)
tmp <- FLQuant(c(tmp), dimnames=list(age=1:7,year=yrsReformat))
catchMisrepOM[,ac(2011:2031)] <- FLQuant(c(tmp[,ac(yrsMC)]),dimnames=list(age=1:7,year=2011:2031,iter=1:1000))[1]
catchMisrepOM[,ac(1978:2010)] <- apply(catchMisrep/catch.n(cod4.2010), c(2,6), mean)
catchMisrepMP <- propagate(apply(catchMisrepOM,2,median),iter=1000)

save(catchMisrepOM, catchMisrepMP, file=paste(my.dir,"/data/condCatch.RData",sep=""))

## M mis-specification
tmp <- window(m(cod4.M), start=1994)
tmp <- FLQuant(c(tmp), dimnames=list(age=1:7,year=yrsReformat))
M_IncrOM <- FLQuant(c(tmp[,ac(yrsMC)]), dimnames=list(age=1:7,year=2011:2031,iter=1:1000))
M_WG <- FLQuant(c(apply(window(m(cod4.Catch), start=1994), 1, mean)), dimnames=list(age=1:7,year=1994:2031))
M_Incr <- M_WG
M_Incr[,ac(1994:2010)] <- apply(window(m(cod4.M), start=1994), 1:2, median)
M_Incr[,ac(2011:2031)] <- apply(M_IncrOM, 1:2, median)

save(M_WG,M_Incr,catchMisrepMP,file=paste(my.dir,"/data/condM.RData",sep=""))

## intermediate year, full recruitment
trgt <- fwdControl(data.frame(year=2011, quantity="f", val=0.75, rel=2010)) #compare EC23/2010 to EC57/2011

## F2011 = F2010*.75
cod4.M <- stf(window(cod4.M, end=2010), nyears=21)
m(cod4.M)[,ac(2011:2031)] <- M_IncrOM
cod4.M <- fwd(cod4.M, trgt, sr=cod4.sr, sr.residuals=sr.deviates)
stock.n(cod4.M)[,ac(2012:2031)] <- NA ## Ugly fudge to solve stock.n being 0 instead of NA.
biol <- as(cod4.M,"FLBiol")
fleet <- as(cod4.M,"FLFleet")
price(fleet,1,1)[] <- c(1.05,1.05,1.40,2.01,2.759,2.759,3.39)

save(fleet, biol, file=paste(my.dir,"/data/OMM_fullR.RData",sep=""))

## F2011 = F2010*.75
cod4.Catch <- fwd(stf(window(cod4.Catch,end=2010),nyears=21), trgt, sr=cod4.sr, sr.residuals=sr.deviates)
stock.n(cod4.Catch)[,ac(2012:2031)] <- NA ## Ugly fudge to solve stock.n being 0 instead of NA.
biol <- as(cod4.Catch,"FLBiol")
fleet <- as(cod4.Catch,"FLFleet")
price(fleet,1,1)[]<-c(1.05,1.05,1.40,2.01,2.759,2.759,3.39)

save(fleet,biol,file=paste(my.dir,"/data/OMCatch_fullR.RData",sep=""))

## intermediate year, low recruitment
trgt <- fwdControl(data.frame(year=2011, quantity="f", value=0.75, rel=2010))

## F2011 = F2010*.75
cod4.M <- stf(window(cod4.M, end=2010), nyears=21)
m(cod4.M)[,ac(2011:2031)] <- M_IncrOM
cod4.M <- fwd(cod4.M, trgt, sr=cod4.sr, sr.residuals=sr.deviates)
stock.n(cod4.M)[,ac(2012:2031)] <- NA ## Ugly fudge to solve stock.n being 0 instead of NA.
biol <- as(cod4.M,"FLBiol")
fleet <- as(cod4.M,"FLFleet")
price(fleet,1,1)[]<-c(1.05,1.05,1.40,2.01,2.759,2.759,3.39)

save(fleet,biol,file=paste(my.dir,"/data/OMM_lowR.RData",sep=""))

## F2011 = F2010*.75
cod4.Catch <- fwd(stf(window(cod4.Catch,end=2010), nyears=21), trgt, sr=cod4.sr, sr.residuals=sr.deviates)
stock.n(cod4.Catch)[,ac(2012:2031)] <- NA ## Ugly fudge to solve stock.n being 0 instead of NA.
biol <- as(cod4.Catch,"FLBiol")
fleet <- as(cod4.Catch,"FLFleet")
price(fleet,1,1)[] <- c(1.05,1.05,1.40,2.01,2.759,2.759,3.39)

save(fleet,biol,file=paste(my.dir,"/data/OMCatch_lowR.RData",sep=""))
save(cod4.Catch,file=paste(my.dir,"/data/Stocks_lowR.RData",sep="")) # EJ: do we need this object ?

#==============================================================================
# Abundance deviates ? 
#==============================================================================
#----------------------------------------------------------
# Why are srDeviates being generated again ?
# some confusion between CV and SD better to sort it out
#----------------------------------------------------------
nits <- 1000
StartYr <- 2011 ## i.e. current year
EndYr <- 2030
srCV <- 0.6
srDeviates <- FLQuant(rlnorm(nits*length(StartYr:(EndYr+1)), mean=0, sd=srCV), dimnames=list(age="rec",year=StartYr:(EndYr+1), iter=1:nits))
## CPUE
idxCV <- 0.3
idxDms <- list(age=1:range(cod4.2010)["plusgroup"], year=1978:2031, iter=1:nits)
idxDeviates <- FLQuant(exp(idxCV*rnorm(prod(unlist(lapply(idxDms,length))))), dimnames=idxDms)

save(idxDeviates,srDeviates,file=paste(my.dir,"/data/Deviates.RData",sep=""))

#==============================================================================
# PLOTS
#==============================================================================

jpeg(file=paste(my.dir,"/data/SR1.jpeg",sep=""))
plot(rec(cod4.sr)~ssb(cod4.sr),pch=19,xlab="SSB",ylab="Recruits at age 1",xlim=c(0,max(ssb(cod4.sr))),ylim=c(0,max(rec(cod4.sr))))
lines(predict(cod4.sr,ssb=FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))~FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))
dev.off()

jpeg(file=paste(my.dir,"/data/SR2.jpeg",sep=""))
plot(rec(cod4.sr)~ssb(cod4.sr),pch=19,xlab="SSB",ylab="Recruits at age 1",xlim=c(0,max(ssb(cod4.sr))),ylim=c(0,max(rec(cod4.sr))))
lines(0.5*predict(cod4.sr,ssb=FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))~FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))
dev.off()

jpeg(file=paste(my.dir,"/data/SR1&2_WOScod.jpeg",sep=""))
plot(rec(cod4.sr)~ssb(cod4.sr),pch=19,xlab="SSB",ylab="Recruits at age 1",xlim=c(0,max(ssb(cod4.sr))),ylim=c(0,max(rec(cod4.sr))))
lines(predict(cod4.sr,ssb=FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))~FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))
lines(0.5*predict(cod4.sr,ssb=FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))~FLQuant(seq(0,max(ssb(cod4.sr)),length.out=100)))
dev.off()

