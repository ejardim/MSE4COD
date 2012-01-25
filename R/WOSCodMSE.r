library(FLCore)
library(FLash)
library(FLBRP)
library(FLXSA)  # Also loads the FLAssess package
#library(DBI)
#library(RSQLite)

packageDescription("FLCore",  fields="Built")
packageDescription("FLAssess",fields="Built")
packageDescription("FLash",   fields="Built")
packageDescription("FLBRP",   fields="Built")
packageDescription("FLXSA",   fields="Built")

memory.limit(size=4000)

# Activate the history function for the plots
windows(180,50,record=T)

my.dir<-"C:/FLR/mse4ej"

## DB
#SQLite(max.con = 16, fetch.default.rec = 500, force.reload = FALSE, shared.cache=FALSE)

## connect DBs
#dbOM  <-paste(my.dir,"/db/hcrEU.dbf", sep="")
#conOM <-dbConnect(dbDriver("SQLite"), dbname=dbOM)

## files & dirs
source(paste(my.dir,"/R/mseFuncs.R",        sep=""))
source(paste(my.dir,"/R/HCR3.R",            sep=""))
load(  paste(my.dir,"/data/condSR.RData",   sep=""))
load(  paste(my.dir,"/data/condM.RData",    sep=""))
load(  paste(my.dir,"/data/condCatch.RData",sep=""))   #no need to adjust misrep (as for IS cod) because HC landings in 2007 close enough to TAC
load(  paste(my.dir,"/data/Deviates.RData", sep=""))

runs <- rbind(expand.grid(OEM=c("wg","m","catch"),OM=c("m","catch"),SR=c(1,0.5),HCR=c("EU"),fcut=c(0.75),minTAC=0.8,maxTAC=1.2),
              expand.grid(OEM=c("catch"),OM=c("catch"),SR=c(1,0.5),HCR=c("EU"),fcut=c(0.75),minTAC=0.001,maxTAC=1000))

nits     <-1:250
StartYr  <-2011 ## i.e. current year
EndYr    <-2030
Bpa      <-22000              ### changed MPH
Blim     <-14000              ### changed MPH
FTarget  <-0.4
plusgroup<-7
biolrun<-list()
fleetrun<-list()

for (i in c(1:2))        #c(20:24,37)
    {
    OM     <-runs[i,"OM", drop=T]
    OEM    <-runs[i,"OEM",drop=T]
    HCR    <-runs[i,"HCR",drop=T]
    SR     <-runs[i,"SR", drop=T]
    minTAC <-runs[i,"minTAC", drop=T]
    maxTAC <-runs[i,"maxTAC", drop=T]
    fcut <-runs[i,"fcut", drop=T]

    print(paste(i, OEM, OM, HCR, SR, fcut, minTAC, maxTAC))
    
    if (SR==1) {
      if (OM=="m")
        load(paste(my.dir,"/data/OMM_fullR.RData",sep="")) else
        load(paste(my.dir,"/data/OMCatch_fullR.RData",sep=""))
    } else {
      if (OM=="m")
        load(paste(my.dir,"/data/OMM_lowR.RData",sep="")) else
        load(paste(my.dir,"/data/OMCatch_lowR.RData",sep=""))
    }

    biol <-iter(biol, nits)

    #fleet<-iter(fleet,nits)
    fleet@metiers[[1]]@catches[[1]]         <-iter(fleet@metiers[[1]]@catches[[1]],nits)
    fleet@metiers[[1]]@effshare[,,,,,nits]<-fleet@metiers[[1]]@effshare[,,,,,nits]
    fleet@effort                            <-fleet@effort[,,,,,nits]

    #### Observation Error Model
    ## OEM: M
    if (OEM=="m")
       MP.M <-M_Incr else
       MP.M <-M_WG

    ## OEM: Catch
    misrep<-FLQuant(1, dimnames=dimnames(catchMisrepMP))
    if (OEM=="catch")
       misrep<-misrep*catchMisrepMP
    if (OM=="catch")
       misrep<-misrep/catchMisrepOM
    misrep<-iter(misrep,nits)
    
    OM.srPar      <-srPar
    OM.srPar[,"a"]<-OM.srPar[,"a"]*SR
    ############################################################################

    #### Setup MP ##############################################################
    ## iYr   = last data year
    ## iYr+1 = now
    ## iYr+2 = TAC year

    #### XSA Options
    xsaCtrl <-FLXSA.control(qage=4,shk.n=F,tsrange=100,tspower=0,fse=2.0,shk.yrs=3,shk.ages=3)  #to weaken shrinkage, set fse=2.0

    #### CPUE
    MPIdx      <-OEMCPUE.RV(biol,iter(idxDeviates,nits),start=dims(biol)$minyear,end=StartYr-2,startf=0.0,endf=0.1)
    MPIdx@range<-MPIdx@range[-c(6:7)] #need this to make FLXSA work
    MPIdx      <-window(MPIdx,end=EndYr)

    #### Stock
    MPStk<-OEMStock(biol,fleet,start=dims(biol)$minyear,end=StartYr-2,M=MP.M,misrep=misrep) #,plusgroup=plusgroup)

    #### HCR summary
    smryHCR<-array(c(0.75,0.4,NA,NA,NA,HCR),dim=c(6,EndYr-StartYr+3,length(nits)),dimnames=list(val=c("lambda","f","ssb","TAC","landings","rule"),year=(StartYr-1):(EndYr+1),iter=nits))
    smryHCR["TAC",ac(2010),]<- 240               ##### EC23/2010
    smryHCR["TAC",ac(2011),]<- 182               ##### EC57/2011
    if (OM=="catch") {
        smryHCR["landings",ac(2010),]<-smryHCR["TAC",ac(2010),]*catchMisrepOM[,ac(2010),,,,nits]
        smryHCR["landings",ac(2011),]<-smryHCR["TAC",ac(2011),]*catchMisrepOM[,ac(2011),,,,nits]
      } else smryHCR["landings",ac(2010:2011),]<-smryHCR["TAC",ac(2010:2011),]

    for (iYr in ((StartYr-1):(EndYr-2)))
       {
         cat("=============",iYr," ============\n")
         #### Misreporting
         if (OM=="catch"){
            misrep[,ac(iYr)]<-pmin(1,c(smryHCR["TAC",ac(iYr),]/smryHCR["landings",ac(iYr),]))
            if (OEM=="catch")
               misrep[,ac(iYr)]<-c(catchMisrepMP[,ac(iYr),,,,nits])*misrep[,ac(iYr)]
            }
  
         ## Assessment
         MPIdx@index[,ac(iYr)]<-OEMCPUE.RV(biol,iter(idxDeviates,nits),start=iYr,startf=0.0,endf=0.01)@index
         MPStk                <-window(MPStk,end=iYr)
         MPStk[,ac(iYr)]      <-OEMStock(biol,fleet,start=iYr,M=MP.M,misrep=misrep) #,plusgroup=plusgroup)
         MPStk                <-MPStk+FLXSA(MPStk,MPIdx,xsaCtrl,diag.flag=FALSE)
  
         ## HCRs #################################################################
         if (HCR=="EU") s.<-hcrEU(iYr,smryHCR,MPStk,Bpa,Blim,0.4,catchMisrepMP,OEM,fcut,minTAC,maxTAC)
  
         MPStk  <-s.$m
         smryHCR<-s.$s

         if (OM=="catch")
            smryHCR["landings",ac(iYr+2),]<-smryHCR["TAC",ac(iYr+2),]*(catchMisrepOM[,ac(iYr+2),,,,nits]) else
            smryHCR["landings",ac(iYr+2),]<-smryHCR["TAC",ac(iYr+2),]
  
         ## Go Fish
         if (HCR %in% c("EU","Norway")){
           trgtVal  <-array(c(smryHCR["landings",ac(iYr+2),,drop=T]),c(1,length(nits)))
           OMtarget <-fwdTarget( year=iYr+2, value=mean(trgtVal), quantity="landings", fleet=1, metier=1, spp=1)
           OMcontrol<-fwdControl(year=iYr+2, fleet=1,metier=1,min=.01,max=1.0)
           } else
         if (HCR == "f") {
           trgtVal  <-array(c(smryHCR["f",ac(iYr+2),,drop=T]),c(1,length(nits)))
           OMtarget <-fwdTarget(year=iYr+2, value=NA, quantity="f", fleet=1, metier=1, spp=1)
           OMcontrol<-fwdControl(year=iYr+2,fleet=1,metier=1,value=NA,min=.01,max=.75)
           }
  
         res<-fwd(biol,fleet,OMtarget,OMcontrol,target.value=trgtVal,sr.model="ricker",sr.params=OM.srPar,sr.residuals=srDeviates)
  
         landings.n(fleet,1,1)[,ac(iYr+2)]<-res$landings.n[,ac(iYr+2)]
         discards.n(fleet,1,1)[,ac(iYr+2)]<-res$discards.n[,ac(iYr+2)]
         landings(  fleet,1,1)[,ac(iYr+2)]<-apply(res$landings.n[,ac(iYr+2)]*landings.wt(fleet,1,1)[,ac(iYr+2)],2:6,sum)
         discards(  fleet,1,1)[,ac(iYr+2)]<-apply(res$discards.n[,ac(iYr+2)]*discards.wt(fleet,1,1)[,ac(iYr+2)],2:6,sum)
  
         effort(fleet)[,ac(iYr+2)]<-res$effort[,ac(iYr+2)]
         n(biol)[      ,ac(iYr+2)]<-res$n[     ,ac(iYr+2)]
  
         if (length(nits)==1) {
            SmryPlot1(fleet,biol,res$f,MPStk,yrs=1990:(iYr),its=1)
            savePlot(paste(my.dir,"/figs/scenario",formatC(i,width=2,flag="0"),"iter",formatC(nits,width=3,flag="00"),"-",iYr,".jpeg",sep=""),type="jpeg")
         } else {
            SmryPlots(fleet,biol,res$f,MPStk,yrs=1990:(iYr))
            savePlot(paste(my.dir,"/figs/scenario",formatC(i,width=2,flag="0"),"-",iYr,".jpeg",sep=""),type="jpeg")
         }
       }

#    if (length(nits)>1)
#       savePlot(paste(my.dir,"figs/Run",i,".jpeg",sep=""),type="jpeg")

    print(apply(smryHCR[,ac(2010:min((iYr+2),EndYr)),],1:2,mean))

    ## Save summaries to DB
    dfOM<-cbind(OEM=OEM,OM=OM,HCR=HCR,SR=SR,fcut=fcut,model.frame(smryStats(fleet,biol,start=StartYr-5,end=EndYr))[,-c(1,3:5)])
    #dbWriteTable(conOM, "OM",  dfOM, append=T)

    dfHCR<-cbind(OEM=OEM,OM=OM,HCR=HCR,SR=SR,fcut=fcut,expand.grid(dimnames(smryHCR)[-1]),
                 MPlambda  =c(smryHCR["lambda",,]),
                 MPTAC     =c(smryHCR["TAC",,]),
                 MPlandings=c(smryHCR["landings",,]),
                 MPssb     =c(smryHCR["ssb",,]),
                 MPf       =c(smryHCR["f",,]),
                 MPrule    =c(smryHCR["rule",,]))

    write.csv(dfHCR,file=paste(my.dir,"/tables/scenario",formatC(i,width=2,flag="0"),".csv",sep=""))
    #dbWriteTable(conOM, "HCR", dfHCR, append=T)
    biolrun[[i]] <- biol
    fleetrun[[i]] <- fleet
    }

## close and clean up databases
#dbListTables(conOM)
#dbDisconnect(conOM)
#file.info(dbOM)

save(biolrun,fleetrun,file=paste(my.dir,"/data/biolfleetWOScod_1-6.RData",sep=""))

