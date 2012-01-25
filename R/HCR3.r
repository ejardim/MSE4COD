hcrEU<-function(iYr,smryHCR,MPStk,Bpa,Blim,FTarget,catchMisrepMP,OEM,fcut,minTAC=0.8,maxTAC=1.2,lambda=TRUE)
   {
   nits<-dims(MPStk)$iter

   #Short-term forecast
   MPStk<-stf(MPStk,nyrs=2)
   mnRec<-FLPar(c(exp(apply(log(stock.n(MPStk)[1,ac(iYr-1:10)]),6,mean))))

   ## SSB at start of year (y+1)
   FCurrent<-array(c(smryHCR["lambda",ac(iYr),]*fbar(MPStk)[,ac(iYr)]),c(1,nits))
   trgt    <-fwdTarget(year=c(iYr+1),quantity="f",value=mean(FCurrent))
   MPStk   <-fwd(MPStk,trgt,value=FCurrent,sr.model="geomean",sr.params=mnRec)
   smryHCR["ssb",ac(iYr+1),]<-c(ssb(MPStk)[,ac(iYr+1)])

   ## 1) Target F in TAC year (y+2)
   smryHCR["f",  ac(iYr+2),                              ]<-c(pmax(0.90*FCurrent,FTarget,na.rm=T))
   smryHCR["f",  ac(iYr+2),c(ssb(MPStk)[,ac(iYr+1)])<Bpa ]<-c(pmax(0.85*FCurrent,FTarget))[c(ssb(MPStk)[,ac(iYr+1)])<Bpa]
   smryHCR["f",  ac(iYr+2),c(ssb(MPStk)[,ac(iYr+1)])<Blim]<-c(     fcut*FCurrent         )[c(ssb(MPStk)[,ac(iYr+1)])<Blim]

   ## 2) TAC next year (y+2) with TAC constraints
   value<-array(NA,c(2,nits))
   min  <-array(NA,c(2,nits))
   max  <-array(NA,c(2,nits))
   value[1,]<-c(smryHCR["f",  ac(iYr+2),])
   min[2,]  <-c(smryHCR["TAC",ac(iYr+1),1:nits])*minTAC
   max[2,]  <-c(smryHCR["TAC",ac(iYr+1),1:nits])*maxTAC

   ## turns off TAC constraint if SSB<Blim
#   tacFlag<-c(ssb(MPStk)[,ac(iYr+1)])<=Blim
#   if (any(tacFlag)){
#      max[2,tacFlag]<-max[2,tacFlag]*1000
#     min[2,tacFlag]<-min[2,tacFlag]*0.001
#      }

   if (OEM=="catch") {
     trgt <- fwdTarget(year=c(iYr+2),quantity=c("f"),value=c(mean(value[1,])))
     MPStk <- fwd(MPStk,trgt,value=array(c(value[1,]),c(1,nits)),sr.model="geomean",sr.params=mnRec)
     value[1,] <- c(computeLandings(MPStk)[,ac(iYr+2)])/c(catchMisrepMP[,ac(iYr+2),,,,1:nits])
     trgt <- fwdTarget(year=rep(iYr+2,2),quantity=c("landings","landings"),value=c(mean(value[1,]),NA),min=c(NA,c(mean(min[2,]))),max=c(NA,c(mean(max[2,]))))
   } else trgt <- fwdTarget(year=rep(iYr+2,2),quantity=c("f","landings"),value=c(mean(value[1,]),NA),min=c(NA,c(mean(min[2,]))),max=c(NA,c(mean(max[2,]))))
  MPStk <- fwd(MPStk,trgt,value=value,min=min,max=max,sr.model="geomean",sr.params=mnRec)

   ## set TAC
   smryHCR["TAC",ac(iYr+2),]<-c(computeLandings(MPStk)[,ac(iYr+2)])
   MPStklambda<-MPStk
   if (OEM=="catch") {
      value[1,] <- c(computeLandings(MPStk)[,ac(iYr+2)])*c(catchMisrepMP[,ac(iYr+2),,,,1:nits])
      trgt <- fwdTarget(year=c(iYr+2),quantity=c("landings"),value=c(mean(value[1,])))
      MPStklambda <- fwd(MPStk,trgt,value=array(c(value[1,]),c(1,nits)),sr.model="geomean",sr.params=mnRec)
   }
   smryHCR["f",  ac(iYr+2),]<-c(           fbar(MPStklambda)[,ac(iYr+2)])

   ## 3) Recalculate lambda
   if (lambda)
      {smryHCR["lambda",ac(iYr+1),smryHCR["f",ac(iYr+2),]==0]<-0.
       smryHCR["lambda",ac(iYr+1),FCurrent[1,]==0&smryHCR["f",ac(iYr+2),]>0]<-10.
       smryHCR["lambda",ac(iYr+1),FCurrent[1,]>0&smryHCR["f",ac(iYr+2),]>0]<-
         smryHCR["f",ac(iYr+2),FCurrent[1,]>0&smryHCR["f",ac(iYr+2),]>0]/FCurrent[1,FCurrent[1,]>0&smryHCR["f",ac(iYr+2),]>0]} else
      smryHCR["lambda",ac(iYr+1),]<-1.0

   #print(sprintf("%5.0f: lambda =%5.3f: F(y+2)=%5.3f: F(y+2)=%5.3f", iYr,c(mean(smryHCR["lambda",ac(iYr+1),])),c(mean(smryHCR["f",ac(iYr+2),])),c(mean(FCurrent[1,]))))

   return(list(s=smryHCR,m=MPStk))
   }

