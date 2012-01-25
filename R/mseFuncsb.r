## Calculates additional stats in function smryStats
## function to apply a linearly increasing trend to an FLQuant
biasLinear<-function(x,obj)
   {
   if (x>0)
      res  <-1-(sort(cumsum(rep(x, dims(obj)$year)),d=T)-x)
   else
      res  <-sort(1-(sort(cumsum(rep(x, dims(obj)$year)),d=T)-x),d=T)

   return(obj*FLQuant(rep(res,each=dims(obj)$age),dimnames=dimnames(obj)))
   }

## Calculate misreporting F
computeCatch.n<-function(x)
   {
   f<-harvest(x)
   m<-m(x)
   z<-m+f
   n<-stock.n(x)

   catchN<-n*f/(z)*(1-exp(-z))

   return(catchN)
   }

## Observation Error Model
OEMStock<-function(biol,fleet,start,M="missing",misrep="missing",end="missing",plusgroup="missing")
    {
    if (missing(end)) end<-start
    yrs<-as.character(end:start)

    units(discards(fleet@metiers[[1]]@catches[[1]])) <- units(landings(fleet@metiers[[1]]@catches[[1]]))
    
    stk <-as(window(fleet@metiers[[1]]@catches[[1]], start=start, end=end),"FLStock")

    catch.n(stk)[,yrs]   <-catch.n(fleet,1,1)[   ,yrs]
    discards.n(stk)[,yrs]<-discards.n(fleet,1,1)[,yrs]
    landings.n(stk)[,yrs]<-landings.n(fleet,1,1)[,yrs]
    if (!missing(misrep) & any(yrs %in% dimnames(misrep)$year))
       {
       yrs.                  <-yrs[yrs %in% dimnames(misrep)$year]
       catch.n(stk)[   ,yrs.]<-sweep(catch.n(   stk)[,yrs.],c(2,6),misrep[,yrs.],"*")
       discards.n(stk)[,yrs.]<-sweep(discards.n(stk)[,yrs.],c(2,6),misrep[,yrs.],"*")
       landings.n(stk)[,yrs.]<-sweep(landings.n(stk)[,yrs.],c(2,6),misrep[,yrs.],"*")
       }
       
    ## biological parameters
    stock.wt(    stk)[,yrs] <-wt(  biol)[dimnames(stk@m)$age,yrs]
    mat(         stk)[,yrs] <-fec( biol)[dimnames(stk@m)$age,yrs]
    harvest.spwn(stk)[,yrs] <-spwn(biol)[dimnames(stk@m)$age,yrs]
    m.spwn(      stk)[,yrs] <-spwn(biol)[dimnames(stk@m)$age,yrs]
    m(           stk)[,yrs]<-m(biol)[dimnames(stk@m)$age,yrs]
    if (!missing(M))
      {
      yrs.<-yrs[yrs %in% dimnames(M)$year]
      m(stk)[,yrs.]<-M[,yrs.]
      }

    if (dims(stk)$iter != dims(discards.wt(stk))$iter) discards.wt(stk)<-propagate(discards.wt(stk), iter=dims(catch.n(stk))$iter)
    if (dims(stk)$iter != dims(landings.wt(stk))$iter) landings.wt(stk)<-propagate(landings.wt(stk), iter=dims(catch.n(stk))$iter)

    if (!missing(plusgroup))
       stk<-setPlusGroup(stk,plusgroup)

    return(stk)
    }

## Research survey indices
OEMCPUE.RV<-function(biol,deviates,start,end="missing",plusgroup="missing",startf="missing",endf="missing")
     {
     ## unbiased population estimates
     
     if (missing(end)) end<-start
     yrs<-start:end

     biol <-window(biol,start=start,end=end)
     idx  <-as(biol,"FLIndex")

     if (!missing(startf)) idx@range["startf"]<-startf
     if (!missing(endf))   idx@range["endf"]  <-endf

     if (!missing(plusgroup))
        index<-setPlusGroup(biol,plusgroup)@n
     else
        index<-biol@n
        
     idx@index<-index*deviates[dimnames(idx@index)$age,ac(yrs)]
 
     return(idx)
     }

## Summary statistics
smryStats<-function(fleet,biol, start=-30, end=0, brp="missing")
   {
   fleet<-window(fleet, start=start, end=end)
   biol <-window(biol , start=start, end=end)

   cost    <-function(fleet) {fleet@fcost+fleet@metiers[[1]]@vcost*fleet@effort}

   revenues=revenue(fleet@metiers[[1]]@catches[[1]])
   #costs   =cost(fleet)
   #profits =revenues-costs
   yield   =computeLandings(fleet,1)[[1]][[1]]
   catch   =apply(catch.n(fleet,1,1)*catch.wt(fleet,1,1),c(2,6),sum)
   discards=computeDiscards(fleet,1)[[1]][[1]]
   rec     =n(biol)[1,]
   ssb     =ssb(biol)
   biomass =apply(biol@n*biol@wt,c(2,6),sum)
   mnsz    =apply(biol@n*biol@wt,c(2,6),sum)/apply(biol@n,c(2,6),sum)
   pmat    =ssb/biomass
   harvest =catch/biomass
   fLbar<-apply(sweep(landings.sel(fleet,1,1),c(2,6),fleet@effort,"*")[ac(2:4)],c(2,6),mean)
   fDbar<-apply(sweep(discards.sel(fleet,1,1),c(2,6),fleet@effort,"*")[ac(2:4)],c(2,6),mean)
   fCbar<-fLbar+fDbar

   res=FLQuants(revenues=revenues,
                #costs   =costs,
                #profits =profits,
                yield   =yield,
                catch   =catch,
                discards=discards,
                rec     =rec,
                ssb     =ssb,
                biomass =biomass,
                mnsz    =mnsz,
                pmat    =pmat,
                harvest =harvest,
                fLbar    =fLbar,
                fDbar    =fDbar,
                fCbar    =fCbar)
   return(res)
   }

perfectAssess<-function(stk,biol,fleet,yrs.assess,bias=0,cv=0)
    {
    stk        <-stk[,yrs.assess]
    stk@stock.n<-setPlusGroup(biol@n[,yrs.assess],stk@range["plusgroup"])
    stk@harvest<-setPlusGroup(calcF(m(biol)[,yrs.assess],
                              catch.n(fleet,1,1)[,yrs.assess],
                              n(biol)[,yrs.assess]),stk@range["plusgroup"])

    if (bias!=0)
       {
       stk@stock.n[,yrs.assess[length(yrs.assess)]]<-stk@stock.n[,yrs.assess[length(yrs.assess)]]*(1-bias)
       stk@harvest[,yrs.assess[length(yrs.assess)]]<-stk@harvest[,yrs.assess[length(yrs.assess)]]*(1-bias)
       }

    if (cv>0)
       {
       stk@stock.n[,yrs.assess[length(yrs.assess)]]<-stk@stock.n[,yrs.assess[length(yrs.assess)]]*exp(rnorm(0,cv))/(cv^2)
       stk@harvest[,yrs.assess[length(yrs.assess)]]<-stk@harvest[,yrs.assess[length(yrs.assess)]]*exp(rnorm(0,cv))/(cv^2)
       }

    return(stk)
    }

SmryPlots<-function(fleet,biol,OMf,stk,yrs,its=c(1,7,9))
     {
 #stk<-MPStk;OMf<-res$f;yrs=2000:(iyr)
 #     fleet<-window(fleet,end=max(yrs)+2)
#     biol <-window(biol, end=max(yrs)+2)
      stk  <-window(stk,  end=max(yrs)+2)

     if (is.numeric(yrs)) yrs<-as.character(min(yrs):(max(yrs)+2))

     c.fleet =apply(catch.n(fleet,1)[[1]]*catch.wt(fleet,1)[[1]],c(2,6),sum)

     #OMf    <-calcF(m(biol),catch.n(fleet,1,1),n(biol))

     fbarRng<-ac(range(biol,"minfbar"):range(biol,"maxfbar"))
     OM.f    =quantile(apply(OMf[fbarRng],         c(2,6),mean),prob=c(0.25,0.50,0.75),na.rm=T)[,yrs]
     fbarRng<-ac(range(stk,"minfbar"):range(stk,"maxfbar"))
     MP.f    =quantile(apply(stk@harvest[fbarRng],c(2,6),mean), prob=c(0.25,0.50,0.75),na.rm=T)[,yrs]
     OM.rec  =quantile(biol@n[1,],                              prob=c(0.25,0.50,0.75),na.rm=T)[,yrs]
     MP.rec  =quantile(stk@stock.n[1,],                         prob=c(0.25,0.50,0.75),na.rm=T)[,yrs]
     OM.ssb  =quantile(ssb(biol),                               prob=c(0.25,0.50,0.75),na.rm=T)[,yrs]
     MP.ssb  =quantile(ssb(stk),                                prob=c(0.25,0.50,0.75),na.rm=T)[,yrs]
     OM.yield=quantile(c.fleet,                                 prob=c(0.25,0.50,0.75),na.rm=T)[,yrs]
     MP.yield=quantile(computeCatch(stk),                       prob=c(0.25,0.50,0.75),na.rm=T)[,yrs]

     worm.f    =apply(OMf[fbarRng],  c(2,6),mean)[,yrs,,,,its]
     worm.rec  =biol@n[1,  yrs,,,,its]
     worm.ssb  =ssb(biol)[,yrs,,,,its]
     worm.yield=c.fleet[,  yrs,,,,its]

     dmns     <-dimnames(OM.ssb)
     dmns$iter<-1:10
     ssb.  <-FLQuant(c(OM.ssb,  MP.ssb,  worm.ssb  ), dimnames=dmns)
     f.    <-FLQuant(c(OM.f,    MP.f,    worm.f    ), dimnames=dmns)
     rec.  <-FLQuant(c(OM.rec,  MP.rec,  worm.rec  ), dimnames=dmns)
     yield.<-FLQuant(c(OM.yield,MP.yield,worm.yield), dimnames=dmns)

     print(xyplot(data~year|qname,groups=iter,data=FLQuants(F        =f.,
                                                            SSB      =ssb.,
                                                            Recruits =rec.,
                                                            Yield    =yield.),
                                                            type="l",col=c(rep("red",3),rep("blue",3),rep("cornsilk3",3)),lwd=c(1,2,1,1,2,1,1,1,1),scale="free"))
     }

SmryPlot1<-function(fleet,biol,OMf,stk,yrs,its=1)
     {
     stk  <-window(stk,  end=max(yrs)+2)

     if (is.numeric(yrs)) yrs<-as.character(min(yrs):(max(yrs)+2))

     c.fleet =apply(catch.n(fleet,1)[[1]]*catch.wt(fleet,1)[[1]],c(2,6),sum)

     fbarRng<-ac(range(biol,"minfbar"):range(biol,"maxfbar"))
     MP.f    =apply(OMf[fbarRng],         c(2,6),mean)[,yrs,,,,its]
     MP.rec  =rec(stk)[, yrs,,,,its]
     MP.ssb  =ssb(stk)[,yrs,,,,its]
     MP.yield=computeCatch(stk)[,  yrs,,,,its]

     fbarRng<-ac(range(stk,"minfbar"):range(stk,"maxfbar"))
     OM.f    =apply(stk@harvest[fbarRng], c(2,6),mean)[,yrs,,,,its]
     OM.rec  =biol@n[1,  yrs,,,,its]
     OM.ssb  =ssb(biol)[,yrs,,,,its]
     OM.yield=c.fleet[,  yrs,,,,its]

     dmns     <-dimnames(OM.ssb)
     dmns$iter<-1:10
     ssb.  <-FLQuant(c(OM.ssb,  MP.ssb  ), dimnames=dmns)
     f.    <-FLQuant(c(OM.f,    MP.f    ), dimnames=dmns)
     rec.  <-FLQuant(c(OM.rec,  MP.rec  ), dimnames=dmns)
     yield.<-FLQuant(c(OM.yield,MP.yield), dimnames=dmns)

     print(xyplot(data~year|qname,groups=iter,data=FLQuants(F        =f.,
                                                            SSB      =ssb.,
                                                            Recruits =rec.,
                                                            Yield    =yield.),
                                                            type="l",col=c("red","blue"),lwd=1,scale="free"))
     }

