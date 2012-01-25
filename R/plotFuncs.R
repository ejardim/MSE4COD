getQuartiles<-function(x,yrs,quant,iters=c(33,17,39,41,89),scale=TRUE){
  x<-x[x[,"year"] %in% yrs,]

  qile<-rbind(cbind(series=.25,year=yrs,data=tapply(x[,quant],x[,c("year")],quantile,prob=c(0.25))),
              cbind(series=.50,year=yrs,data=tapply(x[,quant],x[,c("year")],quantile,prob=c(0.50))),
              cbind(series=.75,year=yrs,data=tapply(x[,quant],x[,c("year")],quantile,prob=c(0.75))))

  xIters<-x[x[,"iter"] %in% iters ,c("iter","year",quant)]

  names(xIters)<-c("series","year","data")

  qile<-rbind(qile,xIters)

  if (scale)
    qile[,"data"]<-qile[,"data"]/qile[qile[,"year"]==2007 & qile[,"series"]==0.5,"data"]

  return(qile)
  }

plotBeer=function(v){
	usr=par("usr")
	ymin=usr[3]
	ymax=usr[4]
	yrng=ymax-ymin
	ymin1=usr[3]+.1*yrng
	ymax1=usr[4]-.2*yrng

	yrng1=ymax1-ymin1
	ymax2=ymin1+abs(v)*yrng1
	xmid=(usr[2]+usr[1])/2
	ymid=(ymax1+ymin1)/2
	xrng=(usr[2]-usr[1])
	xpoly=c(xmid-.15*xrng,xmid-(.15+abs(v)*.1)*xrng,xmid+(.15+abs(v)*.1)*xrng,xmid+.15*xrng,xmid-.15*xrng)
	ypoly=c(ymin1,ymax2,ymax2,ymin1,ymin1)
	polygon(xpoly,ypoly,col="gold",border="burlywood")
	bubblex=runif(round(500*abs(v),0),xmid-(.15+abs(v*.95)*.1)*xrng,xmid+(.15+abs(v*.95)*.1)*xrng)
	bubbley=runif(round(500*abs(v),0),ymax2-.02*yrng1,ymax2+.02*yrng1)
	points(bubblex,bubbley,pch=21,col = "gold", bg = "white",cex=seq(0.1,1,length=10))
	points(c(xmid-.15*xrng,xmid+.15*xrng),c(ymin1,ymin1),type="l",lwd=4)
	points(c(xmid-.15*xrng,xmid-.25*xrng),c(ymin1,ymax1),type="l",lwd=4)
	points(c(xmid+.15*xrng,xmid+.25*xrng),c(ymin1,ymax1),type="l",lwd=4)
	if(v<0){
		text(xmid,ymid,labels=c(paste("-",abs(v))),cex=.5+abs(v)*2)
	}else{
		text(xmid,ymid,labels=c(paste("+",abs(v))),cex=.5+abs(v)*2)
    }
  }

plotHist<-function(x,brp.,title="",refYr,axs=FALSE)
   {
   yrFlag <-x[,"year"]<=ac(refYr)
   hcrFlag<-x[yrFlag,"HCR"]=="EU"
   xTrk.1 <-tapply(x[yrFlag &  hcrFlag,"ssb"],    x[yrFlag &  hcrFlag,"year"],quantile,prob=0.5)/c(brp.[    "ssb"])
   yTrk.1 <-tapply(x[yrFlag &  hcrFlag,"harvest"],x[yrFlag &  hcrFlag,"year"],quantile,prob=0.5)/c(brp.["harvest"])
   xTrk.2 <-tapply(x[yrFlag & !hcrFlag,"ssb"],    x[yrFlag & !hcrFlag,"year"],quantile,prob=0.5)/c(brp.[    "ssb"])
   yTrk.2 <-tapply(x[yrFlag & !hcrFlag,"harvest"],x[yrFlag & !hcrFlag,"year"],quantile,prob=0.5)/c(brp.["harvest"])

   yrFlag<-x[,"year"]==ac(refYr)
   xPts  <-x[yrFlag,    "ssb"]/c(brp.[    "ssb"])
   yPts  <-x[yrFlag,"harvest"]/c(brp.["harvest"])
   colPts<-x[yrFlag,"HCR"]

   plotTapas(xTrk.1,yTrk.1,xTrk.2,yTrk.2,xPts,yPts,"sienna",colPts=colPts,xlab=expression(SSB/Bpa),ylab=expression('Fishing Mortality'),title=title,maxY=1,maxX=5.0,axs)
   }

plotTapas<-function(xTrk.1,yTrk.1,xTrk.2,yTrk.2,xPts,yPts,colTrk,colPts,xlab="X",ylab="Y",title="",maxX="missing",maxY="missing",year,axs=TRUE)
    {
    fish.pg<-function(maxX,maxY)
        {
        polygon(x=c(-0.5,1,1,-0.5),            y=c(0.4,0.4,maxY+.5,maxY+.5),    col="red2")
        polygon(x=c(1.0,maxX+0.5,maxX+0.5,1.0),y=c(-0.5,-0.5,0.4,0.4),          col="lightgreen")
        polygon(x=c(-0.5,1,1,-0.5),            y=c(-0.5,-0.5,0.4,0.4),          col="lightgoldenrod1")
        polygon(x=c(1.0,maxX+0.5,maxX+0.5,1.0),y=c(0.4,0.4,maxY+.5,maxY+.5),    col="lightgoldenrod1")
        }
    if (missing(maxX))
       maxX<-max(xTrk.1,xTrk.2,xPts,na.rm=T)
    else
       maxX<-max(maxX, xTrk.1,xTrk.2,xPts,na.rm=T)
    if (missing(maxY))
       maxY<-max(yTrk.1,yTrk.2,yPts,na.rm=T)
    else
       maxY<-max(maxY,yTrk.1,yTrk.2,yPts,na.rm=T)

    plot(yTrk.1~xTrk.1, col=colTrk, type="l",xlab=xlab,ylab=ylab,main=paste(title),xlim=c(0,5),ylim=c(0,maxY),axes=axs)
    fish.pg(maxX,maxY)
    abline(h=0.4,v=1,col="grey")

    #points(yPts[colPts=="EU"]~xPts[colPts=="EU"], col="grey45",pch=19,cex=0.75)
    #points(yPts[colPts!="EU"]~xPts[colPts!="EU"], col="blue",pch=19, cex=0.75)
    #lines( yTrk.1~xTrk.1,col="grey12",lwd=2, type="b")
    #lines( yTrk.2~xTrk.2,col="blue", lwd=2, type="b")
    
    points(yPts[colPts=="EU"]~xPts[colPts=="EU"], col="black",pch=1,cex=0.75)
    points(yPts[colPts!="EU"]~xPts[colPts!="EU"], col="blue",pch=5, cex=0.75)
    
    points( yTrk.1[length(yTrk.1)]~xTrk.1[length(xTrk.1)],col="black",pch=19, cex=2.0)
    points( yTrk.2[length(yTrk.2)]~xTrk.2[length(xTrk.2)],col="black",bg="blue", pch=23, cex=2.0)    
    }

plotSmry<-function(x,yrs,OM,OEM,SR,minTAC,Title,fileNm,axs=TRUE){     #this produces phase-plots (see AGCREMP report)
    x.<-x[x[,"OEM"]==OEM & x[,"OM"]==OM & x[,"SR"]==SR & x[,"minTAC"]==minTAC & !is.na(x[,"harvest"]) & !is.na(x[,"ssb"]),]
    for (iYr in yrs)
          {
          plotHist(x.,brp.,title=iYr,iYr,axs)

          resultsEU<-x.[x.[,"HCR"]=="EU",c("year","harvest","ssb")]
          grnS<-numeric(dim(resultsEU[resultsEU[,"year"]==iYr,])[1])
          grnF<-grnS
          grnS[resultsEU[resultsEU[,"year"]==iYr,"ssb"]    >=Bpa]<-1
          grnF[resultsEU[resultsEU[,"year"]==iYr,"harvest"]<=0.4]<-1
          EU<-table(grnS,grnF)/250

          resultsNor<-x.[x.[,"HCR"]=="Norway",c("year","harvest","ssb")]
          grnS<-numeric(dim(resultsNor[resultsNor[,"year"]==iYr,])[1])
          grnF<-grnS
          grnS[resultsNor[resultsNor[,"year"]==iYr,"ssb"]    >=Bpa]<-1
          grnF[resultsNor[resultsNor[,"year"]==iYr,"harvest"]<=0.4]<-1
          Nor<-table(grnS,grnF)/250

          if (!("0" %in% dimnames(EU )[[1]] & "0" %in% dimnames(EU )[[2]])) redEU <-as.integer(0) else redEU <-as.integer(EU[ "0","0"]*100)
          if (!("0" %in% dimnames(Nor)[[1]] & "0" %in% dimnames(Nor)[[2]])) redNor<-as.integer(0) else redNor<-as.integer(Nor["0","0"]*100)
          text(-0.2,1,     label=paste("EU",    format(redEU, 2),"%"), col="grey45", cex=1.25, pos=4)
          text(-0.2,0.92,  label=paste("Norway",format(redNor,2),"%"), col="blue",   cex=1.25, pos=4)

          if (!("1" %in% dimnames(EU )[[1]] & "1" %in% dimnames(EU )[[2]])) grnEU <-as.integer(0) else grnEU <-as.integer(EU[ "1","1"]*100)
          if (!("1" %in% dimnames(Nor)[[1]] & "1" %in% dimnames(Nor)[[2]])) grnNor<-as.integer(0) else grnNor<-as.integer(Nor["1","1"]*100)
          text(2.7,0,  label=paste("EU",     format(grnEU, 2), "%"), col="grey45", cex=1.25, pos=4)
          text(3.7,0,  label=paste("Norway", format(grnNor,2), "%"), col="blue",   cex=1.25, pos=4)

          savePlot(paste(fileNm,"_",iYr,".png",sep=""),type="png")
          }
     }

plotRasta<-function(x){
    ssbTrk<-ssb( Stk)/c(refpts(bftBRP.99.90s)[1,"ssb",    "fmax",drop=T])
    fbrTrk<-fbar(Stk)/c(refpts(bftBRP.99.90s)[1,"harvest","fmax",drop=T])

    red   <-FLQuant(0,dimnames=dimnames(ssbTrk))
    green <-red
    yellow<-red

    red[   ssbTrk< 1 & fbrTrk> 1]<-1
    green[ ssbTrk>=1 & fbrTrk<=1]<-1
    yellow[red   ==0 & green ==0]<-1

    print(xyplot(data~year,groups=qname,data=lapply(FLQuants(green=green,red=red,yellow=yellow),function(x) apply(x,2,mean)),
               col=c("green","red","yellow"),lwd=3,type="l",xlab="Year",ylab="Probability of being in the Zone",xlim=c(2007,2023)))
    }
