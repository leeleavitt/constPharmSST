#' Function to display a single calcium imaging trace image and other data
PeakFunc7 <- function(dat,n.names,t.type="t.dat",Plotit.trace=T,Plotit.both=F, info=T,lmain=NULL, bcex=.7, yvar=T, ylim.max=NULL, zf=40, pts=T, lns=T, levs=NULL, underline=T, dat.n=""){
    dat.name<-deparse(substitute(dat))
    if(! dat.name == ""){dat.name<-dat.n
    }else{dat.name<-dat.name}	
    
    if(is.null(lmain)){
        lmain=n.names
    }else{lmain=lmain}
    if(class(t.type)=="character")
    {
        dat.select<-t.type
        dat.t<-dat[[dat.select]]
    }else{
        dat.select<-select.list(names(dat), multiple=T)
        dat.t<-dat[[dat.select]]
    }

    if(yvar){
        ymax<-max(dat.t[,n.names])*1.05
        ymin<-min(dat.t[,n.names])*.95
        yrange<-ymax-ymin
    }else{		
        if(is.null(ylim.max)){ylim.max<-1.4}else{ylim.max<-ylim.max}
        if(Plotit.trace){ylim <- c(-.1,ylim.max)}
        if(Plotit.both){ylim <- c(-.5,ylim.max)}
        ymin<-min(ylim)
        ymax<-max(ylim)
        yrange<-ymax-ymin
    }

    if(Plotit.trace){ylim <- c(ymin,ymax)}
    if(Plotit.both){ymin<- -.5;ylim <- c(ymin,ymax)}
    par(xpd=FALSE)
    xlim <- range(dat.t[,1]) # use same xlim on all plots for better comparison
    
    #   ylim <- range(intensity(s1))
    if(is.null(levs)){levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    }else{levs<-levs}
    par(mar=c(9,6.2,3.5,13), bty="n")
    plot(dat.t[,n.names]~dat.t[,1], main=lmain,xlim=xlim,ylim=ylim,xlab="", ylab="",pch="", cex=.5)

    #axis(3,tick=TRUE, outer=F )
    axis(1, at= seq(0, max(dat.t[,1]),10), tick=TRUE)
    
    # Tool for labeling window regions
    wr<-dat$w.dat[,"wr1"]
    #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
    x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs]
    x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs]
    y1s <- rep(par("usr")[4],length(x1s))
    y2s <- rep(par("usr")[3],length(x1s))
    rect(x1s,y1s,x2s,y2s,col="grey95")
    
    
    # Tool for labeling cellular aspects, gfp.1, gfp.2, tritc, area 
    legend(x=par("usr")[1]-xinch(1.45), y=par("usr")[3]-yinch(.25), xpd=TRUE, inset=c(0,-.14),bty="n", cex=.7, legend=c(
        if(!is.null(dat$c.dat[n.names, "CGRP"])){paste("CGRP","",round(dat$c.dat[n.names,"CGRP"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp"])){paste("GFP","",round(dat$c.dat[n.names,"mean.gfp"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.1"])){paste("GFP.1","",round(dat$c.dat[n.names,"mean.gfp.1"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.2"])){paste("GFP.2","",round(dat$c.dat[n.names,"mean.gfp.2"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.start"])){paste("mean.gfp.start","",round(dat$c.dat[n.names,"mean.gfp.start"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.end"])){paste("mean.gfp.end","",round(dat$c.dat[n.names,"mean.gfp.end"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.gfp.immuno"])){paste("CGRP immunostain","",round(dat$c.dat[n.names,"mean.gfp.immuno"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.dapi"])){paste("DAPI","",round(dat$c.dat[n.names,"mean.dapi"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "IB4"])){paste("IB4","",round(dat$c.dat[n.names,"IB4"],digits=4))},
        if(!is.null(dat$c.dat[n.names, "mean.tritc"])){paste("IB4","",round(dat$c.dat[n.names, "mean.tritc"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.start"])){paste("IB4.start","",round(dat$c.dat[n.names, "mean.tritc.start"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.end"])){paste("IB4.end","",round(dat$c.dat[n.names, "mean.tritc.end"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.tritc.immuno"])){paste("NF200 immunostain","",round(dat$c.dat[n.names, "mean.tritc.immuno"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.cy5.start"])){paste("IB4.start","",round(dat$c.dat[n.names, "mean.cy5.start"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "mean.cy5.end"])){paste("IB4.end","",round(dat$c.dat[n.names, "mean.cy5.end"], digits=4))}, 
        if(!is.null(dat$c.dat[n.names, "area"])){paste("area","", round(dat$c.dat[n.names, "area"], digits=4))},
        if(!is.null(dat$c.dat[n.names, "ROI.Area"])){paste("area","", round(dat$c.dat[n.names, "ROI.Area"], digits=4))},
        #if(!is.null(dat$c.dat[n.names, "perimeter"])){paste("perimeter","", round(dat$c.dat[n.names, "perimeter"], digits=0))},
        if(!is.null(dat$c.dat[n.names, "circularity"])){paste("circularity","", round(dat$c.dat[n.names, "circularity"], digits=4))}
        )
    )
    
    legend(x=par("usr")[2]+xinch(.8), y=par("usr")[3]-yinch(.9), xpd=TRUE, inset=c(0,-.14), bty="n", cex=.7, legend=dat.name)

    
    #Adding binary scoring for labeling to plot
    par(xpd=TRUE)
    if(!is.null(dat$bin[n.names, "gfp.bin"])){text(y=par("usr")[4]+yinch(.5), x=par("usr")[2]+xinch(1.8), paste("GFP:",dat$bin[n.names,"gfp.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "tritc.bin"])){text(y=par("usr")[4]+yinch(.25), x=par("usr")[2]+xinch(1.8), paste("IB4 :",dat$bin[n.names,"tritc.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "cy5.bin"])){text(y=par("usr")[4]+yinch(.25), x=par("usr")[2]+xinch(1.8), paste("IB4 :",dat$bin[n.names,"cy5.bin"]), cex=.7)}
    if(!is.null(dat$bin[n.names, "drop"])){text(y=par("usr")[4]+yinch(0), x=par("usr")[2]+xinch(1.8), paste("Drop :",dat$bin[n.names,"drop"]), cex=.7)}


    # Tool for lableing window region information
    levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
    if(info){
        x.name<-n.names
        #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])), "")
        mtext(c("max","snr"), side=3, at=-max(dat.t[,1])*.05, line=c(0, .7), cex=.6)
        for(i in 1:length(levs)){
            max.name<-paste(levs[i],".max", sep="")
            max.val<-round(dat$scp[x.name, max.name], digits=3)
            mtext(max.val, side=3, at=levs.loc[ levs[i] ], line=0, cex=.6)
            
            tot.name<-paste(levs[i],".snr", sep="")
            tot.val<-round(dat$scp[x.name, tot.name], digits=3)
            mtext(tot.val, side=3, at=levs.loc[ levs[i] ], line=.7, cex=.6)
        }
        
    # Tool for labeling the binary score
        #levs<-setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
        z<-t(dat$bin[n.names,levs])
        zz<-z==1
        zi<-attributes(zz)
        zzz<-which(zz, arr.ind=T)
        #levs<-zi$dimnames[[2]][zzz[,2]]
        levs1<-unique(as.character(row.names(zzz)))
        x1s <- tapply(dat$w.dat[,"Time"],as.factor(wr),min)[levs1]
        x2s <- tapply(dat$w.dat[,"Time"],as.factor(wr),max)[levs1]
        y1s <- rep(par("usr")[4],length(x1s))
        y2s <- rep(par("usr")[3],length(x1s))
        rect(x1s,y1s,x2s,y2s,col="grey80")
        #levs <- setdiff(unique(wr),"")
    }
    
    #text(dat.t[match(levs,wr),"Time"],c(ymin, ymin+(yrange*.2)),levs,pos=4,offset=0,cex=bcex)	
    #text(dat.t[match(levs,wr),"Time"],par("usr")[3],levs,pos=3,offset=-4.2,cex=bcex, srt=90)    
    levs_cex <- nchar(levs)
	levs_cex[ levs_cex <= 12*1.3  ] <- 1
	levs_cex[ levs_cex > 12*1.3  ] <- 12/levs_cex[ levs_cex>12*1.3  ]*1.3

    text(levs.loc,par("usr")[3],levs,pos=3,offset=-4.3,cex=levs_cex, srt=90)	

    if(Plotit.both){
        if(!is.null(dat$der)){lines(dat$der[,n.names]~dat.t[-1,1], lwd=.01, col="paleturquoise4")}
        par(xpd=T)
        abline(h=0)
        if(lns){lines(dat.t[,n.names]~dat.t[,1])
        }else{}
        if(pts){points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        }else{}
        par(xpd=F)
    }
    
    if(Plotit.trace){
        par(xpd=T)
        if(lns){lines(dat.t[,n.names]~dat.t[,1])
        }else{}
        
        if(pts){points(dat.t[,n.names]~dat.t[,1], pch=16, cex=.3)
        }else{}
        
        par(xpd=F)
    }
    
    ##Tool for adding underline to plot
    if(underline){
        par(xpd=F)
        abline(h=min(dat.t[,n.names]), col="black")
        par(xpd=T)
    }else{}

    
    ## Tool for adding rasterImages to plot
    
    ###Finding the picture loaction of the cells
    if(!is.null(dat$img1)){
        if(is.null(zf)){zf<-20
        }else{zf<-zf}

        img.dim<-dim(dat$img1)[1]
        x<-dat$c.dat[n.names,"center.x"]
        left<-x-zf
        if(left<=0){left=0; right=2*zf}
        right<-x+zf
        if(right>=img.dim){left=img.dim-(2*zf);right=img.dim}
        
        y<-dat$c.dat[n.names,"center.y"]
        top<-y-zf
        if(top<=0){top=0; bottom=2*zf}
        bottom<-y+zf
        if(bottom>=img.dim){top=img.dim-(2*zf);bottom=img.dim}
        
        par(xpd=TRUE)
    }
    ### Where to plot pictures
    #ymax<-max(dat.t[,n.names])*1.05
    #ymin<-min(dat.t[,n.names])*.95
    #yrange<-ymax-ymin
    
    ymax<-par("usr")[4]
    xmax<-par("usr")[2]
    if(!is.null(dat$img1)){
        img1<-dat$img1
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax+yinch(.8)
        ybottom<-ymax
        tryCatch(
            rasterImage(img1[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img1[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img2)){
        img2<-dat$img2
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax+yinch(.8)
        ybottom<-ymax
        tryCatch(
            rasterImage(img2[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img2[top:bottom,left:right],xleft,ybottom,xright,ytop))

    }

    if(!is.null(dat$img3)){
        img3<-dat$img3
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        tryCatch(
            rasterImage(img3[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img3[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img4)){
        img4<-dat$img4
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax
        ybottom<-ymax-yinch(.8)
        tryCatch(
            rasterImage(img4[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img4[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }

    if(!is.null(dat$img5)){
        img5<-dat$img5
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        tryCatch(
            rasterImage(img5[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img5[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img6)){
        img6<-dat$img6
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax-yinch(.8)
        ybottom<-ymax-yinch(1.6)
        tryCatch(
            rasterImage(img6[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img6[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img7)){
        img7<-dat$img7
        xleft<-xmax
        xright<-xmax+xinch(.8)
        ytop<-ymax-yinch(1.6)
        ybottom<-ymax-yinch(2.4)
        tryCatch(
            rasterImage(img7[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img7[top:bottom,left:right],xleft,ybottom,xright,ytop))
    }
    
    if(!is.null(dat$img8)){
        img8<-dat$img8
        xleft<-xmax+xinch(.8)
        xright<-xmax+xinch(1.6)
        ytop<-ymax-yinch(1.6)
        ybottom<-ymax-yinch(2.4)
        tryCatch(
            rasterImage(img8[top:bottom,left:right,],xleft,ybottom,xright,ytop),
            error=function(e) rasterImage(img8[top:bottom,left:right],xleft,ybottom,xright,ytop))

    }
}	
