#Software to display the transcriptome next to traces
LinesEvery_ts <- function(exps = exps, exp_info = exp_info, reduce=F, levs=NULL, pdf=T, cell_info = NULL, imgs = NULL){
    
    if(is.null(cell_info)){
        cell_info <- c("area","mean.gfp","mean.cy5","cDNA.conc.ug.mL", "M.Assigned", "qpcr" , genes)
    }
    
    if(is.null(imgs)){
        imgs <- c("img1", "img3", "img4")
    }
    
    cat("\nEnter the names of genes or gene families you would like to observe next to the traces\n")    
    genes <- scan(what='character')
    #This setof functions allow you to select dimentions to view
    genes_tot <- c()
    for(i in 1:length(genes)){
        genes_tot<-c(genes_tot, grep(genes[i], colnames(tmp_rd$c.dat), ignore.case=T, value=T))
    }

    cat("\nPlease reduce your selection\n")
    if(reduce){
        genes <-select.list(genes_tot, multiple=T)
    }else{
        genes <- genes_tot
    }

    #Now we will plot these cells of interest 
    if(pdf){
        pdf("default.pdf", width=15, height= 6*length(exps))
    }else{    
        dev.new(width=20, height= 8*length(exps) )
    }    
    par(mfrow = c(length(exps), 1) )

    for( i in 1:length(exps)){
        exp_i_info <- exp_info[ exp_info$rd.name == exps[i], ]
        exp_i_info <- exp_i_info[order(exp_i_info$label_cellType),]
        exp_cells <- exp_i_info$Cell.name
        tmp_rd <- get(exps[i])
        LinesEvery.6(tmp_rd, exp_cells, values = cell_info, t.type="blc", dat.n=exps[i], plot.new=F, lw=2, underline=F, levs=levs, img=imgs, gnomex_lab = tmp_rd$c.dat[,'Gnomex.Label',drop=F])
        box('figure', lwd=3)
    }   
    if(pdf){dev.off()}  
}

LinesEvery.6 <- function(dat,m.names, img="img1",channel=NULL,pic.plot=TRUE,zf=NULL, t.type="mp.1", snr=NULL,lmain="",cols="black", levs=NULL, levs.cols="grey90", values=NULL,plot.new=T,sf=1,lw=1,bcex=1,p.ht=7,p.wd=10, lns=T, pts=F, underline=T,dat.n=NULL, gnomex_lab=NA){
    #require(RColorBrewer)
    dat.name<-deparse(substitute(dat))
    if(dat.name=="dat" | dat.name == "tmp.rd" | dat.name == "tmp_rd"){
        dat.name<-dat.n
    }else{
        dat.name<-dat.name
    }
    
    #Trace Selector if t.type is empty.  t.type must be character input
    if(class(t.type)=="character"){
        t.dat<-dat[[t.type]]# if trace type is empty select the data, you would like your trace to be
    }else{
        t.type<-menu(names(dat));t.dat<-dat[[t.type]]
    }
    
    m.names <- intersect(m.names,names(t.dat))
    xseq <- t.dat[,1]
    #upper ylimit 
    hbc <- length(m.names)*sf+max(t.dat[,m.names])

    #Selecting multiple images
    if(is.null(img)){
        img.l<-select.list(grep("img",names(dat), value=T), multiple=T)
    }else{
        img.l<-img
    }

    if(length(m.names) > 0){
        #For pdf output
            #if(is.null(pdf.name))
            #	{dev.new(width=14,height=8)}
            #else
                #{if(length(grep("\\.pdf",pdf.name))>0){pdf(pdf.name,width=p.wd,height=p.ht)}else{png(pdf.name,width=1200,height=600)}}
        ## Tool for addind value tags displayed on the right side of trace
        #See line 3016 for where values come into play
        #values<-c("area", "mean.gfp.start", "mean.gfp.end" "mean.tritc.start", "mean.tritc.end")
            if(is.null(values)){
                values<-c("area")
            }else{values<-values}

        ## Tool for color labeleing
        ## Tool for single color labeling
            if(cols=="brew.pal"){
                #cols <- rainbow(length(m.names),start=.55)
                require(RColorBrewer)
                cols <-brewer.pal(8,"Dark2")
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            }
            if(cols=="rainbow"){
                cols<-rainbow(length(m.names),start=.7,end=.1)
            }
            if(cols=="topo"){
                cols<-topo.colors(length(m.names))
            }else{		
                cols<-cols
                cols <- rep(cols,ceiling(length(m.names)/length(cols)))
                cols <- cols[1:length(m.names)]
            }
            if(plot.new){
                if(length(m.names)>10){dev.new(width=10+length(img)+(length(values)*.6),height=12)}
                else(dev.new(width=10+length(img)+(length(values)*.6),height=8))
            }
        ## Begin plotting
            xinch(length(img))
            par(xpd=FALSE,mai=c(1.2, .5+(.5*length(img.l)), .5, 0.4*length(values)), bty="l")
            plot(xseq,t.dat[,m.names[1]],ylim=c(0,hbc),xlab="",main=lmain,type="n",yaxt="n",xlim=c(min(xseq)-1.5,max(xseq)), ylab="")#-sf
            #axis(1, at=seq(floor(min(t.dat[,1])),ceiling(max(t.dat[,1])), 1))
            par(xpd=T)
            text(rep(0,length(m.names))-xinch(.3),seq(1,length(m.names))*sf+t.dat[1,m.names],paste0(m.names,"_", gnomex_lab[m.names,]), cex=.5,col=cols,pos=3)
            par(xpd=F)
        ## Tool for adding window region labeling
            if(is.null(levs)){levs <- setdiff(unique(as.character(dat$w.dat[,"wr1"])),"")
            }else{levs<-levs}
            wr<-dat$w.dat$wr1
            if(length(wr) > 0){
                x1s <- tapply(dat$w.dat[,1],as.factor(wr),min)[levs]
                x2s <- tapply(dat$w.dat[,1],as.factor(wr),max)[levs]
                y1s <- rep(par("usr")[3],length(x1s))
                y2s <- rep(par("usr")[4],length(x1s))
                rect(x1s,y1s,x2s,y2s,col=levs.cols,border="black")
                par(xpd=TRUE)
                #text(x1s-xinch(.1),par("usr")[3]-yinch(1),levs,cex=.8*bcex, srt=90)	
                #dat$t.dat[match(levs,wr),"Time"]
                levs.loc<-tapply(dat$w.dat[,"Time"],as.factor(wr),mean)[levs]
				levs_cex <- nchar(levs)
                zoom_fac <- .7
				levs_cex[ levs_cex <= 12*zoom_fac ] <- 1
				levs_cex[ levs_cex > (12*zoom_fac) ] <- 12/levs_cex[ levs_cex>(12*zoom_fac) ]*zoom_fac
                text(levs.loc,par("usr")[3]+xinch(.05),levs,pos=3,offset=-4.3,cex=levs_cex, srt=90)	
                par(xpd=FALSE)
            }
        
        ## Tool for adding line, point and picture to the plot
            for(i in 1:length(m.names)){
                ypos<-t.dat[,m.names[i]]+i*sf
                if(lns){lines(xseq,ypos, lty=1,col=cols[i],lwd=lw)}
                if(pts){points(xseq,ypos,pch=16,col=cols[i],cex=.3)}
                if(!is.null(snr)){
                    pp1 <- snr[,m.names[i]] > 0 & is.element(wr,levs)
                    pp2 <- snr[,m.names[i]] > 0 & !is.element(wr,levs)
                    points(xseq[pp1],t.dat[pp1,m.names[i]]+i/10,pch=1,col=cols[i])
                    points(xseq[pp2],t.dat[pp2,m.names[i]]+i/10,pch=0,col=cols[i])
                }
                if(underline){abline(h=min(ypos), col="black")}else{}
            }
            par(xpd=TRUE)
            
        ## Tool for adding Value info on right side of trace
            placement<-seq(0,length(values),.4)
            digits<-c(0,rep(4,length(values)))
            text( max(xseq) + xinch(placement[1:length(values)]), par("usr")[4]+yinch(.2), pos=4, values ,cex=bcex*.75, srt=30)
            
            for(i in 1:length(values)){
                if(!is.null(dat$c.dat[m.names, values[i]])){
                    rtag<-values[i]
                    rtag <- round(dat$c.dat[m.names,rtag], digits=1)
                    new_max <- .8
                    new_min <- .5
                    range02 <- function(x){ (x - min(x))/(max(x)-min(x)) * (new_max - new_min) + new_min }
                    rtag_cex <- range02(rtag)

                    text(
                        rep(max(xseq)+xinch(placement[i]),length(m.names)),
                        seq(1,length(m.names))*sf + t.dat[nrow(t.dat),m.names],
                        paste(rtag),
                        #cex=.6*bcex,
                        cex=rtag_cex,
                        col=cols,
                    pos=4)
                }
            }
            
        ##Tool for adding images to the left side of the plot
        if(is.null(zf)){
            zf<-20
        }else{zf<-zf}
        
        pic.pos<-list()
        for(i in 1:length(m.names)){
            ypos<-t.dat[1,m.names[i]]+i*sf
            pic.pos[[i]]<-ypos
        }

        xinchseq1<-seq(1,5,.5)
        xinchseq2<-seq(.5,5,.5)
        
        if(is.null(channel)){channel<-rep(list(c(1:3)),length(img.l))
        }else{channel<-channel}

        for(j in 1:length(img.l)){
            for(i in 1:length(m.names)){		
                img.dim<-dim(dat$img1)[1]
                x<-dat$c.dat[m.names[i],"center.x"]
                left<-x-zf
                if(left<=0){
                    left=0
                    right=2*zf
                }
                
                right<-x+zf
                if(right>=img.dim){
                    left=img.dim-(2*zf)
                    right=img.dim
                }
                    
                y<-dat$c.dat[m.names[i],"center.y"]
                top<-y-zf
                if(top<=0){
                    top=0
                    bottom=2*zf
                }
                bottom<-y+zf
                
                if(bottom>=img.dim){
                    top=img.dim-(2*zf)
                    bottom=img.dim
                }
                
                par(xpd=TRUE)
                xleft<-min(dat$t.dat[,1])-xinch(xinchseq1[j])
                xright<-min(dat$t.dat[,1])-xinch(xinchseq2[j])
                ytop<-pic.pos[[i]]+yinch(.25)
                ybottom<-pic.pos[[i]]-yinch(.25)
                
                tryCatch(
                    rasterImage(dat[[img.l[j]]][top:bottom,left:right,channel[[j]]],xleft,ybottom,xright,ytop),
                    error=function(e) rasterImage(dat[[img.l[j]]][top:bottom,left:right],xleft,ybottom,xright,ytop)
                )
            }
        }
    }

        tryCatch(
			legend(x=par("usr")[2]-xinch(1.2), y=par("usr")[3]-yinch(1.6), xpd=TRUE, inset=c(0,-.14), bty="n", cex=.7, legend=dat.name),
		error=function(e) NULL)


        #if(!is.null(pdf.name))
        #{dev.off()}

    #return(pic.pos)
}

#' Function to Clean up the transcriptome info data.frame
#' 
#' @export
ts_info_cleaner <- function(ts_info){
    # Get rid of any data that has not been sequenced
    runLogic <- !is.na(ts_info$Gnomex.Label)
    ts_RunOnly <- ts_info[runLogic,]

    # Only have data that has an RD file
    rdLogic <- !is.na(ts_RunOnly$rd.name)
    tsInfoClean <- ts_RunOnly[rdLogic, ]

    #remove any DTAM experiments
    tsInfoClean <- tsInfoClean[tsInfoClean$label_experiment != 'DTAM',]

    # Give it some row.names
    row.names(tsInfoClean) <- tsInfoClean$Gnomex.Label

    return(tsInfoClean)
}


#' Function to return Gnomex.labels based on Selections from
#' first the label_experiment Collumn and then the Cell type Collumn
#' @param expSel will allow you to select cells by experiment
#' @param ctSel will allow you to select cells by cell_type
dataSelector <- function(expSel = T, ctSel = T, cellsSelection = NULL){
    cat('Function to Return Gnomex Labels Based on\n1.the experiment type you select\n2.the cell_type you select\nSelecting Cancel will return all values.\n')
    # How I want to First Select Data is based on label_experiment
    experimentTypes <- as.character(unique(tsSuper$ts_info$label_experiment))
    if(expSel){
        cat('\nSelect the label_experiment Types to view your data\n')
        expsToView <- select.list(experimentTypes, multiple = T, title = "label_experiment Types")
    }else{
        expsToView <- vector(mode="numeric", length=0)
    }

    if( length(expsToView) == 0 ){
        tsSuperReduce1 <- tsSuper$ts_info
    }else{
        expsToViewLogic <- tsSuper$ts_info$label_experiment %in% expsToView
        tsSuperReduce1 <- tsSuper$ts_info[expsToViewLogic,, drop=F]
    }
    
    # cat('\nBelow You can See a summary of all the Cell Types Available To you\nfrom you initial Selection\n')
    # print(sort(table(as.character((tsSuperReduce1$label_cellType)))))
    if(ctSel){
        cell_types <- as.character(unique(tsSuperReduce1$label_cellType))
        cell_typesToView <- select.list(sort(cell_types), multiple = T, title = 'Cell Types', preselect=cellsSelection)
    }else{
        cell_typesToView <- vector(mode="numeric", length=0)
    }

    if( length(cell_typesToView) ==0){
        tsSuperReduce2 <- tsSuperReduce1 
    }else{
        tsSuperReduce2Logic <- tsSuperReduce1$label_cellType %in% cell_typesToView
        tsSuperReduce2 <- tsSuperReduce1[tsSuperReduce2Logic,,drop=F]
    }

    gnomexLabs <- as.character(tsSuperReduce2$Gnomex.Label)
    ts_info_reduce <- tsSuper$ts_info[gnomexLabs, , drop=F]
    ts_info_reduce <- ts_info_reduce[order(ts_info_reduce$label_cellType_numeric),,drop=F ]
    # newOrder <- as.character(ts_info_reduce$Gnomex.Label)
    # tsSuper$ts_data <<- tsSuper$ts_data[newOrder,]


    output <- list()
    output[['ts_info_reduce']] <- ts_info_reduce
    output[['preSelect']] <- as.character(ts_info_reduce$label_cellType)
    return(output)

}

#' Function to return genes based on 
Gene.Go.Finder <- function(terms="neuron"){
        i <- terms[1]
        g.names <- unique(tsSuper[['gene.go']][grep(i,tsSuper[['gene.go']][,"GO.term.name"],ignore.case=T),1])
        for(i in terms[-1])
        {
                ihit <- unique(tsSuper[['gene.go']][grep(i,tsSuper[['gene.go']][,"GO.term.name"],ignore.case=T),1])
                g.names <- intersect(g.names,ihit)
        }
        g.names <- as.character(tsSuper[['gene.go']][match(g.names,tsSuper[['gene.go']][,1]),"Gene.name"])
        return(g.names)
}

#' Smart Function to find genes either by gene names or 
#' go terms
geneFinder <- function(toSearch){
    genes <- colnames(tsSuper$ts_data)
    foundGenes <- c()
    for(i in 1:length(toSearch)){
        found <- sort(grep(toSearch[i], genes, value = T, ignore.case = T), decreasing=T)
        foundGenes <- c(found, foundGenes)
    }
    foundGenes <- rev(foundGenes)

    if(length(foundGenes) == 0){
        foundGenes <- sort(Gene.Go.Finder(c('voltage', 'sodium')))
    }

    if(length(foundGenes) == 0){
        stop("Function found no genes.")
    }else{
        return(foundGenes)
    }
}

#' Function to display a dataframe in HTML
heatMapper <- function(geneDF){
    geneDF <- cbind(as.character(row.names(geneDF)), as.data.frame(geneDF))
    colnames(geneDF)[1] <- "Gene Names"
    assign('cf', condformat::condformat(geneDF), envir = .GlobalEnv)

    for(i in 2:dim(geneDF)[2]){
        assign(
            'cf', 
            condformat::rule_fill_gradient2(
                get('cf', .GlobalEnv), 
                !!i, 
                low='white', 
                high='green')
        ,envir = .GlobalEnv)
    }
    condformat::condformat2latex(get('cf', .GlobalEnv))
}

#' Function to select the search term file
#' returns a vector of search terms.
#' Should directly be ported into geneFinder
searchSelector <- function(){
    allFiles <- list.files(recursive=T)
    searchFiles <- grep(".*searches.*", allFiles, value=T)

    filesToSelect <- rapply(strsplit(searchFiles, '[/]'), function(x) rev(x)[[1]])
    fileSelection <- menu(filesToSelect)
    bringToTop(-1)
    fileSelectionPath <- paste0('./',searchFiles[fileSelection])

    toSearch <- scan(fileSelectionPath, 'character', sep='\n', quiet=T)
    return(toSearch)
}

#This peak func allows for multiple t.types to be plotted
#170515: added pts and lns: (logical)
#added dat.n for insertation of the name for the rd file
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

matrixWrangler <- function(geneDF, scale = c('none')){
    # Scale the dataframe the way you define    
    if('log' %in% scale){
        geneDF <- log(geneDF + 1)
    }else if('row' %in% scale){
        geneDF <- scale(geneDF)
    }else if ('column' %in% scale) {
       geneDF <- scale(t(geneDF))
       geneDF <- t(geneDF)
    }
    
    return(geneDF)
}
# hi<- heatMapper(geneDF)

# hi<- knit_print(hi)

# xtable(hi)

