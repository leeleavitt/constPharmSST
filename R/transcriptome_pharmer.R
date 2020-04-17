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
tsInfoCleaner <- function(ts_info){
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
geneGoFinder <- function(terms="neuron"){
    i <- terms[1]
    g.names <- unique(tsSuper[['gene.go']][grep(i,tsSuper[['gene.go']][,"GO.term.name"],ignore.case=T),1])
    for(i in terms[-1])
    {
            ihit <- unique(tsSuper[['gene.go']][grep(i,tsSuper[['gene.go']][,"GO.term.name"],ignore.case=T),1])
            g.names <- union(g.names,ihit)
    }
    g.names <- as.character(tsSuper[['gene.go']][match(g.names,tsSuper[['gene.go']][,1]),"Gene.name"])
    g.names <- sort(g.names, TRUE)
    g.names <- intersect(g.names, colnames(tsSuper$ts_data))
    g.names <- unique(g.names)
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
        for(i in 1:length(toSearch)){
            found <- sort(geneGoFinder(toSearch))
            foundGenes <- c(found, foundGenes)
        }
    }

    if(length(foundGenes) == 0){
        stop("Function found no genes.")
    }else{
        return(foundGenes)
    }
}

#' Function to display a gene data frame in HTML
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
    search <- select.list(filesToSelect,title = 'Select search')
    fileSelectionPath <- paste0('./searches/',search)

    toSearch <- scan(fileSelectionPath, 'character', sep='\n', quiet=T)
    return(toSearch)
}

#' Function to scale the input matrix
matrixWrangler <- function(geneDF, scale = c('none')){
    # Scale the dataframe the way you define    
    if('log' %in% scale){
        geneDF <- log(geneDF + 1)
    }else if('row' %in% scale){
        geneDF <- scale(geneDF, FALSE, scale = colSums(geneDF))
    }else if ('column' %in% scale) {
       geneDF <- scale(t(geneDF), FALSE, scale = colSums(t(geneDF)))
       geneDF <- t(geneDF)
    }
    
    return(geneDF)
}

#' Function to ready the keyboard
readkeygraph <- function(prompt){
    getGraphicsEvent(prompt = prompt, 
                 onMouseDown = NULL, onMouseMove = NULL,
                 onMouseUp = NULL, onKeybd = onKeybd,
                 consolePrompt = "uh")
    Sys.sleep(0.01)
    return(keyPressed)
}

#' Function to read in the keypresses
onKeybd <- function(key){
    keyPressed <<- key
}
# hi<- heatMapper(geneDF)

# hi<- knit_print(hi)

# xtable(hi)

