##################
# Author: Vinicius Henrique da Silva
# Script description: Graphical representations CNV-GWAS package 
# Date: April 11, 2018 
# Code: GRAPCNVGWAS001
###################


# HELPER - Plot QQ-plot using non-log p-values 
# Originally published in: https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R By Matthew Flickinger
# @import lattice

.qqunifPlot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-lengths(pvalues)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(seq_along(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in seq_along(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid::grid.polygon(x=mpts[,1],y=mpts[,2], gp=grid::gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  lattice::xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           lattice::panel.xyplot(x,y, ...);
           lattice::panel.abline(0,1);
         }, par.settings=par.settings, ...
   )
}

#' Plot the manhattan
#' 
#' Function to write a manhattan plot with the p-values from a CNV-GWAS.
#'
#' @param all.paths TODO. 
#' @param regions GRanges from CNV-GWAS run to produce the manhattan plot
#' @param chr.size.order Data-frame with two columns: (i) 'chr' have the chromosome names as character
#' and (ii) 'size' with the length of the chromosome in basepairs
#' @param plot.pdf Produce a pdf file
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples 
#' 
#' folder.package <- system.file('extdata', package = 'cnvAnalyzeR')
#'
#' phen.info <- setupCnvGWAS(name="Example", phen.loc=file.path(folder.package, "Pheno.txt"), ## GEBV values
#'                          cnv.out.loc <- file.path(folder.package, "CNVOut.txt"), ## CNV calls            
#'                          map.loc <- file.path(folder.package, "MapPenn.txt")) ## Genomic positions of the SNPs used in the CNV call
#'
#' all.paths <- phen.info$all.paths
#'
#' segs.pvalue.gr <- cnvGWAS(phen.info)
#' 
#' ## Define the chromosome order in the plot
#' order.chrs <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19",
#'                "20", "21", "22", "23", "24", "25LG1", "25LG2", "26", "27", "28", "LGE22", "1A", "4A")
#'
#' ## Chromosome sizes
#' chr.sizes <- c(114059860,150265477,111636321,68030631,61875929,34342011,37653027, 31324166,25106277,20202851,20315886,20466350,16480340,16193477,
#'               13820886,10486032,11572643,9871655,14661763,7693166,4276343, 6655392,6808513,1092960,809223,6596997,4327975,5101010, 773534,71365269,19934923)
#'
#' chr.size.order <- cbind(as.data.frame(order.chrs, stringsAsFactors = FALSE),as.data.frame(chr.sizes,  stringsAsFactors = FALSE))
#' colnames(chr.size.order) <- c("chr", "sizes")
#'
#' ## Plot a pdf file with a manhatthan within the 'Results' workfolder
#' plotManhattan(all.paths, segs.pvalue.gr, "p-value", chr.size.order)
#'
#' @export
#'  

plotManhattan <- function(all.paths, regions, chr.size.order=NULL, 
                          plot.pdf=TRUE){
  
  ## Produce chromosome limits
  chr.size.order.start <- chr.size.order
  chr.size.order.start$sizes <- 1
  chr.all <- rbind(chr.size.order, chr.size.order.start)
  chr.all$SNP <- paste0("lim", seq_len(nrow(chr.all)))
  colnames(chr.all)[1:2] <- c("CHR", "BP")
  
  ## Produce numeric code for chrs
  
  chr.size.order$chr.numeric <- seq_len(nrow(chr.size.order))
  
    chr.all$P <- 1 
    phen.name <- regions$Phenotype
    phen.name <- unique(phen.name)
    

    gwasResults <- data.frame(  CHR=as.character(GenomicRanges::seqnames(regions)), 
                                BP=GenomicRanges::start(regions),
                                SNP=regions$SegName, 
                                P=regions$MinPvalueAdjusted,
                                stringsAsFactors = FALSE)
    
    gwasResults <- rbind(gwasResults, chr.all)
    
    ### turn chr to numeric
    for(cn in seq_len(nrow(chr.size.order))){
      gwasResults$CHR <- gsub(paste0("\\<",chr.size.order[cn,1], "\\>"),
                              chr.size.order[cn,3], gwasResults$CHR)
    }
    
    gwasResults$CHR <- as.integer(gwasResults$CHR)
    gwasResults$BP <- as.integer(gwasResults$BP)
    gwasResults$P <- as.numeric(gwasResults$P)
    ind <- do.call(order, gwasResults[,c("CHR", "BP")])
    gwasResults <- gwasResults[ind,]
    
    if(plot.pdf){
      pdf(file.path(all.paths[3], paste0(phen.name,"-Manhattan.pdf")))
    qqman::manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P",
                     chrlabs=chr.size.order$chr, suggestiveline =-log10(0.10),
                     genomewideline = -log10(0.05))
    dev.off()
    }else{
      qqman::manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P",
                       chrlabs=chr.size.order$chr, suggestiveline =-log10(0.10),
                       genomewideline = -log10(0.05))
    
    }
    
  
}

