##################
# Author: Vinicius Henrique da Silva
# Script description: Graphical representations CNV-GWAS package 
# Date: April 11, 2018 
# Code: GRAPCNVGWAS001
###################

#' HELPER - The function to plot the violin plot. Copied from easyGgplot2 github package https://github.com/kassambara/easyGgplot2
#' 
#' @param data data.frame or a numeric vector. Columns are variables and rows
#' are observations.
#' @param xName The name of column containing x variable (i.e groups). Default
#' value is NULL.
#' @param yName The name of column containing y variable. If yName=NULL, data
#' should be a numeric vector.
#' @param groupName The name of column containing group variable. This variable
#' is used to color plot according to the group.
#' @param addMean if TRUE, the mean point is added on the plot for each group.
#' Default value is FALSE.
#' @param meanPointShape The shape of mean point.
#' @param meanPointSize The size of mean point
#' @param meanPointColor Border color of the mean point. Default value is
#' "black".
#' @param meanPointFill Fill color of mean point. This parameter is used only
#' when meanPointShape=21 to 25. Default value is "blue"
#' @param addDot If TRUE, dotplot is added on the boxplot. Default value is
#' FALSE.
#' @param dotSize The size of dots.
#' @param dotPosition Possible values are "center" and "jitter". Default value
#' is "center".
#' @param jitter Degree of jitter in x direction. Default value is 0.2.
#' @param groupColors Color of groups. groupColors should have the same length
#' as groups.
#' @param brewerPalette This can be also used to indicate group colors. In this
#' case the parameter groupColors should be NULL. e.g: brewerPalette="Paired".
#' @param \dots Other arguments passed on to ggplot2.customize custom function
#' or to geom_dotplot and to geom_violin functions from ggplot2 package.
#' @return a ggplot
#' @author Alboukadel Kassambara <alboukadel.kassambara@@gmail.com>
#' @seealso \code{\link{ggplot2.dotplot}, \link{ggplot2.boxplot},
#' \link{ggplot2.stripchart}, \link{ggplot2.density}, \link{ggplot2.histogram}}
#' @references http://www.sthda.com
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' df <- ToothGrowth
#' ggplot2.violinplot(data=df, xName='dose',yName='len',
#'                 mainTitle="Plot of length according\n to the dose",
#'                 xtitle="Dose (mg)", ytitle="Length")
#' 
#' #Or use this
#' plot<-ggplot2.violinplot(data=df, xName='dose',yName='len')
#' plot<-ggplot2.customize(plot, mainTitle="Plot of length according\n to the dose",
#'                         xtitle="Dose (mg)", ytitle="Length")
#' print(plot)
#' 
#' @export ggplot2.violinplot

ggplot2.violinplot<-function(data, xName=NULL, yName=NULL, groupName=NULL,
                             addMean=F, meanPointShape=5, meanPointSize=4, meanPointColor="black", meanPointFill="blue",
                             addDot=FALSE, dotSize=1, dotPosition=c("center", "jitter"), jitter=0.2,
                             groupColors=NULL, brewerPalette=NULL,...)
{
  
  params <- list(...)
  trim <- ifelse(!is.null(params$trim), params$trim, TRUE)
  position <- ifelse(!is.null(params$position), params$position, "dodge")
  
  #if yName is missing or null, data should be a numeric vector
  if(is.null(yName) & !is.numeric(data))
    stop("yName is missing or NULL. In this case data should be a numeric vector")
  #data is a numeric vector
  else if(is.numeric(data)){
    data=cbind(y=data, x=rep(1, length(data)))
    xName="x"
    yName="y"
  }
  #xName is missing or  NULL => single boxplot corresponding to the deistribution of the variable
  #bind group column to data
  if(is.null(xName)){
    data=cbind(data, x=rep(1, nrow(data)))
    xName="x"
  }
  
  #data
  data=data.frame(data)
  data[,xName]=factor(data[,xName])
  #plot
  if(is.null(groupName))p<-ggplot(data=data, aes_string(x=xName, y=yName))
  else{
    p<-ggplot(data=data, aes_string(x=xName, y=yName, fill=groupName))
    data[,groupName]=factor(data[,groupName])#transform groupName to factor
  }
  p<-p+geom_violin(trim = trim, position = position)
  #add Mean point
  if(addMean) p<-p+stat_summary(fun.y=mean, geom='point', shape=meanPointShape, 
                                size=meanPointSize, colour=meanPointColor, fill=meanPointFill)
  #add dot
  if(addDot){
    pms <- .dotplot_params(...)
    if(dotPosition[1]=="center") p<-p+geom_dotplot(binaxis='y', stackdir='center', dotsize=dotSize,
                                                   width = pms$width, stackratio = pms$stackratio)
    else p<-p+geom_jitter(position=position_jitter(jitter), cex=dotSize, shape=16)
  }
  #group colors
  if(!is.null(groupColors)){
    p<-p+scale_fill_manual(values=groupColors)
    p<-p+scale_colour_manual(values=groupColors)
  } 
  else if(!is.null(brewerPalette)){
    p<-p+scale_fill_brewer(palette=brewerPalette)
    p<-p+scale_colour_brewer(palette=brewerPalette, guide="none")
  } 
  #ggplot2.customize : titles, colors, background, legend, ....
  p<-ggplot2.customize(p,...)
  
  p
}


#' HELPER - Plot QQ-plot using non-log p-values 
#' Originally published in: https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R By Matthew Flickinger
#' @import lattice

qqunif.plot<-function(pvalues, 
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
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
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
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}

#' HELPER - Plot the distribuition of the response variable among CNV genotypes
#' This function is limited from 0 to 4 copies.
#' @example 
#' .plotFreqViolin(all.paths, segs.pvalue.gr, rank.p=2)

.plotFreqViolin  <- function(all.paths, segs.pvalue.gr, rank.p=1, phen.name=NULL, log.transform=FALSE){
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
  
  probeX <- values(segs.pvalue.gr[rank.p])$NameProbe
  map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")), 
                             chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                               position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")), 
                               probes = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
  ))
  
  snp.id <- subset(map, probes==probeX)$snp.id
  State <- .snpgdsGetGenoCNV(genofile,snp.id)
  
  all.samples <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
  phenotypes <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype"))
  
  if(!identical(phenotypes$samplesPhen, all.samples)){
    stop("Something went wrong with the list of phenotypes. Different order of samples")
  }
  
  sam.cnv <- cbind(State,phenotypes)
  sam.cnv <- as.data.frame(sam.cnv)
  
  if(log.transform==TRUE){
    sam.cnv[,3] <- log2(sam.cnv[,3])
  }
  
  SNPRelate::snpgdsClose(genofile)
  
  ### Define phenotype name 
  if(is.null(phen.name)){
    phen.name <- names(sam.cnv)[3]
  }
  
  ### Define colors to the figure
  
  color <- c("steelblue", "darkseagreen3", "white", "coral", "lavenderblush1")
  all.states <- names(table(sam.cnv$State))
  all.states <- as.numeric(all.states)+1
  color <- color[all.states]

  
  #### Export figure
  pdf(file.path(all.paths[3], paste0(phen.name,"ViolinCNVGenotypes.pdf")))
  
  plot <- ggplot2.violinplot(data=sam.cnv, xName="State", yName=phen.name, addDot=TRUE, dotSize=1.7,
                             dotPosition="jitter", jitter=0.2, groupColors = color, groupName= 'State', showLegend=FALSE)
  
  plot <- plot + stat_summary(fun.data=mean_sdl,
                              geom="pointrange", color="grey", size = 1)
  
  stat_sum_single <- function(fun, geom="point", ...) {
    stat_summary(fun.y=fun, colour="gray48", geom=geom, size = 1, ...)
  }
  
  plot + stat_sum_single(mean, geom="line", aes(group=phen.name)) + theme_classic() + theme(text = element_text(size=10))
  
  print(plot)
  dev.off()
  
  message(paste0("mean 2n = ", mean(sam.cnv[which(sam.cnv$State == 2),]$BT)))
  message(paste0("mean 3n = ", mean(sam.cnv[which(sam.cnv$State == 3),]$BT)))
  
}


#' HELPER - Plot the manhattan
#' @param regions Genomic ranges to produce the manhattan plot
#' @param type.regions 'p-value' or/and 'vst'
#' @param chr.size.order Data-frame with two columns: (i) 'chr' have the chromosome names as character
#' and (ii) 'size' with the length of the chromosome in basepairs

.plotMan <- function(all.paths, regions, type.regions, chr.size.order=NULL){
  
  ## Produce chromosome limits
  chr.size.order.start <- chr.size.order
  chr.size.order.start$sizes <- 1
  chr.all <- rbind(chr.size.order, chr.size.order.start)
  chr.all$SNP <- paste0("lim", seq_len(nrow(chr.all)))
  colnames(chr.all)[1:2] <- c("CHR", "BP")
  
  ## Produce numeric code for chrs
  
  chr.size.order$chr.numeric <- seq_len(nrow(chr.size.order))
  
  if("p-value" %in% type.regions){
    chr.all$P <- 1 
    phen.name <- values(regions)$Phenotype
    phen.name <- unique(phen.name)
    
    pdf(file.path(all.paths[3], paste0(phen.name,"-Manhattan.pdf")))
    gwasResults <- cbind("CHR"=as.character(seqnames(regions)), "BP"=start(regions),
                         "SNP"= values(regions)$SegName, "P"= values(regions)$MinPvalueAdjusted)
    
    gwasResults <- as.data.frame(gwasResults)
    gwasResults <- rbind(gwasResults, chr.all)
    
    ### turn chr to numeric
    for(cn in seq_len(nrow(chr.size.order))){
      gwasResults$CHR <- gsub(paste0("\\<",chr.size.order[cn,1], "\\>"),
                              chr.size.order[cn,3], gwasResults$CHR)
    }
    
    gwasResults$CHR <- as.numeric(gwasResults$CHR)
    gwasResults$BP <- as.numeric(gwasResults$BP)
    gwasResults$P <- as.numeric(gwasResults$P)
    gwasResults <- dplyr::arrange(gwasResults, CHR, BP)
    
    qqman::manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P",
                     chrlabs=chr.size.order$chr, suggestiveline =-log10(0.10),
                     genomewideline = -log10(0.05))
    dev.off()
    
  }
  else if("vst" %in% type.regions){
    chr.all$P <- 0 
    
    values(map.vst)$Populations
    
    gwasResults <- cbind("CHR"=as.character(seqnames(regions)), 
                         "BP"=start(regions),
                         "SNP"= values(regions)$NameProbe, 
                         "P"= values(regions)$vst.result)
    
    gwasResults <- as.data.frame(gwasResults)
    gwasResults <- rbind(gwasResults, chr.all)
    
    ### turn chr to numeric
    for(cn in seq_len(nrow(chr.size.order))){
      gwasResults$CHR <- gsub(paste0("\\<",chr.size.order[cn,1], "\\>"),
                              chr.size.order[cn,3], gwasResults$CHR)
    }
    
    gwasResults$CHR <- as.numeric(gwasResults$CHR)
    gwasResults$BP <- as.numeric(gwasResults$BP)
    gwasResults$P <- as.numeric(gwasResults$P)
    gwasResults <- dplyr::arrange(gwasResults, CHR, BP)
    
    
    pdf(file.path(all.paths[3], paste0("VSTManhattan.pdf")))  
    
    qqman::manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, 
                     ylab = "Vst")
    
    dev.off()
    
    
  }else{
    stop("Unknown type.regions argument")
  }
  
}


#' HELPER Upload all figures to dropbox 
#' 
#' @param path.name Folder name in the Dropbox website
#' @example 
#' 
#' Upload all '.pdf' figures 
#' .uploadAllFigures(all.paths, "CNV_BreedingTime")

.uploadAllFigures <- function(all.paths, drop.path, format=".pdf"){

if(length(list.files(all.paths[1])[list.files(all.paths[1])=="token.rds"]) ==0){
  stop("token.rds is not in the Inputs folder")
}
  
x <- readRDS(file=file.path(all.paths[1], "token.rds"))

all.fig <- list.files(all.paths[3], pattern=format)

for(up in seq_len(length(all.fig))){
rdrop2::drop_upload(file.path(all.paths[3], all.fig[up]), path = drop.path,  dtoken = x)}

}


