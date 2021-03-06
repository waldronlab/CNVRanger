################## 
# Author: Vinicius Henrique da Silva
# Description: Graphical representations for CNV-GWAS
# Date: April 11, 2018   
# Code: GRAPCNVGWAS001
##################

#' Violin plot to show the quantitative phenotype distribution for each CNV genotype 
#' 
#' @param all.paths phen.info$all.paths  
#' @param snp.name SNP name. "Name" returned by \code{\link{cnvGWAS}}  
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples 
#' # Write the pdf plot in disk
#'  all.paths <- phen.info$all.paths
#'  snp.name <- "AX-76818525"
#'  violinAllele(all.paths, snp.name)
#'  
#' @export

violinAllele <- function(all.paths, snp.name, lo.phe=1, ylab=NULL, size.font=15,
                         legend=TRUE, size.axis.font = 15, title=NULL, plot.pdf=TRUE){
  
  ### Open gds
  cnv.gds <- file.path(all.paths[1], "CNV.gds")
  genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork = TRUE, readonly = FALSE)
  
  ### Get the allele
  snp.ids <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
  snp.number <- which(snp.ids==snp.name)
  
  ## Get the genotypes
  col.gen <- length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
  g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(snp.number,1), count=c(1, col.gen))
  
  phen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype"))
  SNPRelate::snpgdsClose(genofile)
  
  phen$gen <- g
  
  phen[,lo.phe+1] <- gsub(-9, NA, phen[,lo.phe+1])
  
  vio.df <- phen[,c(lo.phe+1, ncol(phen))]
  phen.name <- colnames(phen)[lo.phe+1]
  colnames(vio.df) <- c(phen.name, "state")
  vio.df$state <- as.vector(as.character(vio.df$state)) 
  vio.df[, which(colnames(vio.df)==phen.name)] <- as.numeric(vio.df[, which(colnames(vio.df)==phen.name)])
  
  ## Avoid nonexistent states 
  vio.df$state <- as.character(vio.df$state)
  ## Delete NA´s - why are they there?
  vio.df <- na.omit(vio.df)
  ## Call n 
  vio.df$state <- paste0(vio.df$state, "n")
  
  #print(vio.df$state)
  #print(class(vio.df$state))
  #print(table(vio.df$state))
  #print(head(vio.df))
  
if(is.null(ylab)){
  ylab <- colnames(vio.df)[1]}
  
  ## Define colors of states
  states.in.plot <- as.character(names(table(vio.df$state)))
  states.in.plot[which(states.in.plot=="0n")] <- "dodgerblue3"
  states.in.plot[which(states.in.plot=="1n")] <- "darkolivegreen3"
  states.in.plot[which(states.in.plot=="2n")] <- "white"
  states.in.plot[which(states.in.plot=="3n")] <- "coral"
  #print(states.in.plot)
  
  ## Plot the violin 
  p <- ggplot2::ggplot(vio.df, ggplot2::aes_string(x="state", y=phen.name, fill="state")) + 
    ggplot2::geom_violin(trim=FALSE)
  p <- p + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), color = "gray48") +
    ggplot2::stat_summary(fun=median, geom="point", size=2.5, color="black") +
    ggplot2:: stat_summary(fun=median, colour="grey", geom="line", ggplot2::aes(group = 1)) +
    ggplot2::xlab("") + ggplot2::ylab(ylab) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0),
                                                        size=size.font) ,
                   axis.text.y =  ggplot2::element_text(size = size.axis.font),
                   axis.text.x =  ggplot2::element_text(size = size.axis.font)) +
    ggplot2::theme(legend.title=ggplot2::element_blank())

   p <- p + ggplot2::scale_fill_manual(values=states.in.plot)
   #text = element_text(size=20)
   
   #axis.text.y = element_text(face="bold", color="#993333", 
  #                            size=14, angle=45))
   
   if(!legend){
   p <- p + ggplot2::theme(legend.position = "none")
   }
   
   if(!is.null(title)){
    p <- p + ggplot2::ggtitle(title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
   
     
  if(plot.pdf){
    pdf.file <- file.path(all.paths[3], paste0(phen.name, "-violin-", snp.name,".pdf"))
    pdf(pdf.file)
    print(p)
    invisible(dev.off())
  }else{
    return(p)
  }
  
}

if(FALSE){
  all.paths <- phen.info$all.paths
  snp.name <- "AX-76597039"
  CNVRanger::violinAllele(all.paths, snp.name)
}

#' Manhattan Plot
#' 
#' Manhattan plot for p-values of a CNV-GWAS
#'
#' @param all.paths TODO. 
#' @param regions \code{\linkS4class{GRanges}} as returned by \code{\link{cnvGWAS}}  
#' @param chr.size.order \code{\link{data.frame}} with two columns: (i) 'chr': chromosome names (character),
#' and (ii) 'size': length of the chromosomes in bp (integer)
#' @param plot.pdf Logical plot a to pdf file
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples 
#' 
#' # Load phenotype-CNV information
#' data.dir <- system.file("extdata", package="cnvAnalyzeR")
#' 
#' phen.loc <- file.path(data.dir, "Pheno.txt")
#' cnv.out.loc <- file.path(data.dir, "CNVOut.txt")
#' map.loc <- file.path(data.dir, "MapPenn.txt")
#'
#' phen.info <- setupCnvGWAS('Example', phen.loc, cnv.out.loc, map.loc)
#' all.paths <- phen.info$all.paths
#' segs.pvalue.gr <- cnvGWAS(phen.info)
#'
#' # Define the chromosome order in the plot
#' order.chrs <- c(1:24, "25LG1", "25LG2", 26:28, "LGE22", "1A", "4A")
#'
#' # Chromosome sizes
#' chr.sizes <- c(  114059860, 150265477, 111636321, 68030631, 61875929, 34342011,
#'                  37653027, 31324166, 25106277, 20202851, 20315886, 20466350,
#'                  16480340, 16193477, 13820886, 10486032, 11572643, 9871655,
#'                  14661763, 7693166, 4276343, 6655392, 6808513, 1092960, 809223,
#'                  6596997, 4327975, 5101010, 773534, 71365269, rep(19934923,2))
#'
#' chr.size.order <- data.frame(chr=order.chrs, sizes=chr.sizes, stringsAsFactors=FALSE)
#'
#' # Plot Manhatthan to a pdf within the 'Results' workfolder
#' plotManhattan(all.paths, segs.pvalue.gr, chr.size.order)
#'
#' @export
#'  

plotManhattan <- function(all.paths, regions, chr.size.order = NULL, 
                          plot.pdf = TRUE) {
    
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
    
    
    gwasResults <- data.frame(CHR = as.character(GenomicRanges::seqnames(regions)), 
        BP = GenomicRanges::start(regions), SNP = regions$SegName, P = regions$MinPvalueAdjusted, 
        stringsAsFactors = FALSE)
    
    gwasResults <- rbind(gwasResults, chr.all)
    
    ### turn chr to numeric
    for (cn in seq_len(nrow(chr.size.order))) {
        gwasResults$CHR <- gsub(paste0("\\<", chr.size.order[cn, 1], "\\>"), chr.size.order[cn, 
            3], gwasResults$CHR)
    }
    
    gwasResults$CHR <- as.integer(gwasResults$CHR)
    gwasResults$BP <- as.integer(gwasResults$BP)
    gwasResults$P <- as.numeric(gwasResults$P)
    ind <- do.call(order, gwasResults[, c("CHR", "BP")])
    gwasResults <- gwasResults[ind, ]
    
    if (plot.pdf) {
        pdf(file.path(all.paths[3], paste0(phen.name, "-Manhattan.pdf")))
        qqman::manhattan(gwasResults, chr = "CHR", bp = "BP", snp = "SNP", p = "P", 
            chrlabs = chr.size.order$chr, suggestiveline = -log10(0.1), genomewideline = -log10(0.05))
        dev.off()
    } else {
        qqman::manhattan(gwasResults, chr = "CHR", bp = "BP", snp = "SNP", p = "P", 
            chrlabs = chr.size.order$chr, suggestiveline = -log10(0.1), genomewideline = -log10(0.05))
        
    }
}


#' Plot the Manhattan
#' 
#' Function to write a Manhattan plot with the p-values from a CNV-GWAS run or Vst values. 
#' 
#' @param regions GRanges from CNV-GWAS run or Vst values to produce the manhattan plot
#' @param chr.size.order Data-frame with two columns: (i) 'chr' have the chromosome names as character
#' and (ii) 'size' with the length of the chromosome in basepairs
#' @param grob Produce a grob object
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples 
#' 
#' folder.package <- system.file('extdata', package = 'CNVRanger')
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

plotManhattanColor <- function(all.paths, regions, chr.size.order=NULL, 
                              highlight=NULL, write.pdf=FALSE, ylim.man=NULL, cex.lab=1.8, 
                              cex.axis=1.3, keep.y.lab=TRUE, keep.x.lab=TRUE,
                              result.col.name="MinPvalueAdjusted", result.metric="q"){
  
  grob <- !write.pdf
  
  min.qvalue <- min(regions$MinPvalueAdjusted)
  minqvalue.log <- max(-log10(min.qvalue))
  if(is.null(ylim.man)){
  ylim.man = minqvalue.log*1.10
  }
  
  ## Check if overlapped genes exist
  if(is.null(values(regions)$Overlapped.genes)){
    values(regions)$Overlapped.genes <- rep(NA, length(regions))
  }
  
  ## Produce chromosome limits
  colnames(chr.size.order) <- c("chrs", "sizes")
  chr.size.order.start <- as.data.frame(chr.size.order$chrs)
  chr.size.order.start$sizes <- 1
  colnames(chr.size.order.start) <- colnames(chr.size.order)
  chr.all <- rbind(chr.size.order, chr.size.order.start)
  chr.all$SNP <- paste0("lim", seq_len(nrow(chr.all)))
  colnames(chr.all)[1:2] <- c("CHR", "BP")
 
  ## Produce numeric code for chrs
  chr.size.order$chr.numeric <- seq_len(nrow(chr.size.order))
  
  ## Simulate dots to color
  all.p.all <- list()
  for(chrx in seq_len(nrow(chr.size.order))){
    lseq <- function(from=1, to=400000, length.out=6) {
      exp(seq(log10(from), log10(to), length.out = length.out))
    }
    
    pointschrX.x <- seq(from=500000, to=chr.size.order$size[chrx], by=5000000) ## simulated positions
    pointschrX.y <- lseq(from=0.000000000000000000001, to=1, 1000) ## simulated p-values
    all.p <- expand.grid(pointschrX.y, pointschrX.x)
    all.p$CHR <- chr.size.order$chr.numeric[chrx]
    all.p.all[[chrx]] <- as.data.frame(all.p)
  }
  
  backgroundDots <- data.table::rbindlist(all.p.all)
  backgroundDots$SNP <- seq_len(nrow(backgroundDots))
  colnames(backgroundDots)[1:2]  <- c("P", "BP")
  backgroundDots <- subset(backgroundDots, select = c(CHR, BP, SNP, P))
  backgroundDots <- as.data.frame(backgroundDots)
  backgroundDots$Overlapped.genes <- NA #blank column
  backgroundDots$Phen <- NA #blank column

    chr.all$P <- 1 
    phen.name <- values(regions)$Phenotype
    phen.name <- unique(phen.name)
    
    p.defined <- eval(parse(text=paste0("values(regions)$", result.col.name)))
    
    if(grob==FALSE){
    pdf(file.path(all.paths[3], paste0(phen.name,"-Manhattan.pdf")))}
    gwasResults <- cbind("CHR"=as.character(seqnames(regions)), "BP"=start(regions),
                         "SNP"= values(regions)$SegName, 
                         "P"= p.defined,
                         "Overlapped.genes"= values(regions)$Overlapped.genes,
                         "Phen"=values(regions)$Phenotype)
    
    gwasResults <- as.data.frame(gwasResults, stringsAsFactors = FALSE)
    chr.all <- subset(chr.all, select = c(CHR, BP, SNP, P))
    chr.all$Overlapped.genes <- NA #blank column
    chr.all$Phen <- NA #blank column
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
    backgroundDots <- rbind(gwasResults, backgroundDots)
    backgroundDots <- dplyr::arrange(backgroundDots, CHR, BP)
    
    qqman.mod(backgroundDots, col = scales::alpha(c("coral", "azure3"), 0.006),
                     chr="CHR", bp="BP", snp="SNP", p="P",
                     chrlabs=chr.size.order$chr, suggestiveline =-log10(0.10),
                     genomewideline = -log10(0.05), ylim = c(0, ylim.man), 
              cex.axis=cex.axis,
              cex.lab=cex.lab, keep.x.lab=keep.x.lab, keep.y.lab=keep.y.lab,
              result.metric=result.metric);par(new=TRUE)
    
    ### Show only highligh segments
    highlight <- gwasResults[which(gwasResults$SNP %in% highlight),]$Overlapped.genes
    highlight <- highlight[highlight != ""] 
    gwasResults$SNP <- gwasResults$Overlapped.genes
    
    qqman.mod(gwasResults, col = c("blue4"), chr="CHR", bp="BP", snp="SNP", p="P",
                     chrlabs=chr.size.order$chr, suggestiveline =-log10(0.10),
                     genomewideline = -log10(0.05), annotatePval = 0.1, ylim = c(0, ylim.man),
    highlight=highlight, cex.axis=cex.axis,
    cex.lab=cex.lab, keep.x.lab=keep.x.lab, keep.y.lab=keep.y.lab,result.metric=result.metric)
    
    pl <- recordPlot()
    
    if(grob){
      pdf(file.path(all.paths[3], paste0(phen.name,"-Manhattan.pdf")))
      pl
    }
    
    dev.off()
    
  if(grob){
    return(pl)}
}


#' Plot QQ-plot using non-log p-values Originally published in:
#' https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R By
#' Matthew Flickinger 
#' @import lattice
#' @export

qqunifPlot <- function(pvalues, should.thin = TRUE, thin.obs.places = 2, thin.exp.places = 2, 
    xlab = expression(paste("Expected (", -log[10], " p-value)")), 
    ylab = expression(paste("Observed (", -log[10], " p-value)")), 
    draw.conf = TRUE, conf.points = 1000, conf.col = "lightgray", 
    conf.alpha = 0.05, already.transformed = FALSE, pch = 20, aspect = "iso", 
    prepanel = prepanel.qqunif, 
    par.settings = list(superpose.symbol = list(pch = pch)), ...) 
{
    # error checking
    if (length(pvalues) == 0) 
        stop("pvalue vector is empty, can't draw plot")
    if (!(is.numeric(pvalues) || (is.list(pvalues) && all(vapply(pvalues, is.numeric, logical(1)))))) 
        stop("pvalue vector is not numeric, can't draw plot")
    if (any(is.na(unlist(pvalues)))) 
        stop("pvalue vector contains NA values, can't draw plot")
    if (already.transformed == FALSE) {
        if (any(unlist(pvalues) == 0)) 
            stop("pvalue vector contains zeros, can't draw plot")
    } else {
        if (any(unlist(pvalues) < 0)) 
            stop("-log10 pvalue vector contains negative values, can't draw plot")
    }
    
    grp <- NULL
    n <- 1
    exp.x <- c()
    if (is.list(pvalues)) {
        nn <- lengths(pvalues)
        rs <- cumsum(nn)
        re <- rs - nn + 1
        n <- min(nn)
        if (!is.null(names(pvalues))) {
            grp = factor(rep(names(pvalues), nn), levels = names(pvalues))
            names(pvalues) <- NULL
        } else {
            grp = factor(rep(seq_along(pvalues), nn))
        }
        pvo <- pvalues
        pvalues <- numeric(sum(nn))
        exp.x <- numeric(sum(nn))
        for (i in seq_along(pvo)) {
            if (!already.transformed) {
                pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
                exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method = "first") - 
                  0.5)/nn[i])
            } else {
                pvalues[rs[i]:re[i]] <- pvo[[i]]
                exp.x[rs[i]:re[i]] <- -log10((nn[i] + 1 - rank(pvo[[i]], ties.method = "first") - 
                  0.5)/(nn[i] + 1))
            }
        }
    } else {
        n <- length(pvalues) + 1
        if (!already.transformed) {
            exp.x <- -log10((rank(pvalues, ties.method = "first") - 0.5)/n)
            pvalues <- -log10(pvalues)
        } else {
            exp.x <- -log10((n - rank(pvalues, ties.method = "first") - 0.5)/n)
        }
    }
    
    
    # this is a helper function to draw the confidence interval
    panel.qqconf <- function(n, conf.points = 1000, conf.col = "gray", conf.alpha = 0.05, 
        ...) {
        conf.points = min(conf.points, n - 1)
        mpts <- matrix(nrow = conf.points * 2, ncol = 2)
        for (i in seq(from = 1, to = conf.points)) {
            mpts[i, 1] <- -log10((i - 0.5)/n)
            mpts[i, 2] <- -log10(qbeta(1 - conf.alpha/2, i, n - i))
            mpts[conf.points * 2 + 1 - i, 1] <- -log10((i - 0.5)/n)
            mpts[conf.points * 2 + 1 - i, 2] <- -log10(qbeta(conf.alpha/2, i, n - 
                i))
        }
        grid::grid.polygon(x = mpts[, 1], y = mpts[, 2], gp = grid::gpar(fill = conf.col, 
            lty = 0), default.units = "native")
    }
    
    # reduce number of points to plot
    if (should.thin) {
        if (!is.null(grp)) {
            thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places), 
                exp.x = round(exp.x, thin.exp.places), grp = grp))
            grp = thin$grp
        } else {
            thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places), 
                exp.x = round(exp.x, thin.exp.places)))
        }
        pvalues <- thin$pvalues
        exp.x <- thin$exp.x
    }
    gc()
    
    prepanel.qqunif = function(x, y, ...) {
        A = list()
        A$xlim = range(x, y) * 1.02
        A$xlim[1] = 0
        A$ylim = A$xlim
        return(A)
    }
    
    # draw the plot
    lattice::xyplot(pvalues ~ exp.x, groups = grp, xlab = xlab, ylab = ylab, aspect = aspect, 
        prepanel = prepanel, scales = list(axs = "i"), pch = pch, panel = function(x, 
            y, ...) {
            if (draw.conf) {
                panel.qqconf(n, conf.points = conf.points, conf.col = conf.col, conf.alpha = conf.alpha)
            }
            lattice::panel.xyplot(x, y, ...)
            lattice::panel.abline(0, 1)
        }, par.settings = par.settings, ...)
}


## HELPER - Function adapted from qqman package
qqman.mod <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
          "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
          genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
          annotatePval = NULL, annotateTop = TRUE, color.highlight="red", cex.axis,
          cex.lab, keep.x.lab, keep.y.lab, result.metric, ...) 
{
  
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }else{
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    if(keep.x.lab){
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    }else{
    xlabel = paste("", unique(d$CHR), "position")
    }
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                                                             i, ]$pos))/2 + 1)
    }
    
    if(keep.x.lab){
    xlabel = "Chromosome"
    }else{
      xlabel = ""
    }
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  if(keep.y.lab){
  ylab <- expression(-log[10](italic(result.metricW)))
  ylab <- parse(text = gsub("result.metricW", result.metric, ylab, fixed = TRUE))
  }else{
    ylab <- ""  
  }
  
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), 
                   cex.axis=cex.axis, cex.lab=cex.lab,
                   #cex=1.2,  
                   #cex.main=1.2,
                   ylim = c(0, 
                  #ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
                  ceiling(max(d$logp))), xlab = xlabel, ylab = ylab,
                  main=unique(x$Phen))
  dotargs <- list(...)
  par(mar = c(5, 5, 5, 5)) # Set the margin on all sides to 6
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, cex=cex.axis, cex.axis=cex.axis, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, cex=cex.axis, cex.axis=cex.axis, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, 
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    #with(d.highlight, points(pos, logp, col = color.highlight, pch = 20, 
    with(d.highlight, points(pos, logp, col = color.highlight, pch = 19, 
                             ...))
  }
  if (!is.null(annotatePval)) {
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), calibrate::textxy(pos, -log10(P), 
                                              offset = 0.625, labs = topHits$SNP, 
                                              #offset = 0.625, labs = topHits$Overlapped.genes, 
                                                cex = 0.8), 
           ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      calibrate::textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.5, 
             #labs = topSNPs$SNP, cex = 0.8, ...)
             #labs = topHits$Overlapped.genes, cex = 1.3, ...)
             labs = topSNPs$SNP, cex = 0.9, font=2, ...)
    }
  }
  par(xpd = FALSE)
}
