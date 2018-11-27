################## 
# Author: Vinicius Henrique da Silva
# Description: Graphical representations for CNV-GWAS
# Date: April 11, 2018   
# Code: GRAPCNVGWAS001
##################

#' Manhattan Plot
#' 
#' Manhattan plot for p-values of a CNV-GWAS
#'
#' @param all.paths TODO. 
#' @param regions \code{\linkS4class{GRanges}} as returned by \code{\link{cnvGWAS}}  
#' @param chr.size.order \code{\link{data.frame}} with two columns: (i) 'chr': chromosome names (character),
#' and (ii) 'size': length of the chromosomes in bp (integer)
#' @param plot.pdf Logical plot a to pdf file
#' @author Vinicius Henrique da Silva <vinicius.dasilva@wur.nl>
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
plotManhattan <- function(all.paths, regions, chr.size.order = NULL, plot.pdf = TRUE) {
    
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
#' @param type.regions 'p-value' or/and 'vst'
#' @param chr.size.order Data-frame with two columns: (i) 'chr' have the chromosome names as character
#' and (ii) 'size' with the length of the chromosome in basepairs
#' @param grob Produce a grob object
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

plotManhattanColor <- function(all.paths, regions, type.regions, chr.size.order=NULL, 
                               grob=FALSE){
  
  ## Produce chromosome limits
  chr.size.order.start <- chr.size.order
  chr.size.order.start$sizes <- 1
  chr.all <- rbind(chr.size.order, chr.size.order.start)
  chr.all$SNP <- paste0("lim", seq_len(nrow(chr.all)))
  colnames(chr.all)[1:2] <- c("CHR", "BP")
  
  ## Produce numeric code for chrs
  chr.size.order$chr.numeric <- seq_len(nrow(chr.size.order))
  
  ## Simulate dots to color
  all.p.all <- list()
  for(chrx in seq_len(nrow(chr.size.order))){
    lseq <- function(from=1, to=100000, length.out=6) {
      # logarithmic spaced sequence
      # blatantly stolen from library("emdbook"), because need only this
      exp(seq(log10(from), log10(to), length.out = length.out))
    }
    pointschrX.x <- seq(from=500000, to=chr.size.order$size[chrx], by=5000000) ## simulated positions
    pointschrX.y <- lseq(from=0.0000003, to=1, 1000) ## simulated p-values
    all.p <- expand.grid(pointschrX.y, pointschrX.x)
    all.p$CHR <- chr.size.order$chr.numeric[chrx]
    #all.p$Var2 <- all.p$Var2+seq_len(nrow(all.p))
    all.p.all[[chrx]] <- as.data.frame(all.p)
  }
  
  backgroundDots <- data.table::rbindlist(all.p.all)
  backgroundDots$SNP <- seq_len(nrow(backgroundDots))
  colnames(backgroundDots)[1:2]  <- c("P", "BP")
  backgroundDots <- subset(backgroundDots, select = c(CHR, BP, SNP, P))
  backgroundDots <- as.data.frame(backgroundDots)
  backgroundDots$Overlapped.genes <- NA #blank column
  
  if("p-value" %in% type.regions){
    chr.all$P <- 1 
    phen.name <- values(regions)$Phenotype
    phen.name <- unique(phen.name)
    
    if(grob==FALSE){
      pdf(file.path(all.paths[3], paste0(phen.name,"-Manhattan.pdf")))}
    gwasResults <- cbind("CHR"=as.character(seqnames(regions)), "BP"=start(regions),
                         "SNP"= values(regions)$SegName, 
                         "P"= values(regions)$MinPvalueAdjusted,
                         "Overlapped.genes"= values(regions)$Overlapped.genes)
    
    gwasResults <- as.data.frame(gwasResults, stringsAsFactors = FALSE)
    #chr.all <- chr.all[,-3]
    #chr.all <- chr.all[,colnames(chr.all)%in%"chr.numeric"]
    chr.all <- subset(chr.all, select = c(CHR, BP, SNP, P))
    chr.all$Overlapped.genes <- NA #blank column
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
    
    qqman.mod(backgroundDots, col = scales::alpha(c("black", "grey30"), 0.006),
    #qqman::manhattan(backgroundDots, col = scales::alpha(c("black", "grey30"), 0.006),
                     chr="CHR", bp="BP", snp="SNP", p="P",
                     chrlabs=chr.size.order$chr, suggestiveline =-log10(0.10),
                     genomewideline = -log10(0.05), ylim = c(0, 3));par(new=TRUE)
    
    gwasResults$SNP <- gwasResults$Overlapped.genes
    qqman.mod(gwasResults, col = c("blue4", "orange3"), chr="CHR", bp="BP", snp="SNP", p="P",
    #qqman::manhattan(gwasResults, col = c("blue4", "orange3"), chr="CHR", bp="BP", snp="SNP", p="P",
                     chrlabs=chr.size.order$chr, suggestiveline =-log10(0.10),
                     genomewideline = -log10(0.05), annotatePval = 0.1, ylim = c(0, 3))
    
    pl <- recordPlot()
    
    if(grob){
    #  grab_grob <- function(){
    #    gridGraphics::grid.echo()
    ##    grid::grid.grab()
    #  }
    #  dev.control("enable")
    #  pl <- grab_grob()
      
      pdf(file.path(all.paths[3], paste0(phen.name,"-Manhattan.pdf")))
      pl
    }
    
    dev.off()
    
  }
  else if("vst" %in% type.regions){
    chr.all$P <- 0 
    
    populations <- unique(values(regions)$Populations)
    
    gwasResults <- cbind("CHR"=as.character(seqnames(regions)), 
                         "BP"=start(regions),
                         "SNP"= values(regions)$NameProbe, 
                         "P"= values(regions)$vst.result)
    
    gwasResults <- as.data.frame(gwasResults, stringsAsFactors = FALSE)
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
    
    if(grob == FALSE){
      pdf(file.path(all.paths[3], paste0(populations, "-VSTManhattan.pdf"))) }
    
    #qqman::manhattan(gwasResults,col = scales::alpha(c("blue4", "orange3"), 0.01));par(new=TRUE)
    qqman::manhattan(backgroundDots, col = scales::alpha(c("blue4", "orange3"), 0.006),
                     chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, ylab = "Vst");par(new=TRUE)
    
    qqman::manhattan(gwasResults, col = c("blue4", "orange3"),
                     chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, 
                     ylab = "Vst")
    
    #pl <- recordPlot()
    if(grob){
      grab_grob <- function(){
        gridGraphics::grid.echo()
        grid::grid.grab()
      }
      dev.control("enable")
      pl <- grab_grob()
      
      pdf(file.path(all.paths[3], paste0(populations, "-VSTManhattan.pdf")))  
      pl
    }
    
    dev.off()
    
    
  }else{
    stop("Unknown type.regions argument")
  }
  
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
          annotatePval = NULL, annotateTop = TRUE, ...) 
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
  }
  else {
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
    xlabel = paste("Chromosome", unique(d$CHR), "position")
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
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
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
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
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
    with(d.highlight, points(pos, logp, col = "green3", pch = 20, 
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
      calibrate::textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
             #labs = topSNPs$SNP, cex = 0.8, ...)
             #labs = topHits$Overlapped.genes, cex = 1.3, ...)
             labs = topSNPs$SNP, cex = 1, font=2, ...)
    }
  }
  par(xpd = FALSE)
}
