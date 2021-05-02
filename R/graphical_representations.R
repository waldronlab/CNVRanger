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
#' @param all.paths Object returned from \code{CreateFolderTree} function with 
#' the working folder tree
#' @param regions \code{\linkS4class{GRanges}} as returned by \code{\link{cnvGWAS}}  
#' @param chr.size.order \code{\link{data.frame}} with two columns: (i) 'chr': chromosome names (character),
#' and (ii) 'size': length of the chromosomes in bp (integer). A \code{\linkS4class{GRanges}} containing one
#' chromosome per range can be used instead (the chromosomes should be in the expected order).
#' @param plot.pdf Logical plot a to pdf file
#' @return Plots to graphics device.
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples 
#' 
#' # Load phenotype-CNV information
#' data.dir <- system.file("extdata", package="CNVRanger")
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
#' order.chrs <- c(1:24, "25LG1", "25LG2", 27:28, "LGE22", "1A", "4A")
#'
#' # Chromosome sizes
#' chr.size.file <- file.path(data.dir, "Parus_major_chr_sizes.txt")
#' chr.sizes <- scan(chr.size.file)
#' chr.size.order <- data.frame(chr=order.chrs, sizes=chr.sizes, stringsAsFactors=FALSE)
#'
#' # Plot Manhatthan to a pdf within the 'Results' workfolder
#' plotManhattan(all.paths, segs.pvalue.gr, chr.size.order)
#'
#' @export
#'  

plotManhattan <- function(all.paths, regions, chr.size.order, plot.pdf = FALSE) {
    
    ## Check inputs to avoid errors
    stopifnot(length(all.paths)==2)
    stopifnot(is(regions, "GRanges"))
    
    if(is(chr.size.order, "GRanges"))
        chr.size.order <- data.frame(chr=as.character(seqnames(regions)), 
                                   sizes=end(chr.size.order),
                                   stringsAsFactors=FALSE)
    
    
    stopifnot(all(!duplicated(chr.size.order$chr)))
    
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
    
    
    gwasResults <- data.frame(
        CHR = as.character(GenomicRanges::seqnames(regions)), 
        BP = GenomicRanges::start(regions), 
        SNP = regions$SegName, P = regions$MinPvalueAdjusted, 
        stringsAsFactors = FALSE)
    
    gwasResults <- rbind(gwasResults, chr.all)
    
    ### turn chr to numeric
    for (cn in chr.size.order$chr.numeric) {
        chr.str <- paste0("\\<", chr.size.order[cn, 1], "\\>")
        gwasResults$CHR <- gsub(chr.str, chr.size.order[cn,3], gwasResults$CHR)
    }
    
    gwasResults$CHR <- as.integer(gwasResults$CHR)
    gwasResults$BP <- as.integer(gwasResults$BP)
    gwasResults$P <- as.numeric(gwasResults$P)
    ind <- do.call(order, gwasResults[, c("CHR", "BP")])
    gwasResults <- gwasResults[ind, ]
    
    if(plot.pdf) pdf(file.path(all.paths[3], paste0(phen.name, "-Manhattan.pdf")))
    qqman::manhattan(gwasResults, chr = "CHR", bp = "BP", snp = "SNP", p = "P", 
                        chrlabs = chr.size.order$chr, 
                        suggestiveline = -log10(0.1), genomewideline = -log10(0.05))
    if(plot.pdf) dev.off()
}


# HELPER - Plot QQ-plot using non-log p-values Originally published in:
# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R By
# Matthew Flickinger @import lattice

.qqunifPlot <- function(pvalues, should.thin = TRUE, thin.obs.places = 2, thin.exp.places = 2, 
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


