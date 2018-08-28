############################################################
# 
# author: Ludwig Geistlinger
# date: 2018-08-21 12:45:44
# 
# descr: CNV-expression association analysis
# 
############################################################

#' CNV-expression association analysis
#' 
#' Testing CNV regions for effects on the expression level of genes in defined
#' genomic windows.  
#'  
#' Details here.
#'
#'
#' @param cnvrs A \code{\linkS4class{GRanges}} object containing the summarized
#' CNV regions as e.g. obtained with \code{\link{populationRanges}}.
#' @param calls A \code{\linkS4class{GRangesList}} or a 
#' \code{linkS4class{RaggedExperiment}} storing the individual CNV calls for 
#' each sample.
#' @param rcounts A \code{\linkS4class{RangedSummarizedExperiment}} storing the 
#' raw RNA-seq read counts in a rectangular fashion (genes x samples).
#' @param window Defaults to \code{"1Mbp"}.
#' @param multi.calls Character or a function. Defaults to \code{"largest"}.
#' @param min.samples Integer. Defaults to 10.
#' @param min.cpm Integer. Should genes not satisfying a minimum counts-per-million 
#' (cpm) threshold be excluded from the analysis? This is typically recommended. 
#' See the edgeR vignette for details. The default filter is to exclude genes 
#' with cpm < 2 in more than half of the samples.
#' @param de.method Character. Differential expression method. 
#' Defaults to \code{"edgeR"}.
#' @param padj.method Character. Method for adjusting p-values to multiple testing.
#' For available methods see the man page of the function \code{\link{p.adjust}}. 
#' Defaults to \code{"BH"}.
#' @param verbose Logical. Display progress messages? Defaults to \code{FALSE}.
#' @return TODO
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{findOverlaps}} to find overlaps between sets of genomic 
#' regions, \code{\link{qreduceAssay}} to summarize ragged genomic location data 
#' in defined genomic regions, \code{\link{glmQLFit}} and \code{\link{glmQLFTest}}
#' to conduct negative binomial generalized linear models for RNA-seq read 
#' count data.
#' @references Geistlinger et al. (2018) Widespread modulation of gene expression
#' by copy number variation in skeletal muscle. Sci Rep, 8(1):1399.
#' 
#' @examples 
#'
#' states <- sample(c(0,1,3,4), 17, replace=TRUE)
#' 
#' calls <- GRangesList(
#'      sample1 = GRanges( c("chr1:1-10", "chr2:15-18", "chr2:25-34"), state=states[1:3]),
#'      sample2 = GRanges( c("chr1:1-10", "chr2:11-18" , "chr2:25-36"), state=states[4:6] ),
#'      sample3 = GRanges( c("chr1:2-11", "chr2:14-18", "chr2:26-36"), state=states[7:9] ),
#'      sample4 = GRanges( c("chr1:1-12", "chr2:18-35" ), state=states[10:11] ),
#'      sample5 = GRanges( c("chr1:1-12", "chr2:11-17" , "chr2:26-34"), state=states[12:14] ) ,
#'      sample6 = GRanges( c("chr1:1-12", "chr2:12-18" , "chr2:25-35"), state=states[15:17] )
#' ) 
#'
#' cnvrs <- populationRanges(calls, density=0.1)
#'
#' 
#'
#' @export
cnvExprAssoc <- function(cnvrs, calls, rcounts, 
    window="1Mbp", multi.calls="largest", 
    min.samples=10, min.cpm=2, de.method=c("edgeR", "limma"),
    padj.method="BH", verbose=FALSE)
{
    # sanity checks
    stopifnot(is(cnvrs, "GRanges"), 
                is(calls, "GRangesList") || is(calls, "RaggedExperiment"),
                is(rcounts, "RangedSummarizedExperiment"))

    if(is(calls, "GRangesList")) 
        calls <- RaggedExperiment::RaggedExperiment(calls)

    stopifnot(is.integer(RaggedExperiment::assay(calls)), 
                is.integer(SummarizedExperiment::assay(rcounts)))

    # consider samples in cnv AND expression data
    sampleIds <- sort(intersect(colnames(calls), colnames(rcounts)))
    if(verbose) message(paste("Restricting analysis to", 
                    length(sampleIds), "intersecting samples"))
    stopifnot(length(sampleIds) > 0)
    calls <- calls[,sampleIds]
    rcounts <- rcounts[,sampleIds]         
   
    # filter and norm RNA-seq data
    if(verbose) message("Preprocessing RNA-seq data ...")
    y <- .preprocRnaSeq(SummarizedExperiment::assay(rcounts))
    rcounts <- rcounts[rownames(y),]

    # determine states
    if(verbose) message("Summarizing per-sample CN state in each CNV region")
    cnv.states <- .getStates(cnvrs, calls, multi.calls, min.samples) 
    cnvrs <- IRanges::subsetByOverlaps(cnvrs, 
                GenomicRanges::GRanges(rownames(cnv.states)))   
   
    # determine genes to test for each CNV region 
    ecnvrs <- .extendRegions(cnvrs, window=window)
    olaps <- GenomicRanges::findOverlaps(rowRanges(rcounts), ecnvrs) 
    cgenes <- split(S4Vectors::queryHits(olaps), S4Vectors::subjectHits(olaps))
    
    nr.cnvrs <- length(cgenes)
    if(verbose) message(paste("Analyzing", nr.cnvrs, 
                                "regions with >=1 gene in the given window"))
    ind <- as.integer(names(cgenes))
    cnv.states <- cnv.states[ind,]
    cnvrs <- cnvrs[ind]

    #TODO: BiocParallel
    de.method <- match.arg(de.method)
    res <- lapply(seq_len(nr.cnvrs), 
        function(i)
        { 
            if(verbose) message(paste(i, "of", nr.cnvrs))
            r <- .testCnvExpr(y, cgenes[[i]], cnv.states[i,],
                                min.state.freq=min.samples, 
                                de.method=de.method,
                                padj.method=padj.method)
            return(r)
        })
    names(res) <- rownames(cnv.states)
    return(res)
}

# test a single cnv region
.testCnvExpr <- function(y, cgenes, states, min.state.freq=10, 
    de.method=c("limma", "edgeR"), padj.method="BH")
{
    # form groups according to CNV states
    state.freq <- table(states)
    nr.states <- length(state.freq)
    too.less.samples <- state.freq < min.state.freq
    if(sum(too.less.samples))
    {
        too.less.states <- names(state.freq)[too.less.samples]
        too.less.states <- as.integer(too.less.states)
        ind <- states != too.less.states
        nr.states <- length(state.freq) - length(too.less.states) 
        stopifnot(nr.states > 1) 
        states <- states[ind]
        y <- edgeR::`[.DGEList`(y, , ind)
    }
    
    # design
    s <- states - 2
    s <- paste0("s", ifelse(s <= 0, "-", "+"), abs(s))
    group <- as.factor(s)
    y$samples$group <- group
    f <- stats::formula(paste0("~", "group"))
    design <- stats::model.matrix(f)

    # test
    de.method <- match.arg(de.method)
    if(de.method == "limma")
    {
        y <- limma::voom(y, design)
        fit <- limma::lmFit(y, design)

        # restrict to cgenes
        fit <- limma::`[.MArrayLM`(fit, cgenes,)
        fit <- limma::eBayes(fit)
        res <- limma::topTable(fit, number=length(cgenes), coef=2:nr.states,
                                sort.by="none", adjust.method="none")
        colnames(res) <- sub("\\.", "", colnames(res))
    }
    else
    {
        y <- edgeR::estimateDisp(y, design, robust=TRUE)
        fit <- edgeR::glmQLFit(y, design, robust=TRUE)

        # restrict to cgenes
        fit <- edgeR::`[.DGEGLM`(fit, cgenes,)
        res <- edgeR::glmQLFTest(fit, coef=2:nr.states)
        res <- res$table
    }
    fc.cols <- grep("^logFC", colnames(res), value=TRUE)
    rel.cols <- c(fc.cols, "PValue")
    ind <- order(res[,"PValue"])
    de.tbl <- res[ind, rel.cols]    
    
    # multiple testing
    padj <- stats::p.adjust(de.tbl[,"PValue"], method=padj.method)
    de.tbl <- cbind(de.tbl, padj)
    colnames(de.tbl)[ncol(de.tbl)] <- "AdjPValue"
    return(de.tbl)
}

.extendRegions <- function(regions, window="1Mbp")
{
    stopifnot(is.character(window) || is.numeric(window))
    if(is.character(window)) window <- .window2integer(window)
    nstart <- GenomicRanges::start(regions) - window
    nend <- GenomicRanges::end(regions) + window
    GenomicRanges::start(regions) <- pmax(1, nstart) 
    GenomicRanges::end(regions) <- nend
    return(regions)
}

.window2integer <- function(w)
{
    stopifnot(grepl("[0-9]+[GMk]*bp$", w))
    w <- sub("bp$", "", w)
    unit <- 0
    if(grepl("[GMk]$", w))
    {
        n <- nchar(w)
        unit <- substring(w, n, n)
        unit <- switch(unit, k=3, M=6, G=9)
        w <- substring(w, 1, n-1)
    }
    w <- as.integer(w) * 10^unit
    return(w)
}

.getStates <- function(cnvrs, calls, multi.calls="largest", min.samples=10)
{
    #TODO: additional pre-defined options for multi.calls
    stopifnot(is.character(multi.calls) || is.function(multi.calls))
    if(is.character(multi.calls))
    {
        stopifnot(multi.calls == "largest")
        multi.calls <- .largest
    }

    cnv.states <- RaggedExperiment::qreduceAssay(calls, query=cnvrs, 
                    simplifyReduce=multi.calls, background=2)

    # exclude cnv regions with not at least one gain/loss state with >= min.samples 
    tab <- apply(cnv.states, 1, table)
    .isTestable <- function(states)
    { 
        cond1 <- length(states) > 1
        cond2 <- any(states[names(states) != "2"] >= min.samples)
        res <- cond1 && cond2
        return(res)
    } 
    ind <- vapply(tab, .isTestable, logical(1))
    nr.excl <- sum(!ind)
    if(nr.excl) message(paste("Excluding", nr.excl, 
                                "cnvrs not satisfying min.samples threshold"))
    stopifnot(length(cnvrs) > nr.excl)
    cnv.states <- cnv.states[ind,]
    return(cnv.states)
}

.largest <- function(scores, ranges, qranges) 
{
    return.type <- class(scores[[1]])
    default.value <- do.call(return.type, list(1))
    ind <- which.max(width(ranges))
    res <- vapply(seq_along(scores), 
           function(i) scores[[i]][ind[i]], default.value)
    return(res)
}

.weightedmean <- function(scores, ranges, qranges)
{
    isects <- pintersect(ranges, qranges)
    sum(scores * width(isects)) / sum(width(isects))
}

# filter low exprs
.preprocRnaSeq <- function(rcounts, min.cpm=2)
{ 
    rs <- rowSums(edgeR::cpm(rcounts) > min.cpm)
    keep <- rs >= ncol(rcounts) / 2
    rcounts <- rcounts[keep,]   
    y <- edgeR::DGEList(counts=rcounts) 
    y <- edgeR::calcNormFactors(y)
    return(y)
}

