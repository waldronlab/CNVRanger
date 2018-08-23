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
#'
#' @param cnvrs A \code{\linkS4class{GRanges}} object containing the summarized
#' CNV regions as e.g. obtained with \code{\link{populationRanges}}.
#' @param calls A \code{\linkS4class{GRangesList}} or a 
#' \code{linkS4class{RaggedExperiment}} storing the individual CNV calls for 
#' each sample
#' @param rcounts A \code{\linkS4class{RangedSummarizedExperiment}} storing the 
#' raw RNA-seq read counts for each sample.
#'
cnvExprAssoc <- function(cnvrs, calls, rcounts, 
    window="1Mbp", multi.calls="largest", 
    min.samples=10, min.cpm=2, padj.method="BH", verbose=TRUE)
{
    # sanity checks
    stopifnot(is(cnvrs, "GRanges"), 
                is(calls, "GRangesList") || is(calls, "RaggedExperiment"),
                is(rcounts, "RangedSummarizedExperiment"))

    stopifnot(is.integer(RaggedExperiment::assay(calls)), 
                is.integer(SummarizedExperiment::assay(rcounts)))

    if(is(calls, "GRangesList")) 
        calls <- RaggedExperiment::RaggedExperiment(grl)

    # consider samples in cnv AND expression data
    sampleIds <- sort(intersect(colnames(calls), colnames(rcounts)))
    if(verbose) message(paste("Restricting analysis to", 
                    length(sampleIds), "intersecting samples"))
    stopifnot(length(sampleIds) > 0)
    calls <- calls[,sampleIds]
    rcounts <- rcounts[,sampleIds]         
   
    # filter and norm RNA-seq data
    if(verbose)
    { 
        message("Preprocessing RNA-seq data ...")
        message("Calculating normalization factors and estimating dispersion ...")
    }
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
    res <- lapply(seq_len(nr.cnvrs), 
        function(i)
        { 
            if(verbose) message(paste(i, "of", nr.cnvrs))
            if(i %in% c(4,5,8)) return(NULL)
            ind <- cgenes[[i]]
            yi <- edgeR::`[.DGEList`(y, ind, )
            r <- .testCnvExpr(yi, cnv.states[i,],
                                min.state.freq=min.samples, 
                                padj.method=padj.method)
            return(r)
        })
    names(res) <- rownames(cnv.states)
    return(res)
}

# test a single cnv region
.testCnvExpr <- function(y, states, min.state.freq=10, padj.method="BH")
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
    fit <- edgeR::glmQLFit(y, design, robust=TRUE)
    qlf <- edgeR::glmQLFTest(fit, coef=2:nr.states)
    fc.cols <- grep("^logFC", colnames(qlf$table), value=TRUE)
    rel.cols <- c(fc.cols, "PValue")
    ind <- order(qlf$table[,"PValue"])
    de.tbl <- qlf$table[ind, rel.cols]    
    
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
    y <- edgeR::estimateDisp(y, robust=TRUE)
    return(y)
}

