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
    min.samples=10, min.cpm=5, padj.method="BH")
{
    # sanity checks
    stopifnot(is(cnvrs, "GRanges"), 
                is(calls, "GRangesList") || is(calls, "RaggedExperiment"),
                is(is(rcounts, "RangedSummarizedExperiment"))

    stopifnot(is.integer(assay(calls)), is.integer(assay(rcounts)))

    if(is(calls, "GRangesList")) 
        calls <- GenomicRanges::makeGRangesListFromDataFrame(calls, 
                        split.field="sampleId", keep.extra.columns=TRUE)

    # consider samples in cnv AND expression data
    sampleIds <- sort(intersect(colnames(calls), colnames(rcounts)))
    calls <- calls[,sampleIds]
    rcounts <- rcounts[,sampleIds]         
   
    # filter and norm RNA-seq data
    y <- .preprocRnaSeq(assay(rcounts))

    # restrict RNA-seq data to genes in window for each cnv 

    # determine states 
    cnv.states <- RaggedExperiment::qreduceAssay(calls, query=cnvrs, 
                    simplifyReduce=.largest, background=2)

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
    cnv.states <- cnv.states[ind]
    cnvrs <- cnvrs[ind]   
 
    apply(cnv.states, 1, testCnvExpr, y=, min.state.freq=)

}

# test a single cnv region
testCnvExpr <- function(y, states, min.state.freq=10, padj.method="BH")
{
    # form groups according to CNV states
    state.freq <- table(states)
    too.less.samples <- state.freq < min.state.freq
    if(sum(too.less.samples))
    {
        too.less.states <- names(state.freq)[too.less.samples]
        too.less.states <- as.integer(too.less.states)
        ind <- states != too.less.states
        nr.states <- length(state.freq) - length(too.less.states) 
    }
    stopifnot(nr.states > 1) 
    states <- states[ind]
    y <- y[,ind]

    # design
    s2g <- c("B","C","A","D", "E")
    sort(ifelse(s-2 <= 0, paste0("s-", abs(s-2)), paste0("s+", abs(s-2))))
    group <- s2g[states + 1]
    group <- as.factor(group)
    y$group <- group
    f <- stats::formula(paste0("~", "group"))
    design <- stats::model.matrix(f)

    # test
    y <- edgeR::estimateDisp(y, design, robust=TRUE)
    fit <- edgeR::glmQLFit(y, design, robust=TRUE)
    qlf <- edgeR::glmQLFTest(fit, coef=2:nr.states)
    fc.cols <- grep("^logFC", colnames(qlf$table), value=TRUE)
    rel.cols <- c(fc.cols, "PValue")
    ind <- order(de.tbl[,"PValue"])
    de.tbl <- qlf$table[ind, rel.cols]    
    
    # multiple testing
    padj <- stats::p.adjust(de.tbl[,"PValue"], method=padj.method)
    de.tbl <- cbind(de.tbl, padj)
    colnames(de.tbl)[ncol(de.tbl)] <- "AdjPValue"
    return(de.tbl)
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

