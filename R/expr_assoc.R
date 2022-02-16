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
#' @details
#' Association testing between CNV regions and RNA-seq read counts is carried
#' out using edgeR, which applies generalized linear models (GLMs) based on the
#' negative-binomial distribution while incorporating normalization factors for
#' different library sizes.
#'
#' In the case of only one CN state deviating from 2n for a CNV region under
#' investigation, this reduces to the classical 2-group comparison.
#' For more than two states (e.g. 0n, 1n, 2n), edgeR’s ANOVA-like test is applied
#' to test all deviating groups for significant expression differences relative
#' to 2n.
#'
#' To avoid artificial effects due to low expression of a gene or insufficient
#' sample size in deviating groups, it is typically recommended to exclude from
#' the analysis (i) genes with fewer than r reads per million reads mapped
#' (cpm, counts per million) in the maximally expressed sample group,
#' and (ii) CNV regions with fewer than s samples in a group deviating from 2n.
#' Use the \code{min.cpm} and \code{min.samples} arguments, respectively.
#'
#' When testing local effects (adjacent or coinciding genes of a CNV region),
#' suitable thresholds for candidate discovery are r = 3, s = 4, and a nominal
#' significance level of 0.05; as such effects have a clear biological indication
#' and the number of genes tested is typically small.
#'
#' For distal effects (i.e. when testing genes far away from a CNV region)
#' more stringent thresholds such as r = 20 and s = 10 for distal effects in
#' conjunction with multiple testing correction using a conservative adjusted
#' significance level such as 0.01 is typically recommended (due to power
#' considerations and to avoid detection of spurious effects).
#'
#' @param cnvrs A \code{\linkS4class{GRanges}} or character object containing
#' the summarized CNV regions as e.g. obtained with \code{\link{populationRanges}}.
#' Alternatively, the assay name if the 'data' argument is provided.
#' @param calls Either a \code{\linkS4class{GRangesList}} or
#' \code{\linkS4class{RaggedExperiment}} storing the individual CNV calls for
#' each sample. Alternatively, the assay name if 'data' is provided.
#' @param rcounts A \code{\linkS4class{RangedSummarizedExperiment}} or
#' character name storing either the raw RNA-seq read counts in a rectangular
#' fashion (genes x samples). Alternatively, the assay name if 'data' is provided.
#' @param data (optional) A \code{MultiAssayExperiment} object 
#' with `cnvrs`, `calls`, and `rcounts` arguments corresponding to assay names.
#' @param window Numeric or Character. Size of the genomic window in base pairs
#' by which each CNV region is extended up- and downstream. This determines which
#' genes are tested for each CNV region. Character notation is supported for
#' convenience such as "100kbp" (same as 100000) or "1Mbp" (same as 1000000).
#' Defaults to \code{"1Mbp"}. Can also be set to \code{NULL} to test against all
#' genes included in the analysis.
#' @param multi.calls A function. Determines how to summarize the
#' CN state in a CNV region when there are multiple (potentially conflicting)
#' calls for one sample in that region. Defaults to \code{.largest}, which
#' assigns the CN state of the call that covers the largest part of the CNV
#' region tested. A user-defined function that is passed on to
#' \code{\link{qreduceAssay}} can also be provided for customized behavior.
#' @param min.samples Integer. Minimum number of samples with at least one call
#' overlapping the CNV region tested. Defaults to 10. See details.
#' @param de.method Character. Differential expression method.
#' Defaults to \code{"edgeR"}.
#' @param padj.method Character. Method for adjusting p-values to multiple testing.
#' For available methods see the man page of the function \code{\link{p.adjust}}.
#' Defaults to \code{"BH"}.
#' @param filter.by.expr Logical. Include only genes with
#' sufficiently large counts in the DE analysis? If TRUE, excludes genes not 
#' satisfying a minimum number of read counts across samples using the 
#' \code{\link{filterByExpr}} function from the edgeR package.
#' Defaults to TRUE.
#' @param verbose Logical. Display progress messages? Defaults to \code{FALSE}.
#' @return A \code{\linkS4class{DataFrame}} containing measures of association for 
#' each CNV region and each gene tested in the genomic window around the CNV region.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{findOverlaps}} to find overlaps between sets of genomic
#' regions,
#'
#' \code{\link{qreduceAssay}} to summarize ragged genomic location data
#' in defined genomic regions,
#'
#' \code{\link{glmQLFit}} and \code{\link{glmQLFTest}} to conduct negative
#' binomial generalized linear models for RNA-seq read count data.
#'
#' @references Geistlinger et al. (2018) Widespread modulation of gene expression
#' by copy number variation in skeletal muscle. Sci Rep, 8(1):1399.
#'
#' @examples
#'
#' # (1) CNV calls
#' states <- sample(c(0,1,3,4), 17, replace=TRUE)
#' calls <- GRangesList(
#'      sample1 = GRanges( c("chr1:1-10", "chr2:15-18", "chr2:25-34"), state=states[1:3]),
#'      sample2 = GRanges( c("chr1:1-10", "chr2:11-18" , "chr2:25-36"), state=states[4:6] ),
#'      sample3 = GRanges( c("chr1:2-11", "chr2:14-18", "chr2:26-36"), state=states[7:9] ),
#'      sample4 = GRanges( c("chr1:1-12", "chr2:18-35" ), state=states[10:11] ),
#'      sample5 = GRanges( c("chr1:1-12", "chr2:11-17" , "chr2:26-34"), state=states[12:14] ) ,
#'      sample6 = GRanges( c("chr1:1-12", "chr2:12-18" , "chr2:25-35"), state=states[15:17] )
#' )
#'
#' # (2) summarized CNV regions
#' cnvrs <- populationRanges(calls, density=0.1)
#'
#' # (3) RNA-seq read counts
#' genes <- GRanges(c("chr1:2-9", "chr1:100-150", "chr1:200-300",
#'                    "chr2:16-17", "chr2:100-150", "chr2:200-300", "chr2:26-33"))
#' y <- matrix(rnbinom(42,size=1,mu=10),7,6)
#' names(genes) <- rownames(y) <- paste0("gene", 1:7)
#' colnames(y) <- paste0("sample", 1:6)
#'
#' library(SummarizedExperiment)
#' rse <- SummarizedExperiment(assays=list(counts=y), rowRanges=granges(genes))
#'
#' # (4) perform the association analysis
#' res <- cnvEQTL(cnvrs, calls, rse, min.samples=1, filter.by.expr=FALSE)
#'
#' @export
cnvEQTL <- function(cnvrs, calls, rcounts, data,
    window="1Mbp", multi.calls=.largest,
    min.samples=10, de.method=c("edgeR", "limma"),
    padj.method="BH", filter.by.expr=TRUE, verbose=FALSE)
{
    if (!missing(data)) {
        stopifnot(is(data, "MultiAssayExperiment"))
        if (!requireNamespace("MultiAssayExperiment"))
            stop(
                paste0("Please install the 'MultiAssayExperiment'",
                " package to use 'cnvEQTL'")
            )

    }
    # sanity checks
    stopifnot(
        is(cnvrs, "GRanges") || is.character(cnvrs),
        is(calls, "GRangesList") ||
            is(calls, "RaggedExperiment") || is.character(calls),
        is(rcounts, "RangedSummarizedExperiment") || is.character(rcounts)
    )

    if (!is.character(calls))
        calls <- as(calls, "RaggedExperiment")

    if (!missing(data)) {
        ranges <- granges(rowRanges(data[[cnvrs]]))
        data <- data[, , names(data) != cnvrs]
        data <- as(data, "MatchedAssayExperiment")
        calls <- data[[calls]]
        rcounts <- data[[rcounts]]
        cnvrs <- ranges
    } else {
        # consider samples in cnv AND expression data
        sampleIds <- sort(intersect(colnames(calls), colnames(rcounts)))
        if(verbose) message(paste("Restricting analysis to",
                        length(sampleIds), "intersecting samples"))
        stopifnot(length(sampleIds) > 0)
        calls <- calls[,sampleIds]
        rcounts <- rcounts[,sampleIds]
    }
    # filter and norm RNA-seq data
    if(verbose) message("Preprocessing RNA-seq data ...")
    y <- .preprocRnaSeq(SummarizedExperiment::assay(rcounts), filter.by.expr)
    rcounts <- rcounts[rownames(y),]

    # determine states
    if(verbose) message("Summarizing per-sample CN state in each CNV region")
    cnv.states <- .getStates(cnvrs, calls, 
        multi.calls, min.samples, verbose=verbose)
    cnvrs <- IRanges::subsetByOverlaps(cnvrs,
                GenomicRanges::GRanges(rownames(cnv.states)), type="equal")

    # determine genes to test for each CNV region
    if(!is.null(window))
    {
        ecnvrs <- .extendRegions(cnvrs, window = window)
        olaps <- GenomicRanges::findOverlaps(rowRanges(rcounts), ecnvrs)
        cgenes <- split(S4Vectors::queryHits(olaps), S4Vectors::subjectHits(olaps))

        ind <- as.integer(names(cgenes))
        cnv.states <- cnv.states[ind,]
        cnvrs <- cnvrs[ind]
    }
    else
    {
        grid <- seq_len(nrow(rcounts))
        cgenes <- lapply(seq_len(length(cnvrs)), function(i) grid)
    }

    nr.cnvrs <- length(cgenes)
    if(verbose) message(paste("Analyzing", nr.cnvrs,
                                "regions with >=1 gene in the given window"))
    

    #TODO: BiocParallel
    de.method <- match.arg(de.method)
    cnvr.grid <- seq_len(nr.cnvrs)
    res <- lapply(cnvr.grid,
        function(i)
        {
            if(verbose) message(paste(i, "of", nr.cnvrs))
            r <- .testCnvExpr(y, cgenes[[i]], cnv.states[i,],
                                min.state.freq=min.samples,
                                de.method=de.method)
            return(r)
        })
    lens <- rep(cnvr.grid, lengths(cgenes))
    res <- do.call(plyr::rbind.fill, res)
    padj <- stats::p.adjust(res[,"PValue"], method = padj.method)
    res <- data.frame(res, AdjPValue = padj)
    res <- .formatResult(res, rcounts, lens)
    names(res) <- rownames(cnv.states)
    return(res)
}

#' Plot EQTL region
#'
#' Illustrates differential expression of genes in the neighborhood of a CNV.
#'
#' @param cnvr A \code{\linkS4class{GRanges}} of length 1, containing the genomic
#' coordinates of the CNV region of interest.
#' @param genes \code{\linkS4class{GRanges}} containing genes in the neighborhood
#' of the CNV region of interest.
#' @param genome Character. A valid UCSC genome assembly ID such as 'hg19' or 'bosTau6'.
#' @param cn Character. Copy number state of interest.
#' @return None. Plots to a graphics device.
#'
#' @author Ludwig Geistlinger
#' @seealso \code{Gviz::plotTracks}
#' 
#' @examples
#'
#' # CNV region of interest
#' cnvr <- GRanges("chr1:7908902-8336254")
#'
#' # Two genes in the neighborhood
#' genes <- c("chr1:8021714-8045342:+", "chr1:8412464-8877699:-")
#' names(genes) <- c("PARK7", "RERE")
#' genes <- GRanges(genes)
#' 
#' # Annotate differential expression for 1-copy loss
#' genes$logFC.CN1 <- c(-0.635, -0.728)
#' genes$AdjPValue <- c(8.29e-09, 1.76e-08) 
#'
#' # plot
#' plotEQTL(cnvr, genes, genome="hg19", cn="CN1")
#'
#' @export
plotEQTL <- function(cnvr, genes, genome, cn="CN1")
{
    if (!requireNamespace("Gviz", quietly = TRUE))
        stop(paste("Required package \'Gviz\' not found.", 
                    "Use \'BiocManager::install(\"Gviz\") to install it."))

    colM <- "gray24"
    colB <- "gray94"
    chr <- as.character(seqnames(cnvr))

    itrack <- Gviz::IdeogramTrack(genome=genome, chr=chr, fontsize=15)
    gtrack <- Gviz::GenomeAxisTrack(littleTicks=TRUE, fontsize=15)    

    cnvrTrack <- Gviz::AnnotationTrack(cnvr, name="CNVR", group=" ", fill="red",
        col.title=colM, col.axis=colM, background.title=colB, cex.title=0.8)

    pgenes <- subset(genes, strand=="+")
    plusTrack <- Gviz::AnnotationTrack(pgenes, name="", 
            #group=names(genes), just.group="right", 
            cex=0.8, rotation.item=45,
            id=names(pgenes), featureAnnotation="id", fontcolor.feature="darkblue",
            col.title=colM, col.axis=colM, background.title=colB)
    
    mgenes <- subset(genes, strand == "-")
    minusTrack <- Gviz::AnnotationTrack(mgenes, name="", 
            #group=names(genes), just.group="right", 
            cex=0.8, rotation.item=45,
            id=names(mgenes), featureAnnotation="id", fontcolor.feature="darkblue",
            col.title=colM, col.axis=colM, background.title=colB)

    GenomicRanges::strand(genes) <- "*"
    lcol <- paste("logFC", cn, sep = ".")
    lcol <- S4Vectors::mcols(genes)[[lcol]]
    fctrack <- Gviz::DataTrack(genes, data=lcol, type="h",  
                                name="log2FC", cex.title=1, cex.axis=1,
                                font.axis=2, col.title=colM, col.axis=colM, 
                                background.title=colB)           
    
    ps <- -log10(genes$AdjPValue) 
    ptrack <- Gviz::DataTrack(genes, data=ps, type="h",  
                                name="-log10 adjP", cex.title=1, cex.axis=1,
                                font.axis=2, col.title=colM, col.axis=colM, 
                                background.title=colB)           
       
    Gviz::plotTracks(c(itrack, gtrack, cnvrTrack, 
                        plusTrack, minusTrack, fctrack, ptrack))
}

# test a single cnv region
.testCnvExpr <- function(y, cgenes, states, min.state.freq=10,
    de.method=c("limma", "edgeR"))
{
    # form groups according to CNV states
    state.freq <- table(states)
    nr.states <- length(state.freq)
    too.less.samples <- state.freq < min.state.freq
    if(any(too.less.samples))
    {
        too.less.states <- names(state.freq)[too.less.samples]
        too.less.states <- as.integer(too.less.states)
        ind <- !(states %in% too.less.states)
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

    # create result table
    fc.cols <- grep("^logFC", colnames(res), value=TRUE)
    rel.cols <- c(fc.cols, "PValue")
    ind <- order(res[,"PValue"])
    de.tbl <- res[ind, rel.cols]
    de.tbl <- cbind(rownames(de.tbl), de.tbl, stringsAsFactors=FALSE)
    rownames(de.tbl) <- NULL
    colnames(de.tbl)[1] <- "Gene"

    # remap fc cols
    fc.cols <- grep("^logFC", colnames(de.tbl))
    sn <- levels(group)[fc.cols]
    smap <- c("CN0", "CN1", "CN3", "CN4")
    names(smap) <- c("s-2", "s-1", "s+1", "s+2")
    colnames(de.tbl)[fc.cols] <- paste("logFC", smap[sn], sep=".") 
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

.getStates <- function(cnvrs, calls, multi.calls=.largest, min.samples=10, verbose=FALSE)
{
    #TODO: additional pre-defined options for multi.calls
    stopifnot(is.function(multi.calls))

    cnv.states <- RaggedExperiment::qreduceAssay(calls, query=cnvrs,
                    simplifyReduce=multi.calls, background=2)

    # exclude cnv regions with not at least one gain/loss state with >= min.samples
    tab <- apply(cnv.states, 1, table)
    .isTestable <- function(states)
    {
        cond1 <- length(states) > 1
        cond2 <- any(states[names(states) != "2"] >= min.samples)
        cond1 && cond2
    }
    ind <- vapply(tab, .isTestable, logical(1))
    nr.excl <- sum(!ind)
    if (nr.excl && verbose)
        message(paste("Excluding", nr.excl,
            "cnvrs not satisfying min.samples threshold"))
    if(length(cnvrs) == nr.excl)
        stop("No CNV region satisfying min.samples threshold")
    cnv.states[ind, ]
}

.largest <- function(scores, ranges, qranges)
{
    ind <- BiocGenerics::which.max(BiocGenerics::width(ranges))
    vapply(seq_along(scores), function(i) scores[[i]][ind[i]], scores[[1]][1])
}


.weightedmean <- function(scores, ranges, qranges)
{
    isects <- pintersect(ranges, qranges)
    sum(scores * width(isects)) / sum(width(isects))
}

# filter low exprs
.preprocRnaSeq <- function(rcounts, filter.by.expr=TRUE)
{
    if(filter.by.expr)
    {
        keep <- suppressMessages(edgeR::filterByExpr(rcounts))
        rcounts <- rcounts[keep,]
    }
    y <- edgeR::DGEList(counts=rcounts)
    y <- edgeR::calcNormFactors(y)
    return(y)
}
 
# format output as GRangesList
.formatResult <- function(res, rcounts, lens)
{
    colnames(res) <- sub("logFC.NA", "logFC.CN3", colnames(res))
    fc.ind <- grep("^logFC", colnames(res), value=TRUE)
    fc.ind <- sort(fc.ind)
    res <- res[,c("Gene", fc.ind, "PValue", "AdjPValue")] 

    g <- res[,"Gene"]
    g <- SummarizedExperiment::rowRanges(rcounts)[g]
    GenomicRanges::mcols(g) <- res[,2:ncol(res)]
    
    grl <- split(g, lens)
    grl <- GenomicRanges::GRangesList(grl)
    return(grl)
}
