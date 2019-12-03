###########################################################
# 
# author: Ludwig Geistlinger
# date: 2017-08-12 22:31:59
# 
# descr: Functionality for summarizing individual calls
#   across the population under study
#   e.g. defining CNV regions from CNV calls
# 
############################################################

#' Summarizing CNV ranges across a population
#'
#' In CNV analysis, it is often of interest to summarize individual calls across 
#' the population, (i.e. to define CNV regions), for subsequent association 
#' analysis with e.g. phenotype data. 
#'
#' @details
#' \itemize{
#' \item CNVRuler procedure that trims region margins based on regional density
#'
#' Trims low-density areas (usually <10\% of the total contributing individual 
#' calls within a summarized region).
#' 
#' An illustration of the concept can be found here: 
#' https://www.ncbi.nlm.nih.gov/pubmed/22539667 (Figure 1)
#'
#' \item Reciprocal overlap (RO) approach (e.g. Conrad et al., Nature, 2010)
#'
#' Reciprocal overlap of 0.51 between two genomic regions A and B:
#'
#' requires that B overlaps at least 51\% of A, 
#'   *and* that A also overlaps at least 51\% of B
#'
#' Approach:
#'
#' At the top level of the hierarchy, all contiguous bases overlapping at least 
#' 1bp of individual calls are merged into one region. 
#' Within each region, we further define reciprocally overlapping regions with 
#' the following algorithm:
#'
#' \itemize{
#' \item Calculate reciprocal overlap (RO) between all remaining calls.
#' \item Identify pair of calls with greatest RO. If RO > threshold, merge and 
#'       create a new CNV. If not, exit.
#' \item Continue adding unclustered calls to the region, in order of best overlap. 
#'       In order to add a call, the new call must have > threshold to all calls 
#'       within the region to be added. When no additional calls may be added, 
#'       move to next step.
#' \item If calls remain, return to 1. Otherwise exit.
#' }
#'
#' \item  GISTIC procedure (Beroukhim et al., PNAS, 2007) to identify recurrent 
#' CNV regions
#' 
#' GISTIC scores each CNV region with a G-score that is proportional to the total 
#' magnitude of CNV calls in each CNV region.
#' In addition, by permuting the locations in each sample, GISTIC determines 
#' the frequency with which a given score would be attained if the events were 
#' due to chance and therefore randomly distributed.
#' A significance threshold can then be used to determine scores / regions that 
#' are unlikely to occur by chance alone.  
#' }
#'
#' @param grl A \code{\linkS4class{GRangesList}}.
#' @param mode Character. Should population ranges be computed based on regional
#' density ("density") or reciprocal overlap ("RO"). See Details.
#' @param density Numeric. Defaults to 0.1. 
#' @param ro.thresh Numeric. Threshold for reciprocal overlap required for merging
#' two overlapping regions. Defaults to 0.5.
#' @param multi.assign Logical. Allow regions to be assigned to several region clusters?
#' Defaults to \code{FALSE}.
#' @param verbose Logical. Report progress messages? Defaults to \code{FALSE}.
#' @param min.size Numeric. Minimum size of a summarized region to be included. 
#' Defaults to 2 bp.
#' @param classify.ranges Logical. Should CNV frequency (number of samples 
#' overlapping the region) and CNV type (gain, loss, or both) be annotated?
#' Defaults to \code{TRUE}.
#' @param type.thresh Numeric. Required minimum relative frequency of each CNV 
#' type (gain / loss) to be taken into account when assigning CNV type to a region.
#' Defaults to 0.1. That means for a region overlapped by individual gain and 
#' loss calls that both types must be present in >10% of the corresponding samples
#' in order to be typed as 'both'. If gain or loss calls are present below the 
#' threshold they are ignored.    
#' @param est.recur Logical. Should recurrence of regions be assessed via a 
#' permutation test? Defaults to \code{FALSE}. See Details.  
#' @return A \code{\linkS4class{GRanges}} object containing the summarized
#' CNV ranges. 
#'
#' @references 
#' Kim et al. (2012) CNVRuler: a copy number variation-based case-control
#' association analysis tool. Bioinformatics, 28(13):1790-2.
#'
#' Conrad et al. (2010) Origins and functional impact of copy number variation 
#' in the human genome. Nature, 464(7289):704-12.
#'
#' Beroukhim et al. (2007) Assessing the significance of chromosomal aberrations
#' in cancer: methodology and application to glioma. PNAS, 104(50):20007-12.
#'
#' @author Ludwig Geistlinger, Martin Morgan
#' @seealso \code{\link{findOverlaps}}
#' 
#' @examples
#'
#' grl <- GRangesList(
#'      sample1 = GRanges( c("chr1:1-10", "chr2:15-18", "chr2:25-34") ),
#'      sample2 = GRanges( c("chr1:1-10", "chr2:11-18" , "chr2:25-36") ),
#'      sample3 = GRanges( c("chr1:2-11", "chr2:14-18", "chr2:26-36") ),
#'      sample4 = GRanges( c("chr1:1-12", "chr2:18-35" ) ),
#'      sample5 = GRanges( c("chr1:1-12", "chr2:11-17" , "chr2:26-34") ) ,
#'      sample6 = GRanges( c("chr1:1-12", "chr2:12-18" , "chr2:25-35") )
#' )
#'
#' # default as chosen in the original CNVRuler procedure
#' populationRanges(grl, density=0.1, classify.ranges=FALSE)
#'
#' # density = 0 merges all overlapping regions, 
#' # equivalent to: reduce(unlist(grl))
#' populationRanges(grl, density=0, classify.ranges=FALSE) 
#'
#' # density = 1 disjoins all overlapping regions, 
#' # equivalent to: disjoin(unlist(grl))
#' populationRanges(grl, density=1, classify.ranges=FALSE)
#'
#' # RO procedure
#' populationRanges(grl, mode="RO", ro.thresh=0.5, classify.ranges=FALSE)
#'
#' @export
populationRanges <- function(grl, mode=c("density", "RO"), 
    density=0.1, ro.thresh=0.5, multi.assign=FALSE, verbose=FALSE,
    min.size=2, classify.ranges=TRUE, type.thresh=0.1, est.recur=FALSE)
{
    mode <- match.arg(mode)
    if(mode == "density") pranges <- .densityPopRanges(grl, density)
    else pranges <- .roPopRanges(grl, ro.thresh, multi.assign, verbose)

    pranges <- pranges[GenomicRanges::width(pranges) >= min.size]
    if(classify.ranges) pranges <- .classifyRegs(pranges, unlist(grl), type.thresh)
    if(est.recur) pranges <- .estimateRecurrence(pranges, grl)
    return(pranges)
}

#' Plot recurrent CNV regions
#'
#' Illustrates summarized CNV regions along a chromosome.
#'
#' @param regs A \code{\linkS4class{GRanges}}. Typically the result of 
#' \code{\link{populationRanges}} with \code{est.recur=TRUE}.
#' @param genome Character. A valid UCSC genome assembly ID such as 'hg19' or 'bosTau6'.
#' @param chr Character. A UCSC-style chromosome name such as 'chr1'. 
#' @param pthresh Numeric. Significance threshold for recurrence. Defaults to 0.05.
#' @return None. Plots to a graphics device.
#'
#' @author Ludwig Geistlinger
#' @seealso \code{Gviz::plotTracks}
#' 
#' @examples
#'
#' # read in example CNV calls
#' data.dir <- system.file("extdata", package="CNVRanger")
#' call.file <- file.path(data.dir, "Silva16_PONE_CNV_calls.csv")
#' calls <- read.csv(call.file, as.is=TRUE)
#'
#' # store in a GRangesList
#' grl <- GenomicRanges::makeGRangesListFromDataFrame(calls, 
#'    split.field="NE_id", keep.extra.columns=TRUE)
#'
#' # summarize CNV regions
#' cnvrs <- populationRanges(grl, density=0.1, est.recur=TRUE)
#'
#' # plot
#' plotRecurrentRegions(cnvrs, genome="bosTau6", chr="chr1")
#'
#' @export
plotRecurrentRegions <- function(regs, genome, chr, pthresh=0.05)
{
    if (!requireNamespace("Gviz", quietly = TRUE))
        stop(paste("Required package \'Gviz\' not found.", 
                    "Use \'BiocManager::install(\"Gviz\") to install it."))

    colM <- "gray24"
    colB <- "gray94"
    tlevels <- c("gain", "loss", "both")

    chr.regs <- subset(regs, seqnames == chr)
    if(!length(chr.regs)) "No CNV regions to plot on the specified chromosome"

    itrack <- Gviz::IdeogramTrack(genome=genome, chr=chr, fontsize=15)
    gtrack <- Gviz::GenomeAxisTrack(littleTicks=TRUE, fontsize=15)

    gain.regs <- subset(chr.regs, type=="gain")
    tlist <- list()
    if(length(gain.regs))
    {
        gain.track <- Gviz::DataTrack(gain.regs, data=gain.regs$freq, 
                                type="h", groups=factor("gain", levels=tlevels), 
                                name="#samples", cex.title=1, cex.axis=1,
                                font.axis=2, col.title=colM, col.axis=colM, 
                                background.title=colB, legend=TRUE)
        tlist <- c(tlist, gain.track)
    }

    loss.regs <- subset(chr.regs, type=="loss")
    if(length(loss.regs))
    {
        loss.track <- Gviz::DataTrack(loss.regs, data=loss.regs$freq, 
                                type="h", groups=factor("loss", levels=tlevels), 
                                name="#samples", cex.title=1, cex.axis=1,
                                font.axis=2, col.title=colM, col.axis=colM, 
                                background.title=colB, legend=TRUE)
        tlist <- c(tlist, loss.track)
    }

    both.regs <- subset(chr.regs, type=="both")
    if(length(both.regs))
    {
        both.track <- Gviz::DataTrack(both.regs, data=both.regs$freq, 
                                type="h", groups=factor("both", levels=tlevels), 
                                name="#samples", cex.title=1, cex.axis=1,
                                font.axis=2, col.title=colM, col.axis=colM, 
                                background.title=colB, legend=TRUE)
        tlist <- c(tlist, both.track)
    }

    otrack <- Gviz::OverlayTrack(trackList = tlist, background.title=colB)
    ylim <- vapply(tlist, function(tr) range(S4Vectors::values(tr)), numeric(2))
    ylim <- extendrange(range(as.vector(ylim)))

    tracklist <- list(itrack, gtrack, otrack)

    # significant regions
    sig.regs <- subset(chr.regs, pvalue < pthresh)
    if(length(sig.regs))
    {
        atrack <- Gviz::AnnotationTrack(sig.regs, name="recur", cex.title=1, 
                                    cex.axis=1, font.axis=2, col.title=colM,
                                    col.axis=colM, background.title=colB)
        tracklist <- c(tracklist, atrack)
    }

    Gviz::plotTracks(tracklist, ylim=ylim)
}

#' OncoPrint plot for CNV regions
#'
#' Illustrates overlaps between CNV calls and genomic features across a 
#' sample population.
#'
#' @param calls Either a \code{\linkS4class{GRangesList}} or
#' \code{\linkS4class{RaggedExperiment}} storing the individual CNV calls for
#' each sample.
#' @param features A \code{\linkS4class{GRanges}} object containing
#' the genomic features of interest, typically genes. Feature names
#' are either expected as a meta-column \code{symbol} or as the \code{names}
#' of the object. 
#' @param multi.calls A function. Determines how to summarize the
#' CN state in a CNV region when there are multiple (potentially conflicting)
#' calls for one sample in that region. Defaults to \code{.largest}, which
#' assigns the CN state of the call that covers the largest part of the CNV
#' region tested. A user-defined function that is passed on to
#' \code{\link{qreduceAssay}} can also be provided for customized behavior.
#' @param top.features integer. Restricts the number of features for plotting to
#' features experiencing highest alteration frequency. Defaults to 25. 
#' Use \code{-1} to display all features.
#' @param top.samples integer. Restricts the number of samples for plotting to
#' samples experiencing highest alteration frequency. Defaults to 100.
#' Use \code{-1} to display all samples.
#' @param ... Additional arguments passed on to \code{ComplexHeatmap::oncoPrint}
#' @return None. Plots to a graphics device.
#'
#' @author Ludwig Geistlinger
#' @seealso \code{ComplexHeatmap::oncoPrint}
#' 
#' @examples
#'
#' # read in example CNV calls
#' data.dir <- system.file("extdata", package="CNVRanger")
#' call.file <- file.path(data.dir, "Silva16_PONE_CNV_calls.csv")
#' calls <- read.csv(call.file, as.is=TRUE)
#'
#' # store in a GRangesList
#' calls <- makeGRangesListFromDataFrame(calls, 
#'    split.field="NE_id", keep.extra.columns=TRUE)
#'
#' # three example genes
#' genes <- c(  "chr1:140368053-140522639:-", 
#'              "chr2:97843887-97988140:+",
#'              "chr2:135418586-135422028:-")
#' names(genes) <- c("ATP2C1", "MAP2", "ACTL8")
#' genes <- GRanges(genes)
#' 
#' # plot
#' cnvOncoPrint(calls, genes, top.samples = 25)
#'
#' @export
cnvOncoPrint <- function(calls, features, multi.calls=.largest,
                            top.features = 25, top.samples = 100, ...)
{
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
        stop(paste("Required package \'ComplexHeatmap\' not found.", 
                    "Use \'BiocManager::install(\"ComplexHeatmap\") to install it."))

    # summarize CNV calls in feature regions
    calls <- as(calls, "RaggedExperiment")
    qassay <- RaggedExperiment::qreduceAssay(calls, features, 
                                    simplifyReduce=multi.calls, background=2)

    if("symbol" %in% colnames(mcols(features))) rownames(qassay) <- features$symbol
    else rownames(qassay) <- names(features)

    qassay <- qassay[rownames(qassay) != "",]
    qassay <- qassay[!duplicated(rownames(qassay)),]
    
    # select most frequently altered features
    indr <- order(rowSums(qassay != 2), decreasing=TRUE)
    has.top.features <- top.features > 0 && top.features < nrow(qassay)
    if(has.top.features) indr <- indr[seq_len(top.features)]

    # select most frequently altered samples
    indc <- order(colSums(qassay != 2), decreasing=TRUE)
    has.top.samples <- top.samples > 0 && top.samples < ncol(qassay)
    if(has.top.samples) indc <- indc[seq_len(top.samples)]
    
    qassay <- qassay[indr, indc] 

    # recode genotypes
    qassay[qassay == 2] <- " "
    qassay[qassay == "0"] <- "HOMDEL"
    qassay[qassay == "1"] <- "HETDEL"
    qassay[qassay == "3"] <- "GAIN1"
    qassay[qassay == "4"] <- "GAIN2+"


    # assign colors to genotypes
    cb.red <- "#D55E00"
    cb.blue <- "#0072B2"
    cb.lightblue <- "#56B4E9"
    cb.orange <- "#E69F00"

    colors = c( "HOMDEL" = cb.blue, 
                "HETDEL" = cb.lightblue, 
                "GAIN1" = cb.orange, 
                "GAIN2+" = cb.red)

    # map functions to CNV types
    .mf <- function(couleur) 
    { 
        args <- alist(x =, y =, w =, h =)
        args <- as.pairlist(args)
        body <- substitute({
            grid::grid.rect(x, y, 
                            w - grid::unit(0.5, "mm"), 
                            h - grid::unit(0.5, "mm"), 
                            gp = grid::gpar(fill = z, col = NA))
        }, list(z = couleur))
        eval(call("function", args, body))
    }
    mutfuns <- lapply(colors, .mf)

    background <- function(x, y, w, h) 
        grid::grid.rect(x, y, 
                        w - grid::unit(0.5, "mm"), 
                        h - grid::unit(0.5, "mm"),
            gp = grid::gpar(fill = "#FFFFFF", col = "#FFFFFF"))
    mutfuns2 <- c(background = background, mutfuns)
   
    # plot
    heatmap_legend_param <- list(title = "CNV type", 
        at = c("HOMDEL", "HETDEL", "GAIN1", "GAIN2"), 
        labels = c("2-copy loss", 
                    "1-copy loss", 
                    "1-copy gain",
                    ">= 2-copy gain"))

    suppressMessages(
        ComplexHeatmap::oncoPrint(
            qassay, alter_fun = mutfuns2, col = colors, show_pct = FALSE,
            heatmap_legend_param = heatmap_legend_param, ...) 
    )
}


## (1) Density approach
.densityPopRanges <- function(grl, density)
{
    gr <- unlist(grl)
    cover <- GenomicRanges::reduce(gr)
    disjoint <- GenomicRanges::disjoin(gr)

    olaps <- GenomicRanges::findOverlaps(disjoint, cover)
    covered.hits <- S4Vectors::subjectHits(olaps)
    cover.support <- GenomicRanges::countOverlaps(cover, gr)[covered.hits]
    ppn <- GenomicRanges::countOverlaps(disjoint, gr) / cover.support
    
    pranges <- GenomicRanges::reduce(disjoint[ppn >= density])
    return(pranges)    
}

.classifyRegs <- function(regs, calls, type.thresh=0.1)
{
    if(!("state" %in% colnames(mcols(calls))))
        stop("Required column \'state\' storing integer copy number state not found")

    olaps <- GenomicRanges::findOverlaps(regs, calls)
    qh <- S4Vectors::queryHits(olaps)
    sh <- S4Vectors::subjectHits(olaps)
    
    # number of samples
    sampL <- split(names(calls)[sh], qh)
    .nrSamples <- function(x) length(unique(x))
    samples <- vapply(sampL, .nrSamples, numeric(1), USE.NAMES=FALSE)
    regs$freq <- samples 
    
    # type: gain, loss, both
    stateL <- split(calls$state[sh], qh)
    types <- vapply(stateL, .getType, character(1), 
                        type.thresh=type.thresh, USE.NAMES=FALSE)    
    regs$type <- types

    return(regs)
}

# @param states integer vector
# @param state.thresh numeric 
# @return character
.getType <- function(states, type.thresh=0.1)
{
    fract.amp <- mean(states > 2)
    fract.del <- mean(states < 2)    
    type <- c(fract.amp, fract.del) > type.thresh
    type <- c("gain", "loss")[type]
    if(length(type) > 1) type <- "both"
    return(type)
}   

## (2) RO approach
#
# @param grl GenomicRanges::GRangesList
# @param ro.tresh numeric
# @param multi.assign logical
# @param verbose logical
# @return GenomicRanges::GRanges
.roPopRanges <- function(grl, 
    ro.thresh=0.5, multi.assign=FALSE, verbose=FALSE)
{
    gr <- unlist(grl)
    S4Vectors::mcols(gr) <- NULL  

    # build initial clusters
    init.clusters <- GenomicRanges::reduce(gr)
    
    if(verbose) message(paste("TODO:", length(init.clusters)))

    # cluster within each initial cluster
    olaps <- GenomicRanges::findOverlaps(init.clusters, gr)
    qh <- S4Vectors::queryHits(olaps)
    sh <- S4Vectors::subjectHits(olaps)
    cl.per.iclust <- lapply(seq_along(init.clusters), 
        function(i)
        {
            if(verbose) message(i)
            # get calls of cluster
            ind <- sh[qh==i]
            ccalls <- gr[ind]
            if(length(ccalls) < 2) return(ccalls) 
            clusters <- .clusterCalls(ccalls, ro.thresh, multi.assign)
            if(is.list(clusters)) 
                clusters <- IRanges::extractList(ccalls, clusters) 
            clusters <- range(clusters)
            if(is(clusters, "GRangesList")) clusters <- sort(unlist(clusters))  
            return(clusters)
    })
    ro.ranges <- unname(unlist(GenomicRanges::GRangesList(cl.per.iclust)))
    return(ro.ranges)
}


# the clustering itself then goes sequentially through the identified RO hits, 
# touching each hit once, and checks whether this hit could be merged to 
# already existing clusters
#
# @param calls GenomicRanges::GRanges
# @param ro.thresh numeric
# @param multi.assign logical
# @return list of integer vectors
.clusterCalls <- function(calls, ro.thresh=0.5, multi.assign=FALSE)
{
    hits <- .getROHits(calls, ro.thresh)        
    
    # exit here if not 2 or more hits
    if(length(hits) < 2) return(calls)       
  
    # worst case: there are as many clusters as hits
    cid <- seq_along(hits)
    qh <- S4Vectors::queryHits(hits)   
    sh <- S4Vectors::subjectHits(hits)   
 
    # touch each hit once and check whether ... 
    # ... it could be merged to a previous cluster
    for(i in 2:length(hits))
    {
        # has this hit already been merged?
        if(cid[i] != i) next 
        curr.hit <- hits[i]
       
        # check each previous cluster
        for(j in seq_len(i-1))
        {
            # has this hit already been merged?
            if(cid[j] != j) next 
            
            # if not, check it
            prev.cluster <- hits[cid == j]
            mergeIndex <- .getMergeIndex(curr.hit, prev.cluster, hits)
            
            if(!is.null(mergeIndex))
            {
                if(length(mergeIndex) == 1) cid[i] <- j
                else cid[mergeIndex] <- j
                break 
            }
        }
    }
   
    # compile hit clusters 
    hit.clusters <- unname(S4Vectors::splitAsList(hits, cid))

    # extract call clusters
    call.clusters <- lapply(hit.clusters, 
        function(h) union(S4Vectors::queryHits(h), S4Vectors::subjectHits(h)))
    
    # can calls be assigned to more than one cluster?
    if(!multi.assign) call.clusters <- .pruneMultiAssign(call.clusters)
    
    return(call.clusters)
}

# given a set individual calls, returns overlaps (hits) between them 
# that satisfy the RO threshold
#
# @param calls GenomicRanges::GRanges
# @param ro.thresh numeric
# @return S4Vectors::Hits
.getROHits <- function(calls, ro.thresh=0.5)
{
    # calculate pairwise ro
    hits <- GenomicRanges::findOverlaps(calls, drop.self=TRUE, drop.redundant=TRUE)
    
    x <- calls[S4Vectors::queryHits(hits)]
    y <- calls[S4Vectors::subjectHits(hits)]
    pint <- GenomicRanges::pintersect(x, y)
    rovlp1 <- GenomicRanges::width(pint) / GenomicRanges::width(x)
    rovlp2 <- GenomicRanges::width(pint) / GenomicRanges::width(y)

    # keep only hits with ro > threshold
    ind <- rovlp1 > ro.thresh & rovlp2 > ro.thresh
    hits <- hits[ind]

    # exit here if not 2 or more hits
    if(length(hits) < 2) return(hits)       

    rovlp1 <- rovlp1[ind]
    rovlp2 <- rovlp2[ind]
    
    # order hits by RO 
    pmins <- pmin(rovlp1, rovlp2)
    ind <- order(pmins, decreasing=TRUE)
    
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)       
    hits <- S4Vectors::Hits(qh[ind], sh[ind], 
                S4Vectors::queryLength(hits), S4Vectors::subjectLength(hits))
    S4Vectors::mcols(hits)$RO1 <- rovlp1[ind]
    S4Vectors::mcols(hits)$RO2 <- rovlp2[ind]

    return(hits)
}

# decides whether a given hit can be merged to an already existing cluster
# mergeability requires that all cluster members satisfy the pairwise RO 
# threshold
#
# @param hit S4Vectors::Hits
# @param cluster S4Vectors::Hits
# @param hits S4Vectors::Hits
# @return integer vector
.getMergeIndex <- function(hit, cluster, hits)
{
    # (1) check whether query / S4Vectors::subject of hit is part of cluster
    curr.qh <- S4Vectors::queryHits(hit)
    curr.sh <- S4Vectors::subjectHits(hit)
       
    prev.qh <- S4Vectors::queryHits(cluster)
    prev.sh <- S4Vectors::subjectHits(cluster)
    prev.members <- union(prev.qh, prev.sh)
    is.part <- c(curr.qh, curr.sh) %in% prev.members

    # (2) can it be merged?    
    mergeIndex <- NULL

    # (2a) query *and* subject of hit are part of cluster
    if(all(is.part)) mergeIndex <- 1
    
    # (2b) query *or* subject of hit are part of cluster
    else if(any(is.part))
    {
        # check whether the call which is not part of the cluster
        # has sufficient RO with all others in the cluster 
        npart <- c(curr.qh, curr.sh)[!is.part]
        len <- length(prev.members)
        req.hits <- S4Vectors::Hits(   rep.int(npart, len), 
                            prev.members, 
                            S4Vectors::queryLength(hits), 
                            S4Vectors::subjectLength(hits))        
        is.gr <- S4Vectors::queryHits(req.hits) > S4Vectors::subjectHits(req.hits) 
        req.hits[is.gr] <- t(req.hits[is.gr])
        mergeIndex <- match(req.hits, hits)        
        if(any(is.na(mergeIndex))) mergeIndex <- NULL
    }
    return(mergeIndex)
}

# as the outlined procedure can assign a call to multiple clusters 
# (in the most basic case a call A that has sufficient RO with a call B and 
# a call C, but B and C do not have sufficient RO), this allows to optionally 
# strip away such multi-assignments
#
# @param clusters list of integer vectors
# @return list of integer vectors 
.pruneMultiAssign <- function(clusters)
{
    cid <- seq_along(clusters)
    times <- lengths(clusters)
    cid <- rep.int(cid, times)
    ind <- unlist(clusters)
    ndup <- !duplicated(ind)
    ind <- ind[ndup]
    cid <- cid[ndup]
       
    pruned.clusters <- split(ind, cid)
    return(pruned.clusters)
}

.estimateRecurrence <- function(regs, calls, mode=c("approx", "perm"),
    perm=1000, genome="hg19")
{
    mode <- match.arg(mode)
    calls <- unlist(calls)
    obs.scores <- lapply(seq_along(regs), 
                function(i) .scoreRegion(regs[i], calls))
    
    if(mode == "approx") pvals <- .approxP(obs.scores, regs$type)
    else pvals <- .permP(obs.scores, regs, calls, perm, genome)
    regs$pvalue <- pvals
    return(regs)
}

.approxP <- function(obs.scores, types)
{
    both.scores <- obs.scores[types == "both"]
    both.scores <- unlist(both.scores)
    gain.ind <- seq(1, length(both.scores) - 1, by=2)
    loss.ind <- seq(2, length(both.scores), by=2)
    gain.scores <- obs.scores[types == "gain"]
    gain.scores <- c(unlist(gain.scores), both.scores[gain.ind])
    loss.scores <-  obs.scores[types == "loss"]  
    loss.scores <- c(unlist(loss.scores), both.scores[loss.ind])

    gain.ecdf <- ecdf(gain.scores)
    loss.ecdf <- ecdf(loss.scores)
    
    .getP <- function(i)
    {
        s <- obs.scores[[i]]
        ty <- types[i]
        if(ty == "both")
        {     
            m <- which.max(s)
            p <- ifelse(m == 1, gain.ecdf(s), loss.ecdf(s))
        }
        else p <- ifelse(ty == "gain", gain.ecdf(s), loss.ecdf(s))
        return(1 - p)
    } 
    
    pvals <- vapply(seq_along(obs.scores), .getP, numeric(1))
    return(pvals)
}

.permP <- function(obs.scores, regs, calls, perm, genome)
{
    # genome specified?
    g <- GenomeInfoDb::genome(regs)
    ind <- !is.na(g)
    if(any(ind)) genome <- g[ind] 
    rregs <- regioneR::randomizeRegions(calls, genome=genome, mask=NA, 
                allow.overlaps=TRUE, per.chromosome=TRUE)
    # TODO:
    pvals <- vapply(seq_along(regs), function(i) i, numeric(1))
    return(pvals)
}

# @param r GenomicRanges::GRanges
# @param calls GenomicRanges::GRanges
# @param fract.thresh numeric
# @return numeric
.scoreRegion <- function(r, calls, fract.thresh=0.1)
{
    chr <- GenomicRanges::seqnames(r)
    chr <- as.character(chr)

    ind <- GenomicRanges::seqnames(calls) == chr
    rcalls <- GenomicRanges::restrict(calls[ind],
                        start=GenomicRanges::start(r), 
                        end=GenomicRanges::end(r))

    type <- .getType(rcalls$state, fract.thresh) 
    if(type != "both") score <- .getFreq(rcalls, type)
    else score <- vapply(c("gain", "loss"), 
                    function(type) .getFreq(rcalls, type), numeric(1))
    return(score) 
}

# @param calls GenomicRanges::GRanges
# @param type character
# @return numeric 
.getFreq <- function(calls, type=c("gain", "loss"))
{   
    type <- match.arg(type)
    is.amp <- type == "gain"

    chr <- GenomicRanges::seqnames(calls)[1]
    chr <- as.character(chr)
    
    ind <- if(is.amp) calls$state > 2 else calls$state < 2
    scalls <- calls[ind]
    
    w <- if(is.amp) scalls$state - 2 else 2 - scalls$state
    x <- GenomicRanges::coverage(scalls, weight=w)[[chr]]
    cs <- cumsum(S4Vectors::runLength(x))

    len <- length(cs)
    if(len == 1) return(0)
 
    grid <- seq_len(len - 1)
    starts <-  cs[grid] + 1 
    ends <- cs[grid + 1]
    cov <- S4Vectors::runValue(x)[grid + 1]
    
    max.score <- max(cov)
    return(max.score)
}
