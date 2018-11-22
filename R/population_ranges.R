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

#' Summarizing genomic ranges across a population
#'
#' In CNV analysis, it is often of interest to summarize invidual calls across 
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
#' 1bp of indvidual calls are merged into one region. 
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
#' }}
#'
#' @param grl A \code{\linkS4class{GRangesList}}.
#' @param mode Character. Should population ranges be computed based on regional
#' density ("density") or reciprocal overlap ("RO"). See Details.
#' @param density Numeric. Defaults to 0.1. 
#' @param ro.thresh Numeric. Threshold for reciprocal overlap required for merging
#' two overlapping regions. Defaults to 0.5.
#' @param multi.assign Logical. Allow regions to be assigned to several region clusters?
#' Defaults to FALSE.
#' @param verbose Logical. Report progress messages? Defaults to FALSE.
#'
#' @return A \code{\linkS4class{GRanges}} object containing the summarized
#' genomic ranges. 
#'
#' @references 
#' Kim et al. (2012) CNVRuler: a copy number variation-based case-control
#' association analysis tool. Bioinformatics, 28(13):1790-2.
#'
#' Conrad et al. (2010) Origins and functional impact of copy number variation 
#' in the human genome. Nature, 464(7289):704-12.
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
#' populationRanges(grl, density=0.1)
#'
#' # density = 0 merges all overlapping regions, 
#' # equivalent to: reduce(unlist(grl))
#' populationRanges(grl, density=0) 
#'
#' # density = 1 disjoins all overlapping regions, 
#' # equivalent to: disjoin(unlist(grl))
#' populationRanges(grl, density=1)
#'
#' # RO procedure
#' populationRanges(grl, mode="RO", ro.thresh=0.5)
#'
#' @export
populationRanges <- function(grl, mode=c("density", "RO"), 
    density=0.1, ro.thresh=0.5, multi.assign=FALSE, verbose=FALSE)
{
    mode <- match.arg(mode)
    if(mode == "density") .densityPopRanges(grl, density)
    else .roPopRanges(grl, ro.thresh, multi.assign, verbose)
    #TODO: classify regs
    #TODO: filter for freq and size
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
    pranges <- pranges[GenomicRanges::width(pranges) > 1]
    return(pranges)    
}

.classifyRegs <- function(regs, calls, state.thresh=0.1)
{
    olaps <- GenomicRanges::findOverlaps(regs, calls)
    qh <- S4Vectors::queryHits(olaps)
    sh <- S4Vectors::subjectHits(olaps)
    
    # number of samples
    sampL <- split(names(calls)[sh], qh)
    .nrSamples <- function(x) length(unique(x))
    samples <- vapply(sampL, .nrSamples, numeric(1), USE.NAMES=FALSE)
    regs$freq <- samples 
    
    # type: gain, loss, both
    stateL <- split(calls$State[sh], qh)
    types <- vapply(stateL, .getType, character(1), 
                        state.thresh=state.thresh, USE.NAMES=FALSE)    
    regs$type <- types

    return(regs)
}

.getType <- function(states, state.thresh=0.1)
{
    fract.amp <- mean(states > 2)
    fract.del <- mean(states < 2)    
    type <- c(fract.amp, fract.del) > state.thresh
    type <- c("gain", "loss")[type]
    if(length(type) > 1) type <- "both"
    return(type)
}   

## (2) RO approach
.roPopRanges <- function(grl, 
    ro.thresh=0.5, multi.assign=FALSE, verbose=FALSE)
{
    gr <- unlist(grl)
    S4Vectors::mcols(gr) <- NULL  

    # build initial clusters
    init.clusters <- GenomicRanges::reduce(gr)
    # check 1269
    # init.clusters <- init.clusters[2500:2652]
    
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
    hit.clusters <- unname(IRanges::splitAsList(hits, cid))

    # extract call clusters
    call.clusters <- lapply(hit.clusters, 
        function(h) union(S4Vectors::queryHits(h), S4Vectors::subjectHits(h)))
    
    # can calls be assigned to more than one cluster?
    if(!multi.assign) call.clusters <- .pruneMultiAssign(call.clusters)
    
    return(call.clusters)
}

# given a set individual calls, returns overlaps (hits) between them 
# that satisfy the RO threshold
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


## (3) Recurrence approach
.recurrentPopRanges <- function(regs, calls, mode=c("approx", "perm"),
    perm=1000, genome="hg19")
{
    mode <- match.arg(mode)
    calls <- unlist(calls)
    obs.scores <- lapply(seq_along(regs), 
                function(i) .scoreRegion(regs[i], calls))
    

    # genome specified?
    g <- GenomeInfoDb::genome(regs)
    ind <- !is.na(g)
    if(any(ind)) genome <- g[ind] 
    
    
    rregs <- regioneR::randomizeRegions(calls, genome=genome, mask=NA, 
                    allow.overlaps=TRUE, per.chromosome=TRUE)
    pvals <- vapply(seq_along(regs)) 
}

.scoreRegion <- function(r, calls, fract.thresh=0.1)
{
    chr <- GenomicRanges::seqnames(r)
    chr <- as.character(chr)

    ind <- GenomicRanges::seqnames(calls) == chr
    rcalls <- GenomicRanges::restrict(calls[ind],
                        start=GenomicRanges::start(r), 
                        end=GenomicRanges::end(r))

    type <- .getType(rcalls$State, fract.thresh) 
    if(type != "both") score <- .getFreq(rcalls, type)
    else score <- vapply(c("gain", "loss"), 
                    function(type) .getFreq(rcalls, type), numeric(1))
    return(score) 
}

.getFreq <- function(calls, type=c("gain", "loss"))
{   
    type <- match.arg(type)
    is.amp <- type == "gain"

    chr <- GenomicRanges::seqnames(calls)[1]
    chr <- as.character(chr)
    
    ind <- if(is.amp) calls$State > 2 else calls$State < 2
    scalls <- calls[ind]
    
    w <- scalls$State + ifelse(is.amp, -2, 2)
    x <- GenomicRanges::coverage(scalls, weight=w)[[chr]]
    cs <- cumsum(S4Vectors::runLength(x))

    len <- length(cs)
    grid <- seq_len(len-1)
    starts <-  cs[grid] + 1 
    ends <- cs[grid + 1]
    cov <- S4Vectors::runValue(x)[grid + 1]
    
    max.score <- max(cov)
    return(max.score)
}
