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
#' CNVRuler procedure that trims region margins based on regional density
#'
#' In CNV analysis, it is often of interest to summarize invidual calls across 
#' the population, (i.e. to define CNV regions), for subsequent association 
#' analysis with e.g. phenotype data.
#'
#' In the simplest case, this just merges overlapping individual calls into 
#' summarized regions.
#' However, this typically inflates CNVR size and trimming low-density areas
#' (usually <10\% of the total contributing individual calls within a summarized 
#' region) is advisable.
#' 
#' An illustration of the concept can be found here: 
#' https://www.ncbi.nlm.nih.gov/pubmed/22539667 (Figure 1)
#'
#' @param grl A \code{\linkS4class{GRangesList}}.
#' @param density Numeric. Defaults to 0.1. 
#'
#' @return A \code{\linkS4class{GRanges}} object containing the summarized
#' genomic ranges. 
#'
#' @references 
#' Kim JH et al. (2012) CNVRuler: a copy number variation-based case-control
#' association analysis tool. Bioinformatics, 28(13):1790-2.
#'
#' @author Martin Morgan
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
#' @export
populationRanges <- function(grl, density=0.1)
{
    gr <- unlist(grl)
    cover <- GenomicRanges::reduce(gr)
    disjoint <- GenomicRanges::disjoin(gr)

    dj_covered_hits <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(disjoint, cover))
    cover_support <- GenomicRanges::countOverlaps(cover, gr)[dj_covered_hits]
    ppn <- GenomicRanges::countOverlaps(disjoint, gr) / cover_support
    GenomicRanges::reduce(disjoint[ppn >= density])
}


#' Summarizing genomic ranges based on reciprocal overlap across a population
#'
#' Reciprocal overlap (RO) approach (e.g. Conrad et al., Nature, 2010)
#'
#' Reciprocal overlap of 0.51 between two genomic regions A and B:
#'
#' requires that B overlap at least 51% of A, 
#'   *and* that A also overlaps at least 51% of B
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
#' }
#'
#' @param grl A \code{\linkS4class{GRangesList}}.
#' @param ro.thresh Numeric. Threshold for reciprocal overlap required for merging
#' two overlapping regions. Defaults to 0.5.
#' @param multi.assign Logical. Allow regions to be assigned to several region clusters?
#' Defaults to FALSE.
#' 
#' @return A \code{\linkS4class{GRanges}} object containing the summarized
#' genomic ranges. 
#'
#' @references 
#' Conrad et al. (2010) Origins and functional impact of copy number variation 
#' in the human genome. Nature, 464(7289):704-12.
#'
#' @author Ludwig Geistlinger
#' @seealso \code{\link{subsetByOverlaps}}
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
#' populationRangesRO(grl, ro.thresh=0.5)
#'
#' @export
populationRangesRO <- function(grl, ro.thresh=0.5, multi.assign=FALSE)
{
    gr <- unlist(grl)

    # build initial clusters
    init.clusters <- GenomicRanges::reduce(gr)
    message(paste("TODO:", length(init.clusters)))

    # cluster within each initial cluster
    cl.per.iclust <- lapply(seq_along(init.clusters), 
        function(i)
        {
            message(i)
            ic <- init.clusters[i]
            # get calls of cluster
            ccalls <- IRanges::subsetByOverlaps(gr, ic)
            clusters <- .clusterCalls(ccalls, ro.thresh, multi.assign)
            clusters <- range(IRanges::extractList(ccalls, clusters))
            clusters <- sort(unlist(clusters))  
    })
    ro.ranges <- unname(unlist(GenomicRanges::GRangesList(cl.per.iclust)))
    return(ro.ranges)
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
    hits <- S4Vectors::Hits(qh[ind], sh[ind], S4Vectors::queryLength(hits), S4Vectors::subjectLength(hits))
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
        req.hits <- S4Vectors::Hits(   rep(npart, length(prev.members)), 
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
    cid <- rep(cid, times)
    ind <- unlist(clusters)
    ndup <- !duplicated(ind)
    ind <- ind[ndup]
    cid <- cid[ndup]
       
    pruned.clusters <- split(ind, cid)
    return(pruned.clusters)
}

# the clustering itself then goes sequentially through the identified RO hits, 
# touching each hit once, and checks whether this hit could be merged to 
# already existing clusters
.clusterCalls <- function(calls, ro.thresh=0.5, multi.assign=FALSE)
{
    hits <- .getROHits(calls, ro.thresh)        
    
    # exit here if not 2 or more hits
    if(length(hits) < 2) return(hits)       
  
    # worst case: there as many clusters as hits
    cid <- seq_along(hits)
    
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

