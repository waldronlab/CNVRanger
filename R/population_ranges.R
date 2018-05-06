############################################################
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
#' @param density Numeric. 
#'
#' @return A \code{\linkS4class{GRanges}} object containing the summarized
#' genomic ranges. 
#'
#' @references 
#' Kim JH et al. (2012) CNVRuler: a copy number variation-based case-control
#' association analysis tool. Bioinformatics, 28(13):1790-2.
#'
#' @author Martin Morgan
#' @seealso
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
    cover <- reduce(gr)
    disjoint <- disjoin(gr)

    dj_covered_hits <- subjectHits(findOverlaps(disjoint, cover))
    cover_support <- countOverlaps(cover, gr)[dj_covered_hits]
    ppn <- countOverlaps(disjoint, gr) / cover_support
    reduce(disjoint[ppn >= density])
}


#
# Reciprocal overlap procedure
#

# reciprocal overlap (RO) approach (e.g. Conrad et al, Nature, 2010)
#
# reciprocal overlap of 0.51 between two CNV calls A and B:
#
# requires that B overlap at least 51% of A, 
#   *and* that A also overlaps at least 51% of B

# Approach:
#
# At the top level of the hierarchy, all contiguous bases overlapping at least 
# 1bp of a CNV call are merged into a “CNV region” (CNVR). 
# Within each CNVR we further define CNVs with the following algorithm:
#
# 1.   Calculate reciprocal overlap (RO) between all remaining calls.
#
# 2.   Identify pair of calls with greatest RO. If RO > threshold, merge and 
#       create a new CNV. If not, exit.
#
# 3.    Continue adding unclustered calls to the CNV, in order of best overlap. 
#       In order to add a call, the new call must have > threshold to all calls 
#       within CNV to be added. When no additional calls may be added, move to 
#       next step.
#
# 4.   If calls remain, return to 1. Otherwise exit.
#
# A typical reciprocal overlap threshold value is 0.5 for constructing CNVs

# this clustering procedure is then invoked on each initial cluster, which are, 
# according to the procedure outlined by Conrad et al., constructed by merging 
# all calls with >= 1 bp overlap.
populationRangesRO <- function(grl, ro.thresh=0.5, multi.assign=FALSE)
{
    gr <- unlist(grl)

    # build initial clusters
    init.clusters <- reduce(gr)
    message(paste("TODO:", length(init.clusters)))

    # cluster within each initial cluster
    cl.per.iclust <- sapply(init.clusters, 
        function(ic)
        {
            message(which(init.clusters == ic))
            # get calls of cluster
            ccalls <- subsetByOverlaps(gr, ic)
            clusters <- .clusterCalls(ccalls, ro.thresh, multi.assign)
            clusters <- range(extractList(ccalls, clusters))
            clusters <- sort(unlist(clusters))  
    })
    ro.ranges <- unname(unlist(GRangesList(cl.per.iclust)))
    return(ro.ranges)
}


# given a set individual calls, returns overlaps (hits) between them 
# that satisfy the RO threshold
.getROHits <- function(calls, ro.thresh=0.5)
{
    # calculate pairwise ro
    hits <- findOverlaps(calls, drop.self=TRUE, drop.redundant=TRUE)
    
    x <- calls[queryHits(hits)]
    y <- calls[subjectHits(hits)]
    pint <- pintersect(x, y)
    rovlp1 <- width(pint) / width(x)
    rovlp2 <- width(pint) / width(y)

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
    
    qh <- queryHits(hits)
    sh <- subjectHits(hits)       
    hits <- Hits(qh[ind], sh[ind], queryLength(hits), subjectLength(hits))
    mcols(hits)$RO1 <- rovlp1[ind]
    mcols(hits)$RO2 <- rovlp2[ind]

    return(hits)
}

# decides whether a given hit can be merged to an already existing cluster
# mergeability requires that all cluster members satisfy the pairwise RO 
# threshold
.getMergeIndex <- function(hit, cluster, hits)
{
    # (1) check whether query / subject of hit is part of cluster
    curr.qh <- queryHits(hit)
    curr.sh <- subjectHits(hit)
       
    prev.qh <- queryHits(cluster)
    prev.sh <- subjectHits(cluster)
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
        req.hits <- Hits(   rep(npart, length(prev.members)), 
                            prev.members, 
                            queryLength(hits), 
                            subjectLength(hits))
        is.gr <- queryHits(req.hits) > subjectHits(req.hits) 
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
    times <- sapply(clusters, length)
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
    hit.clusters <- unname(splitAsList(hits, cid))

    # extract call clusters
    call.clusters <- lapply(hit.clusters, 
        function(h) union(queryHits(h), subjectHits(h)))
    
    # can calls be assigned to more than one cluster?
    if(!multi.assign) call.clusters <- .pruneMultiAssign(call.clusters)
    
    return(call.clusters)
}

