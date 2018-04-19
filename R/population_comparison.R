##################
# Author: Vinicius Henrique da Silva
# Script description: Functions related to population comparison and evolutionary analysis
# Date: April 19, 2018
# Code: CNVEVOPACK001
###################

#' HELPER - replace the genofile at the gds by the CNV genotype at each probe
#' @param lo loop number
#' @param genofile loaded gds file
#' @param n is the object from \code{gdsfmt::add.gdsn}
#' @param ranges.gr is the CNVs of each sample

.replaceSNPtoCNV <- function(lo, genofile, n, all.sample, ranges.gr, map.cnv){
  
  sampleX <- all.sample[[lo]]
  g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id=sampleX))
  
  #2n 
  
  g[seq_len(length(g))] <- "1"
  
  # 0n 
  ranges.gr0n <- ranges.gr[which(values(ranges.gr)$sample == sampleX)]
  ranges.gr0n <- ranges.gr0n[which(values(ranges.gr0n)$type == "state1,cn=0")]
  map.cnv.0n <- subsetByOverlaps(map.cnv, ranges.gr0n)
  
  g[values(map.cnv.0n)$tag.snp] <- "0"
  
  # 1n
  ranges.gr1n <- ranges.gr[which(values(ranges.gr)$sample == sampleX)]
  ranges.gr1n <- ranges.gr1n[which(values(ranges.gr1n)$type == "state2,cn=1")]
  map.cnv.1n <- subsetByOverlaps(map.cnv, ranges.gr1n)
  
  g[values(map.cnv.1n)$tag.snp] <- "0"
  
  # 3n
  ranges.gr3n <- ranges.gr[which(values(ranges.gr)$sample == sampleX)]
  ranges.gr3n <- ranges.gr3n[which(values(ranges.gr3n)$type == "state5,cn=3")]
  map.cnv.3n <- subsetByOverlaps(map.cnv, ranges.gr3n)
  
  g[values(map.cnv.3n)$tag.snp] <- "2"
  
  # 4n
  ranges.gr4n <- ranges.gr[which(values(ranges.gr)$sample == sampleX)]
  ranges.gr4n <- ranges.gr4n[which(values(ranges.gr4n)$type == "state6,cn=4")]
  map.cnv.4n <- subsetByOverlaps(map.cnv, ranges.gr4n)
  
  g[values(map.cnv.4n)$tag.snp] <- "2"
  
  # 2n
  #g[c(-(values(map.cnv.0n)$tag.snp),-(values(map.cnv.1n)$tag.snp),-(values(map.cnv.3n)$tag.snp),-(values(map.cnv.4n)$tag.snp))] <- "1"
  
  g <- as.numeric(g)
  
  ### Replace the genotype in the gds file
  gdsfmt::write.gdsn(n, g, start=c(1,lo), count=c(length(g),1))
  
}


#' Construct the CNV gds as SNPs (i.e. 0, 1 or 2)
#' @param path.gds to the gds to be imported
#' @param final.name of the gds file after inclusion of CNV genotypes
#' @param ranges.gr is the CNVs of each sample
#' @param n.cor is the number of cores to be used - one is default
#' @param name String with a project code or name
#' @param folder Choose manually the project folder (i.e. path as the root folder)
#' @return List with paths placed to store the files produced by subsequent analysis 
#' @export

produceSNPlikeCNVgds <- function(path.gds, final.name, ranges.gr, n.cor=1, name, folder=NULL){
  
  all.paths <- .createFolderTree(name, folder)  
  
  file.copy(file.path(path.gds), file.path(all.paths[1], "CNVtoAdmixture.gds"), overwrite=TRUE)
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNVtoAdmixture.gds"), readonly = FALSE, allow.fork=TRUE))
  
  # Produce map 
  map <- as.data.frame(cbind(snp.id = read.gdsn(index.gdsn(genofile, "snp.id")), chr = read.gdsn(index.gdsn(genofile, "snp.chromosome")), 
                             position = read.gdsn(index.gdsn(genofile, "snp.position")), 
                             alleles = read.gdsn(index.gdsn(genofile, "snp.allele")),
                             probes = read.gdsn(index.gdsn(genofile, "snp.rs.id"))
  ))
  
  map.gr <- makeGRangesFromDataFrame(map, seqnames.field = "chr", start.field = "position",
                                     end.field = "position", keep.extra.columns = TRUE)
  map.cnv <- subsetByOverlaps(map.gr, ranges.gr)
  values(map.cnv)$tag.snp <- seq_len(length(map.cnv))
  
  gdsfmt::closefn.gds(genofile)
  
  # Select only SNPs overlapping CNVs
  GWASTools::gdsSubset(file.path(all.paths[1], "CNVtoAdmixture.gds"), file.path(all.paths[1], "CNVtoAdmixtureM.gds"), 
            snp.include=as.character(values(map.cnv)$snp.id))
  file.rename(file.path(all.paths[1], "CNVtoAdmixtureM.gds"), file.path(all.paths[1], "CNVtoAdmixture.gds"))
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNVtoAdmixture.gds"), readonly = FALSE, allow.fork=TRUE)) 
  
  genoFile <- (g <- gdsfmt::read.gdsn(index.gdsn(genofile, "genotype"))) 
  
  n <- gdsfmt::add.gdsn(genofile, "genotype", genoFile, replace=TRUE)
  gdsfmt::read.gdsn(n)
  
  all.snp <- gdsfmt::read.gdsn(index.gdsn(genofile, "snp.id"))
  all.sample <- gdsfmt::read.gdsn(index.gdsn(genofile, "sample.id"))
  
  if(!all(values(ranges.gr)$sample %in% all.sample)){
    stop("CNV ranges and GDS does not have the same sample set")
  }
  
  if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
    multicoreParam <-  BiocParallel::MulticoreParam(workers = n.cor)
    BiocParallel::bplapply(1:length(all.sample), .replaceSNPtoCNV, BPPARAM = multicoreParam, genofile=genofile, n=n,
                           all.sample=all.sample, ranges.gr=ranges.gr, map.cnv=map.cnv)
  }
  
  if(rappdirs:::get_os() == "win"){
    param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
    suppressMessages(BiocParallel::bplapply(1:length(all.sample), .replaceSNPtoCNV, BPPARAM = multicoreParam, genofile=genofile, n=n,
                                            all.sample=all.sample, ranges.gr=ranges.gr, map.cnv=map.cnv))
  }
  
  gdsfmt::closefn.gds(genofile)
  
  file.rename(file.path(all.paths[1], "CNVtoAdmixture.gds"), file.path(all.paths[1], final.name))
  
  return(all.paths)
  
}

