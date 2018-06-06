##################
# Author: Vinicius Henrique da Silva
# Script description: Functions related to population comparison and evolutionary analysis
# Date: April 19, 2018
# Code: CNVEVOPACK001
###################

#' HELPER - Formula to infer Vst
#' 
#' 

.vstForm <- function(lo, value, pops.names){
  message(paste0("Vst ", lo, " of ", nrow(value)))
  value <- value[lo,]
  pops.unique  <- unique(pops.names)
  names(value) <- pops.names
  popc <- as.numeric(value)
  popc <- na.omit(popc)
  popa <- as.numeric(value[which(names(value)==pops.unique[1])])
  popa <- na.omit(popa)
  popb <- as.numeric(value[which(names(value)==pops.unique[2])])
  popb <- na.omit(popb)
  
  VarTotal=var(popc)
  VarPopa=var(popa)
  VarPopb=var(popb)
  VarWithin=VarPopa*length(popa)+VarPopb*length(popb)
  AvVarWithin=VarWithin/length(popc)
  VarDiff=VarTotal-AvVarWithin
  vst=VarDiff/VarTotal
  
  return(vst)
  
}

#' Vst
#' Estimate the Vst among two populations from a gds previouly loaded with LRR (with \link{importLRR_BAF} function)
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree 
#' @param pops Two populations to be compared 
#' @param chr.code.name A data-frame with the integer name in the first column and the original name for each chromosome  
#' @return Vst for all probes with CNV genotypes listed in the gds

calcVst <- function(all.paths, pops, chr.code.name=NULL, n.cor=1){
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  #(genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNVLRRBAFAllSNPs.gds"), allow.fork=TRUE, readonly=FALSE))
  
  pops.names <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Pops.names")))  
  
  if(all(pops %in% pops.names)){
    
    value <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "LRR")))
    value <- value[,which(pops.names %in% pops)]
    
    pops.names <- pops.names[pops.names %in% pops]
    
    
    if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
      multicoreParam <-  BiocParallel::MulticoreParam(workers = n.cor)
      vst.result <- BiocParallel::bplapply(1:nrow(value), .vstForm, BPPARAM = multicoreParam, value=value, pops.names=pops.names)
      
    }
    
    if(rappdirs:::get_os() == "win"){
      param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
      vst.result <- BiocParallel::bplapply(1:nrow(value), .vstForm, BPPARAM = param, value=value, pops.names=pops.names)
    }
    
    #vst.result <- apply(value, 1, FUN=.vstForm, pops.names=pops.names)
    vst.result <- unlist(vst.result)
    vst.result[is.na(vst.result)] <- 0
    
  }else{
    stop("Population names not listed in the gds file")
  }
  
  
  map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")), 
                             chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                             position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")), 
                             probes = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
  ))
  
  map$vst.result <- vst.result
  map$NameProbe <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")))
  
  if(!is.null(index.gdsn(genofile, "Chr.names", silent=TRUE))){
    
    chr.code.name < gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Chr.names"))
    
    for(lopN in seq_len(nrow(chr.code.name))){
      map$chr <- gsub(chr.code.name[lopN,1],chr.code.name[lopN,2], map$chr)
    }
    
    
  }
  
  map.vst <- makeGRangesFromDataFrame(map, seqnames.field = "chr", start.field = "position", end.field = "position",
                                      keep.extra.columns = TRUE)
  
  values(map.vst)$Populations <- paste(unlist(pops), sep = "", collapse ="-")
  
  SNPRelate::snpgdsClose(genofile)
  
  return(map.vst)
  
  
}


#' HELPER - PCA
#' 
#' 

PCAcnv <-  function(all.paths, segs.pvalue.gr, n.cor=1, pops.names=NULL, method.ld="dprime"){
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  
  ### change nodes to perform the PCA analysis
  genoBack <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "genotype")))  
  CNVgeno <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotypeSNPlike")))  
  
  n <- gdsfmt::add.gdsn(genofile, "genotype", CNVgeno, replace=TRUE)
  gdsfmt::read.gdsn(n)
  
  ### Selectinformative  probes with CNV genotype based on LD
  #snpset <- SNPRelate::snpgdsLDpruning(genofile, method=method.ld, remove.monosnp = FALSE)
  map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")), 
                             chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                             position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")), 
                             probes = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
  ))
  
  snp.id <- subset(map, probes %in% values(segs.pvalue.gr)$NameProbe)$snp.id
  
  
  ### Run PCA with all samples
  pca <- SNPRelate::snpgdsPCA(genofile, num.thread=n.cor, autosome.only=F, snp.id=snp.id)
  
  n <- gdsfmt::add.gdsn(genofile, "genotype", genoBack, replace=TRUE)
  gdsfmt::read.gdsn(n)
  
  # variance proportion (%)
  pc.percent <- pca$varprop*100
  
  # make a data.frame
  tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1],    # the first eigenvector
                    EV2 = pca$eigenvect[,2],    # the second eigenvector
                    stringsAsFactors = FALSE)
  
  
  if(!is.null(pops.names)){
    pops.names <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Pops.names")))  
    tab$pop <- pops.names
    
    pdf(file.path(all.paths[3], "PCA_analysis.pdf"))
    p <- ggplot2::ggplot(tab, ggplot2::aes(EV1, EV2))
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = factor(pop)), size = 4, alpha = 0.5, shape="o") + 
      ggplot2::guides(fill=ggplot2::guide_legend(title="Downstream inversion boundary"))
    print(p)
    dev.off()
  }else{
    
    pdf(file.path(all.paths[3], "PCA_analysis.pdf"))
    p <- ggplot2::ggplot(tab, ggplot2::aes(EV1, EV2))
    p <- p + ggplot2::geom_point(ggplot2::aes(size = 4, alpha = 0.5, shape="o"))
    print(p)
    dev.off()
    
  }
  
  SNPRelate::snpgdsClose(genofile)
}


#' HELPER - Admixture - TODO
#' 
#' 


