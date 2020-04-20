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
#' @param n.cor Number of cores to be used
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



#' PCA with CNVs
#' 
#' 
#'
#' @param all.paths Object returned from \code{CreateFolderTree} function with 
#' the working folder tree 
#' @param segs.pvalue.gr Result from \code{\link{cnvGWAS}} 
#' @param n.cor Number of cores
#' @param pops.names Name of the populations
#' @param method.ld LD method, "r2" or "dprime"
#' @param attach.pca Create a node in the CNV.gds containing the PCA analysis 
#' @author Vinicius Henrique da Silva <vinicius.dasilva@wur.nl>
#' @examples 
#' 
#' @export
#'  


PCAcnv <-  function(all.paths, segs.pvalue.gr, n.cor=1, pops.names=NULL, 
                    method.ld="dprime", attach.pca=TRUE, use.seg=TRUE){
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
  if(use.seg){
  pca <- SNPRelate::snpgdsPCA(genofile, num.thread=n.cor,
                              autosome.only=FALSE, snp.id=snp.id)
  }else{
    pca <- SNPRelate::snpgdsPCA(genofile, num.thread=n.cor,
                                autosome.only=FALSE)
  }
  
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
    
    if(attach.pca){
      gdsfmt::add.gdsn(genofile, name = "PCA", val = tab, replace = TRUE) 
    }
    
    pdf(file.path(all.paths[3], "PCA_analysis.pdf"))
    p <- ggplot2::ggplot(tab, ggplot2::aes(EV1, EV2))
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = factor(pop)), size = 4, alpha = 0.5, shape="o") + 
      ggplot2::guides(fill=ggplot2::guide_legend(title="Downstream inversion boundary"))
    print(p)
    dev.off()
  }else{
    
    if(attach.pca){
      gdsfmt::add.gdsn(genofile, name = "PCA", val = tab, replace = TRUE) 
    }
    
    pdf(file.path(all.paths[3], "PCA_analysis.pdf"))
    p <- ggplot2::ggplot(tab, ggplot2::aes(EV1, EV2))
    p <- p + ggplot2::geom_point(ggplot2::aes(size = 4, alpha = 0.5, shape="o"))
    print(p)
    dev.off()
    
  }
  
  SNPRelate::snpgdsClose(genofile)
  return(p)
}


#' HELPER - Admixture - TODO
#' 
#' 


#' Compare CNV segment frequencies among populations and draw Fst for the CNVgenotypeSNPlike node
#' 
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree
#' @param segs.pvalue.gr Object returned from \code{\link{setupCnvGWAS}} function
#' @param pops Two populations to be compared
#' @param p.value.min Minimum p-value to consider normal distribution  
#' @param n.cor Number of cores to be used
#' 
#' @export


compareCNVsegFreq <- function(all.paths, segs.pvalue.gr, pops, p.value.min=0.05, n.cor=1){
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  
  map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")), 
                             chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                             position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")), 
                             probes = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
  ))
  
  snp.id <- subset(map, probes %in% values(segs.pvalue.gr)$NameProbe)$snp.id
  snp.id <- as.numeric(as.character(snp.id))
  
  pops.names <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Pops.names")))  
  
  genos <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotypeSNPlike"))) 
  genos <- genos[snp.id, ] ## Subset the genotypes
  
  genos.pop1 <- genos[,which(pops.names== pops[1])] #pop1
  index.pop1 <- pops.names[which(pops.names== pops[1])]
  genos.pop2 <- genos[,which(pops.names== pops[2])] #pop2
  index.pop2 <- pops.names[which(pops.names== pops[2])]
  index.all <- c(index.pop1,index.pop2)
  
  SnpsAll <- cbind(genos.pop1, genos.pop2)
  SnpsAll <- SnpsAll+1
  SnpsAll <- new("SnpMatrix", t(SnpsAll))
  SnpsAll
  fstValues <- snpStats::Fst(SnpsAll, index.all, pairwise=FALSE)
  fst.map <- cbind("snp.id"=snp.id,"fst"=fstValues$Fst)
  fst.map <- as.data.frame(fst.map)
  fst.map <- merge(fst.map, map, by="snp.id", all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  fst.map.gr <- makeGRangesFromDataFrame(fst.map, start.field = "position", end.field = "position", keep.extra.columns = TRUE)
  
  node <- gdsfmt::index.gdsn(genofile, "Chr.names", silent=TRUE)
  
  if(!is.null(chr.code.name)){
    chr.code.name <- gdsfmt::read.gdsn(node)
    for(lopN in seq_len(nrow(chr.code.name))){
      seqlevels(fst.map.gr)[which(match(seqlevels(fst.map.gr), chr.code.name[lopN,1], nomatch=0) == 1)] <- as.character(chr.code.name[lopN,2])
    }
  }
  
  #### Check frequency distribution of the CNV segments - Normally distributed? TODO - Use general comparsion
  freqs <- as.numeric(values(segs.pvalue.gr)$Frequency)
  res.norm <- nortest::ad.test(freqs)
  
  avoid <- TRUE
  if(avoid==FALSE){
    if(res.norm$p.value<=p.value.min){
      
      #### Compare frequencies in parallel
      if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
        multicoreParam <-  BiocParallel::MulticoreParam(workers = n.cor)
        freq.result <- BiocParallel::bplapply(1:nrow(value), .compCNfreqs, BPPARAM = multicoreParam, 
                                              genofile=genofile,snp.id=snp.id, 
                                              pops.names=pops.names, pops=pops)
        #value=value, pops.names=pops.names)
        
      }
      
      if(rappdirs:::get_os() == "win"){
        param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
        freq.result <- BiocParallel::bplapply(1:nrow(value), .compCNfreqs, BPPARAM = param, 
                                              genofile=genofile,snp.id=snp.id, 
                                              pops.names=pops.names, pops=pops)
        #value=value, pops.names=pops.names)
      }
      
      
    }else(  
      stop("Not normally distributed") ## Include non-parametric alternative
    )
  } ## INCOMPLETE TODO - put methods from Maulik paper here
  
  SNPRelate::snpgdsClose(genofile)
  
  return(fst.map.gr)
  
}

### Tests

test.population <- function(x){
  ### Get Fst values and check genes around
  fst.map.gr <- compareCNVsegFreq(all.paths, segs.pvalue.gr, pops)
  values(fst.map.gr[is.na(values(fst.map.gr)$fst )])$fst <- 0
  
  ### Define threshould to the Fst
  quant.c <- quantile(values(fst.map.gr)$fst, probs = seq(0, 1, by=0.05))
  quant.c <- quant.c[(length(quant.c)-1)]
  quant.c
  
  fst.map.gr <- fst.map.gr[(elementMetadata(fst.map.gr)[, "fst"] >= as.numeric(quant.c))] ## Select by minimum fst value
  #fst.map.gr <- fst.map.gr[(elementMetadata(fst.map.gr)[, "fst"] >= 0.005)] ## Select by minimum fst value
  fst.map.gr
  
  ### Run PCA with higher fst values
  segs.pvalue.gr.fst <- subsetByOverlaps(segs.pvalue.gr, fst.map.gr)
  PCAcnv(all.paths, n.cor=1, pops.names=c("NL", "UK"), segs.pvalue.gr=segs.pvalue.gr.fst)
  .uploadAllFigures(all.paths, "CNVBreedingTime")
  
  
  load(file="/home/NIOO/viniciusd/great_tit/cnv/cnvrs/genes11NCBIGFFGr.rda")
  seqlevels(genes11Gr) <- gsub("chr", "", seqlevels(genes11Gr))
  subsetByOverlaps(genes11Gr, fst.map.gr)
  seqnames(genes11Gr)
  
  #' HELPER - compare CNV segment frequencies - t.test - INCOMPLETE - Implement strategy from Maulik's paper
  
  .compCNfreqs <- function(lo, genofile,snp.id, pops.names, pops){ ## INCOMPLETE
    genos <- .snpgdsGetGenoCNV(genofile, snp.id[lo])
    genos.pop1 <- genos[which(pops.names== pops[1])] #pop1
    genos.pop2 <- genos[which(pops.names== pops[2])] #pop2
    #t.test(y~x)   
    genoBack <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "genotype")))  
    CNVgeno <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotypeSNPlike"))) 
    n <- gdsfmt::add.gdsn(genofile, "genotype", CNVgeno, replace=TRUE)
    gdsfmt::read.gdsn(n)
    SnpsAll <- cbind(snps_inv, snps_noinv)
    SnpsAll <- SnpsAll+1
    SnpsAll1 <- SnpsAll
    SnpsAll <- new("SnpMatrix", t(SnpsAll))
    fstValues <- snpStats::Fst(genofile, num.thread=n.cor, autosome.only=F, snp.id=snp.id)
    n <- gdsfmt::add.gdsn(genofile, "genotype", genoBack, replace=TRUE)
    gdsfmt::read.gdsn(n)
  }
  
  
  #' HELPER - Gene enrichement
  #' 
  #' 
  nv
  ### Extend the fst genomic window
  start(fst.map.gr)<- start(fst.map.gr)-10000
  end(fst.map.gr) <- end(fst.map.gr)+10000
  genesEn <- values(genes11Gr)$Name ###### Extract overlapped gene names
  ids <- bitr(genesEn, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") ## human
  genes <- as.character(ids[,2])
  GTuniverseEgo <- as.character(ids[,2])
  eg2np <- bitr_kegg(genes, fromType='kegg', toType='uniprot', organism='hsa') ## human
  GTuniversekk <- as.character(eg2np[,2]) ## All proteins liked with GT genes
  
  
  genesEn <- values(subsetByOverlaps(genes11Gr, fst.map.gr))$Name ###### Extract overlapped gene names
  
  ids <- bitr(genesEn, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") ## human
  genes <- as.character(ids[,2])
  genesego <- as.character(ids[,2])
  eg2np <- bitr_kegg(genes, fromType='kegg', toType='uniprot', organism='hsa') ## human
  genes <- as.character(eg2np[,2])
  genes
  
  ## Erich KEGG pathways
  kk <- enrichKEGG(gene   = as.vector(genes),
                   universe = GTuniversekk,
                   organism     = "hsa", ## human
                   qvalueCutoff = 1,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   keyType ="uniprot",
                   minGSSize = 5,
                   use_internal_data = FALSE)
  kkdf <- as.data.frame(summary(kk))
  kkdf[1:5,1:6]
  
  
  ego <- enrichGO(gene          = genesego,
                  universe = GTuniverseEgo,
                  OrgDb         = org.Hs.eg.db, ## human
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 1, 
                  pvalueCutoff  = 1, 
                  readable      = TRUE)
  
  egodfBP=as.data.frame(summary(ego))
  egodfBP[1:10,1:6]
  
}


#### HELPER - Overlap CNV segments with genes

.assocGeneWithCNVSeg <- function(segs.pvalue.gr, gff.file, exclude.loc=TRUE){
  genes <- ape::read.gff(gff.file)
  genes <- genes[which(genes$type=="gene"),]
  genes$gene.name <- gsub(".*(gene=.*?);.*", "\\1", genes$attributes)
  genes$gene.name <- gsub("gene=", "", genes$gene.name)
  #setwd("/home/NIOO/viniciusd/great_tit/cnv/genomic_wave/refgen")
  chnam <- read.table("/home/NIOO/viniciusd/great_tit/geneinf/chrs_refseqs_GTBuild1.1.txt", sep="\t", header=T)
  genes <- merge(genes, chnam, by.x="seqid", by.y="RefSeq")
  
  genes.gr<- makeGRangesFromDataFrame(genes, seqnames.field = "Name", start.field = "start",
                                      end.field = "end", keep.extra.columns = TRUE)
  
  if(exclude.loc){
  genes.gr <- genes.gr[!grepl("LOC", genes.gr$gene.name)]
  }
  
  for(seg in 1:length(segs.pvalue.gr)){
  values(segs.pvalue.gr)$Overlapped.genes[seg] <- paste0(values(subsetByOverlaps(genes.gr, 
                                            segs.pvalue.gr[seg]))$gene.name, collapse=",")
  }
  return(segs.pvalue.gr)
}
