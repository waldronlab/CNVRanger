#' Import LRR and BAF from text files used in the CNV analysis
#'
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree 
#' @param path.files Folder containing the input CNV files used for the CNV calling (i.e. one text file with 5 collumns 
#' for each sample). Columns should contain (i) probe name, (ii) Chromosome, (iii) Position, (iv) LRR and (v) BAF
#' @param list.of.files Data-frame with two columns where the (i) is the file name with signals and (ii) is the 
#' correspondent name of the sample in the gds file
#' @param n.cor Number of cores to be used
importLRR_BAF <- function(all.paths, path.files, list.of.files, n.cor){
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  
  ### Create file path for all text files
  file.paths <- file.path(path.files, list.of.files[,1])
  gds.names <- list.of.files[,2]
  list.filesLo <- cbind("add"=file.paths, "gds"=gds.names)
  list.filesLo <- as.data.frame(list.filesLo)
  
  snps.included <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")))
  snps.included <- as.character(snps.included)
  
  all.samples <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
  all.samples <- as.character(all.samples)
  
  ## Check if all files 
  list.filesLo.back <- list.filesLo
  list.filesLo <- list.filesLo[which(list.filesLo$gds %in% all.samples),]
  if(nrow(list.filesLo) != nrow(list.filesLo.back)){
    message("Warning: list.of.files has different length of the list of samples from gds")
  }
  
  LRR.matrix <- matrix(0, nrow = length(snps.included), 
                      ncol = length(all.samples))
  
  nLRR <- gdsfmt::add.gdsn(genofile, "LRR", LRR.matrix, replace=TRUE)
  gdsfmt::read.gdsn(nLRR)
  
  
  BAF.matrix <- matrix(0.5, nrow = length(snps.included), 
                       ncol = length(all.samples))
  
  nBAF <- gdsfmt::add.gdsn(genofile, "BAF", BAF.matrix, replace=TRUE)
  gdsfmt::read.gdsn(nBAF)
  
  if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
    multicoreParam <-  BiocParallel::MulticoreParam(workers = n.cor)
    message("Start the parallel import of LRR/BAF values")
    BiocParallel::bplapply(1:nrow(list.filesLo), .freadImport, BPPARAM = multicoreParam, list.filesLo=list.filesLo,
                           genofile=genofile, all.samples=all.samples, nLRR=nLRR, nBAF=nBAF, snps.included=snps.included)
                           
  }
  
  if(rappdirs:::get_os() == "win"){
    message("Start the parallel import of LRR/BAF values")
    param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
    BiocParallel::bplapply(1:nrow(list.filesLo), .freadImport, BPPARAM = param, list.filesLo=list.filesLo,
                           genofile=genofile, all.samples=all.samples, nLRR=nLRR, nBAF=nBAF, snps.included=snps.included)
  }
  
  SNPRelate::snpgdsClose(genofile)
  
} ## BUG - Parallel not working


# HELPER - Function to apply during LRR/BAF import
# @param lo loop number, based on the number of SNPs per turn
# @param list.filesLo Data-frame with two columns where the (i) is the path + file name with signals and (ii) is the 
# correspondent name of the sample in the gds file
# @param genofile loaded gds file
# @param all.samples All samples in the gds file
# @param nLRR Connection to write LRR values
# @param nBAF Connection to write BAF values
# @param snps.included All SNP probe names to be included
.freadImport <- function(lo, list.filesLo, genofile, all.samples, nLRR=NULL, nBAF=NULL, snps.included){
  
  message(paste0("sample ", lo, " of ", length(all.samples)))
  
  file.x <- c(as.character(list.filesLo[lo,1]), as.character(list.filesLo[lo,2]))
  
  ### Import the sample file with fread
  sig.x <- data.table::fread(file.x[1], skip=1L, header=FALSE)
  sig.x <- as.data.frame(sig.x)
  colnames(sig.x) <- c("name", "chr", "position", "lrr", "baf")
  
  ### Order the signal file as in the gds
  sig.x <- subset(sig.x, name %in% snps.included)
  
  ### Include missing SNPs
  missing.snps <- snps.included[-which(snps.included %in% sig.x$name)]
  if(length(missing.snps)>0){
  missing.snps <- as.data.frame(missing.snps)
  missing.snps <- cbind(missing.snps, "chr"= "NoN", "position"= "NoN", "lrr"=NA, "baf"=NA)
  colnames(missing.snps)[1] <- "name"
  sig.x <- rbind(sig.x, missing.snps)
  }
  
  sig.x <- sig.x[order(match(sig.x[,1],snps.included)),]
  
  ### Extract the LRR
  LRR.x  <- sig.x$lrr
  
  ### Extract the BAF
  BAF.x <- sig.x$baf
  
  ### Find the number of the sample 
  sam.num <- which(all.samples == file.x[2])
  
  if(!is.null(nLRR)){
  ### Write LRR value for sample x
  gdsfmt::write.gdsn(nLRR, LRR.x, start=c(1,sam.num), count=c(length(snps.included),1))
  }
  else if(!is.null(nBAF)){
  ### Write BAF value for sample x
  gdsfmt::write.gdsn(nBAF, BAF.x, start=c(1,sam.num), count=c(length(snps.included),1))}
  
} ## BUG - Parallel not working


# HELPER - Formula to infer Vst
.vstForm <- function(value, pops.names){

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
calcVst <- function(all.paths, pops, chr.code.name=NULL){
  
  #(genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNVLRRBAFAllSNPs.gds"), allow.fork=TRUE, readonly=FALSE))
  
  pops.names <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Pops.names")))  
  
  if(all(pops %in% pops.names)){
  
  value <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "LRR")))
  value <- value[,which(pops.names %in% pops)]
  
  pops.names <- pops.names[pops.names %in% pops]
  
  vst.result <- apply(value, 1, FUN=.vstForm, pops.names=pops.names)
  
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
                    EV1 = pca$eigenvect[,1],
                    EV2 = pca$eigenvect[,2],
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

##################
# Author: Vinicius Henrique da Silva
# Script description: Functions related to population comparison and evolutionary analysis
# Date: April 19, 2018
# Code: CNVEVOPACK001
###################

#' HELPER - Function to apply during LRR/BAF import
#' @param lo loop number, based on the number of SNPs per turn
#' @param list.filesLo Data-frame with two columns where the (i) is the path + file name with signals and (ii) is the 
#' correspondent name of the sample in the gds file
#' @param genofile loaded gds file
#' @param all.samples All samples in the gds file
#' @param nLRR Connection to write LRR values
#' @param nBAF Connection to write BAF values
#' @param snps.included All SNP probe names to be included

.freadImport <- function(lo, list.filesLo, genofile, all.samples, nLRR=NULL, nBAF=NULL, snps.included){
  
  message(paste0("sample ", lo, " of ", length(all.samples)))
  
  file.x <- c(as.character(list.filesLo[lo,1]), as.character(list.filesLo[lo,2]))
  
  ### Import the sample file with fread
  sig.x <- data.table::fread(file.x[1], skip=1L, header=FALSE)
  sig.x <- as.data.frame(sig.x)
  colnames(sig.x) <- c("name", "chr", "position", "lrr", "baf")
  
  ### Order the signal file as in the gds
  sig.x <- subset(sig.x, name %in% snps.included)
  
  ### Include missing SNPs
  missing.snps <- snps.included[-which(snps.included %in% sig.x$name)]
  if(length(missing.snps)>0){
  missing.snps <- as.data.frame(missing.snps)
  missing.snps <- cbind(missing.snps, "chr"= "NoN", "position"= "NoN", "lrr"=NA, "baf"=NA)
  colnames(missing.snps)[1] <- "name"
  sig.x <- rbind(sig.x, missing.snps)
  }
  
  sig.x <- sig.x[order(match(sig.x[,1],snps.included)),]
  
  ### Extract the LRR
  LRR.x  <- sig.x$lrr
  
  ### Extract the BAF
  BAF.x <- sig.x$baf
  
  ### Find the number of the sample 
  sam.num <- which(all.samples == file.x[2])
  
  if(!is.null(nLRR)){
  ### Write LRR value for sample x
  gdsfmt::write.gdsn(nLRR, LRR.x, start=c(1,sam.num), count=c(length(snps.included),1))
  }
  else if(!is.null(nBAF)){
  ### Write BAF value for sample x
  gdsfmt::write.gdsn(nBAF, BAF.x, start=c(1,sam.num), count=c(length(snps.included),1))}
  
} ## BUG - Parallel not working

#' Import LRR and BAF from text files used in the CNV analysis
#'
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree 
#' @param path.files Folder containing the input CNV files used for the CNV calling (i.e. one text file with 5 collumns 
#' for each sample). Columns should contain (i) probe name, (ii) Chromosome, (iii) Position, (iv) LRR and (v) BAF
#' @param list.of.files Data-frame with two columns where the (i) is the file name with signals and (ii) is the 
#' correspondent name of the sample in the gds file
#' @param n.cor Number of cores to be used

importLRR_BAF <- function(all.paths, path.files, list.of.files, n.cor){
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  
  ### Create file path for all text files
  file.paths <- file.path(path.files, list.of.files[,1])
  gds.names <- list.of.files[,2]
  list.filesLo <- cbind("add"=file.paths, "gds"=gds.names)
  list.filesLo <- as.data.frame(list.filesLo)
  
  snps.included <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")))
  snps.included <- as.character(snps.included)
  
  all.samples <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
  all.samples <- as.character(all.samples)
  
  ## Check if all files 
  list.filesLo.back <- list.filesLo
  list.filesLo <- list.filesLo[which(list.filesLo$gds %in% all.samples),]
  if(nrow(list.filesLo) != nrow(list.filesLo.back)){
    message("Warning: list.of.files has different length of the list of samples from gds")
  }
  
  LRR.matrix <- matrix(0, nrow = length(snps.included), 
                      ncol = length(all.samples))
  
  nLRR <- gdsfmt::add.gdsn(genofile, "LRR", LRR.matrix, replace=TRUE)
  gdsfmt::read.gdsn(nLRR)
  
  
  BAF.matrix <- matrix(0.5, nrow = length(snps.included), 
                       ncol = length(all.samples))
  
  nBAF <- gdsfmt::add.gdsn(genofile, "BAF", BAF.matrix, replace=TRUE)
  gdsfmt::read.gdsn(nBAF)
  
  if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
    multicoreParam <-  BiocParallel::MulticoreParam(workers = n.cor)
    message("Start the parallel import of LRR/BAF values")
    BiocParallel::bplapply(1:nrow(list.filesLo), .freadImport, BPPARAM = multicoreParam, list.filesLo=list.filesLo,
                           genofile=genofile, all.samples=all.samples, nLRR=nLRR, nBAF=nBAF, snps.included=snps.included)
                           
  }
  
  if(rappdirs:::get_os() == "win"){
    message("Start the parallel import of LRR/BAF values")
    param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
    BiocParallel::bplapply(1:nrow(list.filesLo), .freadImport, BPPARAM = param, list.filesLo=list.filesLo,
                           genofile=genofile, all.samples=all.samples, nLRR=nLRR, nBAF=nBAF, snps.included=snps.included)
  }
  
  SNPRelate::snpgdsClose(genofile)
  
} ## BUG - Parallel not working


#' HELPER - Formula to infer Vst
#' 
#' 

.vstForm <- function(value, pops.names){

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

calcVst <- function(all.paths, pops, chr.code.name=NULL){
  
  #(genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNVLRRBAFAllSNPs.gds"), allow.fork=TRUE, readonly=FALSE))
  
  pops.names <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Pops.names")))  
  
  if(all(pops %in% pops.names)){
  
  value <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "LRR")))
  value <- value[,which(pops.names %in% pops)]
  
  pops.names <- pops.names[pops.names %in% pops]
  
  vst.result <- apply(value, 1, FUN=.vstForm, pops.names=pops.names)
  
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


