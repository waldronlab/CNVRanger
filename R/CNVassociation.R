##################
# Author: Vinicius Henrique da Silva
# Script description: Functions related to CNV-GWAS
# Date: March 21, 2018
# Code: CNVASSOPACK002
###################

#' HELPER - Create the folder tree to keep necessary files during and after the analysis 
#' @param name String with a project code or name
#' @param folder Choose manually the project folder (i.e. path as the root folder)
#' @return List with paths placed to store the files produced by subsequent analysis 
#' @examples 
#' all.paths <- createFolderTree("Project_name")

.createFolderTree <- function(name, folder=NULL){
  
  name <- paste0("/", name, "/")
  
  if(is.null(folder)){
    dir.create(paste0(user_data_dir("CNVAnalyzer"),name, "/Inputs"), recursive = TRUE)
    Inputs <- paste0(user_data_dir("CNVAnalyzer"), name, "/Inputs")
    dir.create(paste0(user_data_dir("CNVAnalyzer"), name, "/PLINK"), recursive = TRUE)
    PLINK <- paste0(user_data_dir("CNVAnalyzer"), name, "/PLINK")
    dir.create(paste0(user_data_dir("CNVAnalyzer"), name, "/Results"), recursive = TRUE)
    Results <- paste0(user_data_dir("CNVAnalyzer"), name, "/Results")}
  
  if(!is.null(folder)){
    dir.create(paste0(folder, "/CNVAnalyzer/", name, "/Inputs"), recursive = TRUE)
    Inputs <-  paste0(folder, "/CNVAnalyzer/", name, "/Inputs")
    dir.create(paste0(folder, "/CNVAnalyzer/", name, "/PLINK"), recursive = TRUE)
    PLINK <-  paste0(folder, "/CNVAnalyzer/", name, "/PLINK")
    dir.create(paste0(folder, "/CNVAnalyzer/", name, "/Results"), recursive = TRUE)
    Results <-  paste0(folder, "/CNVAnalyzer/", name, "/Results")
  }
  
  return(unlist(list(Inputs, PLINK, Results)))
}

#' HELPER - Download and test PLINK 1.07
#' @import utils
#' @import rappdirs
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree
#' @return A message of success (TRUE) of fail (FALSE) in running PLINK

.getPLINK <- function(all.paths){
  
  setwd(all.paths[2]) # Change to PLINK directory
  
  if(as.character(rappdirs:::get_os()) == "unix"){
    download.file("http://zzz.bwh.harvard.edu/plink/dist/plink-1.07-x86_64.zip", "PLINK.zip")
  }
  
  if(as.character(rappdirs:::get_os()) == "mac"){
    download.file("http://zzz.bwh.harvard.edu/plink/dist/plink-1.07-mac-intel.zip", "PLINK.zip")
  }
  
  if(as.character(rappdirs:::get_os()) == "win"){
    download.file("http://zzz.bwh.harvard.edu/plink/dist/plink-1.07-dos.zip", "PLINK.zip")
  }
  
  unzip("PLINK.zip")
  setwd(paste0(getwd(),"/", dir(path = ".", pattern = "plink*")))
  
  plinkPath <- paste0("\'", getwd(), "/plink\'")
  Sys.chmod(paste0(getwd(), "/plink"), mode = "0777", use_umask = TRUE)
  
  testP <- suppressWarnings(suppressMessages(as.character(grep("v1.07", system(paste(plinkPath, "--noweb"), wait=TRUE, 
                                                                               intern=TRUE, ignore.stderr = TRUE), value=TRUE))))
  return(0<grep("v1.07", testP))
}

#' HELPER - Load CNVs and phenotypes for each of the samples to be analyzed
#'
#' @param file.nam Name of the file to be imported
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree
#' @return List phen.info with samplesPhen, phenotypes, phenotypesdf, phenotypesSam, FamID and SexIds
#' @examples
#' phen.info <- loadCNVPhen("/home/.../CNVPhen.txt")
#' sapply(phen.info, class)

.loadCNVPhen <- function(file.nam, all.paths){
  setwd(all.paths[1]) # Change to PLINK directory
  all <- data.table::fread(file.nam, header=TRUE, sep="\t") ### Import all info
  
  if(!all(colnames(all)[1:3]==c("sample_id", "fam", "sex"))){
    base::stop("Unexpected first three column names")
  }
  
  samplesPhen <- base::unique(all$sample_id) ## Greb sample names
  phenotypes <- base::names(all[, 4:ncol(all)]) ## Greb phenotype names
  
  phenotypesdf <- base::as.data.frame(all[, 4:base::ncol(all), drop = FALSE]) ### Greb phenotypes values
  phenotypesdf <- base::as.data.frame(phenotypesdf[!base::duplicated(phenotypesdf),,drop = FALSE])
  phenotypesSam <- base::as.data.frame(base::cbind(samplesPhen, phenotypesdf))
  phenotypesSam$samplesPhen <- base::as.character(phenotypesSam$samplesPhen)
  
  FamID <- base::as.data.frame(base::cbind(samplesPhen, all$fam[base::as.numeric(base::rownames(phenotypesSam))]))
  SexIds <- base::as.data.frame(base::cbind(samplesPhen,all$sex[base::as.numeric(base::rownames(phenotypesSam))]))
  
  phen.info <- base::list("samplesPhen" = samplesPhen, "phenotypes" = phenotypes, "phenotypesdf" = phenotypesdf,
                          "phenotypesSam" = phenotypesSam, "FamID" = FamID, "SexIds" = SexIds)
  return(phen.info)
}

#' Setup the folders and files to apply the CNV-GWAS analysis
#'
#' @param name String with a project code or name
#' @param file.nam Name of the file to be imported
#' @param phen.loc Path to the tab separated text file cointaining phenotype and sample info
#' @param map.loc Path to the probe map used in PennCNV analysis
#' @param penn.out.loc Path to the PennCNV analysis output (i.e. CNVs in each sample)
#' @param folder Choose manually the project folder (i.e. path as the root folder)
#' @return List phen.info with samplesPhen, phenotypes, phenotypesdf, phenotypesSam, FamID, SexIds and all.paths
#' @examples
#' name <- "Project1"
#' phen.loc <- ".../Pheno.txt"
#' map.loc <- ".../Map.txt"
#' penn.out.loc <- ".../PennCNVOutput.txt"
#' phen.info <- setupCnvAna(name, phen.loc, map.loc, penn.out.loc)
#' @export

setupCnvAna <- function(name, phen.loc, map.loc, penn.out.loc, folder=NULL){
  
  ## Create the folder structure for all subsequent analysis
  all.paths <- .createFolderTree(name, folder)
  
  ## Import the phenotype and sample info from external folder
  file.copy(phen.loc, paste0(all.paths[1], "/Pheno.txt"))
  ## Import the probe map from external folder
  file.copy(map.loc,paste0(all.paths[1], "/MapPenn.txt"))
  ## Import the PennCNV output from external folder
  file.copy(penn.out.loc,paste0(all.paths[1], "/PennOut.txt"))
  
  file.nam <- "Pheno.txt"
  phen.info <- .loadCNVPhen(file.nam, all.paths)
  
  PlinkVer <- .getPLINK(all.paths)
  
  if(PlinkVer != "TRUE"){
    stop("PLINK setup failed")}
  
  setwd(all.paths[2])
  NewPLINK <- paste0(getwd(),"/", dir(path = ".", pattern = "plink*"))
  all.paths[2] <- NewPLINK
  
  phen.info$all.paths <- all.paths
  return(phen.info)
}

#' HELPER - Write CNV-genotype in a probe-wise fashion
#' 

.writeProbesCNV <- function(lo, all.samples, genofile, CNVsGr, probesCNVGr, n){
  sampleX <- all.samples[[lo]]
  g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id=sampleX))
  CNVSamByState <- split(CNVsGr[S4Vectors::values(CNVsGr)$V5 == all.samples[lo]], 
                         as.character(S4Vectors::elementMetadata(CNVsGr[S4Vectors::values(CNVsGr)$V5 == all.samples[lo]])$V4))
  for(slo in 1:length(CNVSamByState)){
    state <- unique(substr(x=as.character(S4Vectors::values(CNVSamByState[[slo]])$V4),
                           start=nchar(as.character(S4Vectors::values(CNVSamByState[[slo]])$V4)),
                           stop=nchar(as.character(S4Vectors::values(CNVSamByState[[slo]])$V4))))
    SamCNVSNP <- S4Vectors::values(subsetByOverlaps(probesCNVGr, CNVSamByState[[slo]]))$snp.id
    g[SamCNVSNP] <- state
  }
  g <- as.integer(g)
  gdsfmt::write.gdsn(n, g, start=c(1,lo), count=c(length(g),1))
}

#' HELPER Wait x seconds
#' from https://stat.ethz.ch/R-manual/R-devel/library/base/html/Sys.sleep.html

testit <- function(x)
{
  p1 <- proc.time()
  Sys.sleep(x)
  proc.time() - p1 
}

#' Produce CNV-GDS for the phenotyped samples
#'
#' @param phen.info Returned by \code{setupCnvAna}
#' @param lo.phe The phenotype to be analyzed in the PhenInfo$phenotypesSam data-frame
#' @param freq.cn Minimum frequency - 1 = 100% # Default 0.01 = 1%
#' @param snp.matrix Only FALSE implemented - If TRUE the SNP genotypes would be used to reconstruct CNV-SNP genotypes
#' @param n.cor Number of cores to be used
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree
#' @return
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples
#'  penn.out <- "PennCNVOutput.txt"
#'  map.prob <- "ProbeMap.txt"
#'  freq.cn <- 0.01
#'  snp.matrix <- "FALSE"
#'  n.cor <- 2
#'
#'  prodGdsCnv(phen.info, lo.phe, freq.cn, snp.matrix, n.cor, all.paths)
#'  @export

prodGdsCnv <- function(phen.info, freq.cn, snp.matrix, n.cor=1, lo.phe=1){  
  packageStartupMessage("initializing ...", appendLF = FALSE)
  
  phenotypesSam <- phen.info$phenotypesSam
  samplesPhen <- phen.info$samplesPhen
  FamID <- phen.info$FamID
  SexIds <- phen.info$SexIds
  all.paths <- phen.info$all.paths
  
  phenotypesSamX <- phenotypesSam[,c(1,(lo.phe+1))]
  phenotypesSamX <- na.omit(phenotypesSamX) ### Exclude samples without phenotypes. i.e NA
  
  ######################  Import CNVs to data-frame
  setwd(all.paths[1])
  cnvs <- read.table("PennOut.txt", sep="", header=F) ### PennCNV output
  label <- cnvs[,1]
  write.table(label, "del.txt", sep=";", col.names=F, row.names=F, quote=F)
  label <- read.table("del.txt", header=F, sep=":")
  write.table(label, "del.txt", sep="-", col.names=T, row.names=F, quote=F)
  label <- read.table("del.txt", header=F, sep="-", skip=1)
  cnvs <- cbind(label, cnvs)
  cnvs$V1 <- gsub("chr", "", cnvs$V1)
  colnames(cnvs)[1:3] <- c("chr", "start", "end")
  CNVs <- cnvs
  CNVs$start <- as.integer(as.character(CNVs$start))
  CNVs$end <- as.integer(as.character(CNVs$end))
  CNVsGr <- GenomicRanges::makeGRangesFromDataFrame(CNVs, keep.extra.columns = T)
  CNVsGr <- CNVsGr[values(CNVsGr)$V5 %in% samplesPhen] ### Subset CNVs in phenotyped samples
  
  ######################  Import SNP map to data-frame
  probes <- data.table::fread("MapPenn.txt", header=T, sep="\t")
  probes <- as.data.frame(probes)
  probes$Position <- as.numeric(as.character(probes$Position))
  probes <- probes[complete.cases(probes), ]
  probesGr <- GenomicRanges::makeGRangesFromDataFrame(probes, seqnames.field = "Chr", start.field = "Position", end.field = "Position",
                                       keep.extra.columns = T)
  
  all.samples <- unique(as.character(S4Vectors::values(CNVsGr)$V5))
  
  ###################### Select probes  within CNVs
  probesCNV <- suppressWarnings(as.character(S4Vectors::values(GenomicRanges::subsetByOverlaps(probesGr, 
                                                                     CNVsGr))$Name))
  probesCNV <- unique(unlist(probesCNV))
  probesCNVGr <- probesGr[S4Vectors::values(probesGr)$Name %in% probesCNV]
  
  counts <- as.integer(GenomicRanges::countOverlaps(probesCNVGr, CNVsGr))
  S4Vectors::values(probesCNVGr)$freq <- counts
  
  ##### Subset by frequency 
  NumSam <- freq.cn*length(all.samples)
  probesCNVGr <- probesCNVGr[S4Vectors::values(probesCNVGr)$freq>=NumSam]
  
  ##### Order the probes
  probesCNVGr <- GenomeInfoDb::sortSeqlevels(probesCNVGr)
  probesCNVGr <- sort(probesCNVGr)
  
  S4Vectors::values(probesCNVGr)$snp.id <- seq(1:length(probesCNVGr))
  
  ###################### SNP genotype matrix not available 
  if(snp.matrix == "FALSE"){
    CNVBiMa <- matrix(2, nrow = length(probesCNVGr), 
                      ncol = length(all.samples))}
  
  ###################### SNP genotype matrix available - TODO
  if(snp.matrix == "TRUE"){
    ...
  }
  
  ###################### Create a GDS with chr and SNP names to numeric
  SNPRelate::snpgdsCreateGeno("CNV.gds", 
                   genmat = CNVBiMa,
                   sample.id = all.samples, 
                   snp.id = as.character(values(probesCNVGr)$snp.id),
                   snp.rs.id = as.character(values(probesCNVGr)$Name),
                   snp.chromosome = as.character(seqnames(probesCNVGr)),
                   snp.position = start(probesCNVGr)
  )
  
  testit(15)
  
  ###################### Replace genotype matrix with CNV genotypes
  (genofile <- SNPRelate::snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE))
  
  CNVBiMaCN <- matrix(2, nrow = length(probesCNVGr), 
                      ncol = length(all.samples))
  n <- gdsfmt::add.gdsn(genofile, "CNVgenotype", CNVBiMaCN, replace=TRUE)
  read.gdsn(n)
  
  gback <- as.numeric(CNVBiMaCN[,1])
  
  if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
    multicoreParam <-  BiocParallel::MulticoreParam(workers = n.cor)
    suppressMessages(BiocParallel::bplapply(1:length(all.samples), .writeProbesCNV, BPPARAM = multicoreParam, all.samples=all.samples,
                              genofile=genofile, CNVsGr=CNVsGr, probesCNVGr=probesCNVGr, n=n))
  }
  
  if(rappdirs:::get_os() == "win"){
    param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
    suppressMessages(BiocParallel::bplapply(1:length(all.samples), .writeProbesCNV, BPPARAM = param, all.samples=all.samples, 
                              genofile=genofile, CNVsGr=CNVsGr, probesCNVGr=probesCNVGr, n=n))
  }
  
  testit(15) ## Wait to make sure that gds is writen
  
  packageStartupMessage("GDS writen...", appendLF = FALSE)
  
  ###################### Include the phenotype 'lo' in the GDS
  phenotypesSamX <- phenotypesSamX[match(all.samples, phenotypesSamX$samplesPhen),] ### Check order correspondence
  gdsfmt::add.gdsn(genofile, name="phenotype", val=phenotypesSamX, replace=TRUE)
  FamID <- FamID[match(all.samples, FamID$samplesPhen),] ### Check order correspondence
  gdsfmt::add.gdsn(genofile, name="FamID", val=FamID[,2], replace=TRUE)
  FamID <- SexIds[match(all.samples, SexIds$samplesPhen),] ### Check order correspondence
  gdsfmt::add.gdsn(genofile, name="Sex", val=SexIds[,2], replace=TRUE)
  SNPRelate::snpgdsClose(genofile)
  ########################################### END #######################################
  return(probesCNVGr)
  packageStartupMessage("done", appendLF = FALSE)
}

#' HELPER - Produce the .gvar file requested for the CNV-GWAS in PLINK
#' @examples
#' Gens <- .prodGvar(lo, genofile, sam.gen, snps, fam.id)

.prodGvar <- function(lo, genofile, sam.gen, snps, fam.id){
  g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id=sam.gen[[lo]]))
  A <- rep("A", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
  B <- rep("B", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
  CNVg <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))
  CNVg <- CNVg-1
  Dose1 <- rep(1, length(CNVg))
  Dose1[CNVg==-1] <- 0
  CNVg[CNVg==-1] <- 0
  B[CNVg==1] <- "A"
  
  if(SNPmatrix=="FALSE"){
    df <- as.data.frame(cbind(as.character(fam.id[[lo]]), sam.gen[[lo]], snps, 
                              A, CNVg, B, Dose1))
    colnames(df) <- c("FID", "IID", "NAME", "ALLELE1", "DOSAGE1", "ALLELE2", "DOSAGE2")}
  
  if(SNPmatrix=="TRUE"){ #### TODO - To implement
    ...
  }
  return(df)
}

#' Run the CNV-GWAS 
#' 
#' (i) Produces the requested inputs for a PLINK analysis and (ii) run a 
#' CNV-GWAS analysis using the GDS wrote by the last \code{\link{prodGdsCnv}}
#' run.
#'
#' @param phen.info Returned by \code{\link{setupCnvAna}}
#' @param n.cor Number of cores to be used
#' @param method.m.test Correction for multiple tests to be used. FDR is default, see \code{\link(p.adjust)} for
#' other methods.
#' @param min.sim Minimum CNV genotype distribuition similarity among subsequent probes. Default is 0.95 = 95%
#' @return
#' @examples
#' @export

cnvGWAS <- function(phen.info, n.cor, probesCNVGr, min.sim = 0.95, method.m.test = "fdr", lo.phe=1){
  phenotypesSam <- phen.info$phenotypesSam
  phenotypesSamX <- phenotypesSam[,c(1,(lo.phe+1))]
  phenotypesSamX <- na.omit(phenotypesSamX) ### Exclude samples without phenotypes. i.e NA
  all.paths <- phen.info$all.paths
  ##################################### (iii) Produce PLINK map ##################
  setwd(all.paths[1])
  (genofile <- SNPRelate::snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE))
  Pmap <- as.data.frame(cbind(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                              gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
                              0, gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position"))))
  
  colnames(Pmap) <- c("Chr", "NAME", "GD", "Position")
  setwd(all.paths[2])
  write.table(Pmap, "mydata.map", sep="\t", col.names=T, row.names=F, quote=F)
  SNPRelate::snpgdsClose(genofile)
  ########################################### END #######################################
  
  ##################################### (iv) Produce gvar to use as PLINK input #########
  setwd(all.paths[1])
  (genofile <- SNPRelate::snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE)) ## read GDS
  sam.gen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")) ## Extract samples
  snps <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")) ## Extract SNPs
  fam.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "FamID")) ## Extract family
  
  ## Produce gvar file in parallel processing
  # Identify the SO
  if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
    multicoreParam <- BiocParallel::MulticoreParam(workers = n.cor)
    Gens <- BiocParallel::bplapply(1:length(sam.gen), .prodGvar, BPPARAM = multicoreParam, genofile=genofile, sam.gen=sam.gen,
                     snps=snps, fam.id=fam.id)
  }
  
  if(rappdirs:::get_os() == "win"){
    param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
    Gens <- BiocParallel::bplapply(1:length(sam.gen), .prodGvar, BPPARAM = param, genofile=genofile, sam.gen=sam.gen,
                     snps=snps, fam.id=fam.id)
  }
  
  gentype.all <- data.table::rbindlist(Gens)
  nrow(gentype.all)
  gentype.all1 <- as.data.frame(gentype.all)
  
  setwd(all.paths[2])
  write.table(gentype.all1, "mydata.gvar", sep="\t", col.names=T, row.names=F, quote=F)
  SNPRelate::snpgdsClose(genofile)
  ########################################### END #######################################
  
  ##################################### (v) Produce fam (phenotype) to use as PLINK input #############
  setwd(all.paths[1])
  (genofile <- SNPRelate::snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE)) ## read GDS
  SamGen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")) ## Extract samples
  SNPs <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")) ## Extract SNPs
  FAMID <- as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "FamID"))) ## Extract family
  Sex <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Sex")) ## Extract sex
  Phen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")) ## Extract phenotype
  
  setwd(all.paths[2])
  famdf <- as.data.frame(cbind(FAMID, SamGen, SamGen, "NC", Sex, Phen[,2]))
  colnames(famdf) <- c("FAID", "IID", "PID", "MID", "Sex", "Phenotype")
  
  write.table(famdf, "mydata.fam", sep = "\t", col.names = T, row.names=F, quote=F)
  SNPRelate::snpgdsClose(genofile)
  ########################################### END #######################################
  
  ##################################### (vi) Run PLINK #############
  setwd(all.paths[2])
  plinkPath <- paste0("\'", getwd(), "/plink\'")
  system(paste(plinkPath, "--gfile mydata --noweb"), wait=TRUE)
  
  setwd(all.paths[1])
  (genofile <- SNPRelate::snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE)) ## read GDS
  sumphen <- paste0("plink.gvar.summary.",names(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")))[[2]])
  setwd(all.paths[2])
  file.rename("plink.gvar.summary", sumphen)
  SNPRelate::snpgdsClose(genofile)
  ########################################### END #######################################
  
  ##################################### (vii) Produce CNV segments #############
  setwd(all.paths[1])
  (genofile <- SNPRelate::snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE)) ## read GDS
  
  chrx <- as.data.frame(cbind(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                              seq(from=1, to=length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome"))), by=1) ))
  g1 <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"))
  g2 <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"))
  
  chrx$V1 <- factor(chrx$V1, levels=unique(chrx$V1))
  chrx.split <- split(chrx, chrx$V1)
  
  ### Find similarity between subsequent probes 
  SimiLar.all <- NULL
  for(lo in 1:length(chrx.split)){
    snpindex <- as.numeric(as.character(chrx.split[[lo]]$V2))
    g1x <- g1[snpindex,]
    g1xM <- g1x[-(nrow(g1x)),]
    g1CN <- apply(g1xM, 1, function (x) sum(table(x)[names(table(x))!=2]))
    g1x <- as.data.frame(t(sapply(1:nrow(g1x), function(i) replace(g1x[i,], g1x[i,] == 2, paste0(i, 'R')))))
    g2x <- g2[snpindex,]
    g2x <- as.data.frame(t(sapply(1:nrow(g2x), function(i) replace(g2x[i,], g2x[i,] == 2, paste0(i, 'R')))))
    g1x <- g1x[-(nrow(g1x)),]
    g2x <- g2x[-1,]
    ftma <- g1x==g2x
    SimiLar <- apply(ftma, 1, function(x) as.numeric(table(x)["TRUE"]))
    SimiLar[is.na(SimiLar)] <- 0
    SimiLar.all[[lo]] <- as.numeric(SimiLar)/g1CN ### CNV genotype similarity between subsequent SNPs
  }
  
  all.segs <- NULL
  for(ch in 1:length(SimiLar.all)){
    SimiLarX <- SimiLar.all[[ch]]
    nextP <- NULL
    all.segsch <- NULL
    for(se in 1:length(SimiLarX)){
        if(SimiLarX[[se]]>=min.sim){
        nextP[[se]] <- "TRUE"
      }
      if(SimiLarX[[se]]<min.sim){
        nextP[[se]] <- "FALSE"  
      }
      
      if(se ==1){
        all.segsch[[se]] <- "start"
      }
      
      if(se >1){
        if(nextP[[se-1]]=="TRUE"){
          all.segsch[[se]] <- "cont"
        }
        if(nextP[[se-1]]=="FALSE"){
          all.segsch[[se]] <- "start"
        }
      }
      
      if(se ==length(SimiLarX)){
        if(SimiLarX[[se]]>=min.sim){
          all.segsch[[se+1]] <- "start"}
        if(SimiLarX[[se]]<min.sim){
          all.segsch[[se+1]] <- "cont"}
      }
      
    }
    all.segs[[ch]] <- all.segsch
  }
  
  S4Vectors::values(probesCNVGr)$starts.seg <- unlist(all.segs)
  
  pstarts <- probesCNVGr[S4Vectors::values(probesCNVGr)$starts.seg == "start"]
  pstarts.split = split(pstarts, GenomeInfoDb::seqnames(pstarts))
  
  all.seqsGr <- GenomicRanges::GRangesList()
  for(chx in 1:length(pstarts.split)){
    pstartsx <- pstarts.split[[chx]]
    if(length(pstartsx)==0){
      all.seqsGr[[chx]] <- pstartsx
      next
    }
    all.seqsGrX <- GenomicRanges::GRangesList()
    
    for(st in 1:length(pstartsx)){
      prx <- pstartsx[st]
      if(st < length(pstartsx)){
        LimPrx <- BiocGenerics::start(pstartsx[st+1])-1}
      if(st == length(pstartsx)){
        LimPrx <-  max(BiocGenerics::start(probesCNVGr[GenomeInfoDb::seqnames(probesCNVGr) == GenomeInfoDb::seqnames(prx)]))
      }
      seqLar <- GenomicRanges::GRanges(as.character(GenomeInfoDb::seqnames(prx)), IRanges::IRanges(BiocGenerics::start(prx), LimPrx))
      seqPrbs <- GenomicRanges::subsetByOverlaps(probesCNVGr,seqLar)
      seqx <- GenomicRanges::GRanges(as.character(GenomeInfoDb::seqnames(prx)), IRanges::IRanges(BiocGenerics::start(prx), BiocGenerics::start(seqPrbs[length(seqPrbs)])))
      all.seqsGrX[[st]] <- seqx
    }
    all.seqsGr[[chx]] <- unlist(all.seqsGrX)
  }
  
  snpgdsClose(genofile)
  
  ########################################### END #######################################
  
  ##################################### (viii) Associate SNPs with CNV segments #############
  SegsGr <- unlist(all.seqsGr)
  S4Vectors::values(SegsGr)$SegName <- seq(from=1, to=length(SegsGr), by=1)
  
  setwd(all.paths[1])
  (genofile <- SNPRelate::snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE))
  
  setwd(all.paths[2])
  results <- read.table(paste0("plink.gvar.summary.",names(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")))[[2]]), header=T) 
  mydata.map <- read.table("mydata.map", header=T)
  resultsp <- merge(results, mydata.map, by.x="NAME", by.y="NAME", all.x=T, all.y=F)
  
  #### Choose the method - TODO - Change if SNPMatrix implemented - Allow to get other p-values
  resultsp <- resultsp[ grep("P\\(CNP\\)", resultsp$FIELD),]
  resultsp$VALUE <- as.numeric(as.character(resultsp$VALUE))
  resultsp <- resultsp[order(resultsp$VALUE, decreasing=F),]
  resultspGr <- GenomicRanges::makeGRangesFromDataFrame(resultsp, seqnames.field = "Chr",
                                         start.field = "Position", end.field = "Position",
                                         keep.extra.columns = T)
  
  ######## Assign the lowest p-value to the CNV segment
  
  values.all <- vector("list", length(SegsGr))
  names.all <- vector("list", length(SegsGr))
  for(se in 1:length(SegsGr)){
    values.all[[se]] <- min(S4Vectors::values(GenomicRanges::subsetByOverlaps(resultspGr, SegsGr[se]))$VALUE)
    Prbx <- GenomicRanges::subsetByOverlaps(resultspGr, SegsGr[se])
    names.all[[se]] <- as.character(Prbx[values(Prbx)$VALUE %in% min(values(GenomicRanges::subsetByOverlaps(resultspGr, SegsGr[se]))$VALUE)]$NAME[[1]])
  }
  
  S4Vectors::values(SegsGr)$MinPvalue <- unlist(values.all)
  S4Vectors::values(SegsGr)$NameProbe <- unlist(names.all)
  
  S4Vectors::values(SegsGr)$MinPvalueAdjusted <- stats::p.adjust(values(SegsGr)$MinPvalue, method =method.m.test, n = length(SegsGr))
  SegsGr$Phenotype <- names(phenotypesSamX)[[2]]
  SNPRelate::snpgdsClose(genofile)
  ######################################################## END ###########################################
  
  return(SegsGr)
}

