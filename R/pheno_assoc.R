##################
# Author: Vinicius Henrique da Silva
# Script description: Functions related to CNV-GWAS
# Date: April 12, 2018
# Code: CNVASSOPACK002
###################

#' HELPER - Create the folder tree to keep necessary files during and after the analysis 
#' @param name String with a project code or name
#' @param folder Choose manually the project folder (i.e. path as the root folder)
#' @return List with paths placed to store the files produced by subsequent analysis 
#' @examples 
#' all.paths <- createFolderTree("Project_name")

.createFolderTree <- function(name, folder=NULL){
  
  # data dir
  if(is.null(folder)) folder <- rappdirs::user_data_dir("cnvAnalyzeR")
  if(!file.exists(folder)) dir.create(folder, recursive=TRUE)
  
  # project directory
  proj.dir <- file.path(folder, name)
  if(!file.exists(proj.dir)) dir.create(proj.dir)
  
  # subdirs
  dirs <- c("Inputs", "PLINK", "Results")
  all.paths <- file.path(proj.dir, dirs)
  for(d in all.paths) if(!file.exists(d)) dir.create(d)
  
  return(all.paths)
}

#' HELPER - Download and test PLINK 1.07
#' @import utils
#' @import rappdirs
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree
#' @return A message of success (TRUE) of fail (FALSE) in running PLINK

.getPLINK <- function(plink.path, version="1.07"){
  
  plink.url <- "http://zzz.bwh.harvard.edu/plink/dist/plink-"
  plink.url <- paste0(plink.url, version, "-")
  
  plink.file <- file.path(plink.path, "PLINK.zip")
  
  # get os-specific binary 
  os <- rappdirs:::get_os()
  os <- switch(os,
               unix = "x86_64",
               mac = "mac-intel",
               win = "dos",
               os)
  
  plink.url <- paste0(plink.url, os, ".zip")
  
  # download & unzip
  download.file(plink.url, plink.file)
  unzip(plink.file, exdir=plink.path)
  
  # remove zip
  file.remove(plink.file)
  
  # path to plink binary
  plink.file <- dir(plink.path, pattern = "plink*")
  plink.path <- file.path(plink.path, plink.file, "plink")
  Sys.chmod(plink.path, mode = "755")
  
  suppressWarnings(
    res <- system2(plink.path, args="--noweb", stdout=TRUE, stderr=FALSE)
  )
  
  res <- any( grepl(version, res) )
  return(res)
}

#' HELPER - Load phenotypes for each of the samples to be analyzed
#'
#' @param file.nam Name of the file to be imported
#' @param pheno.path First item of the list returned by \code{CreateFolderTree} function
#' @return List phen.info with samplesPhen, phenotypes, phenotypesdf, phenotypesSam, FamID and SexIds
#' @examples
#' phen.info <- loadPhen("/home/.../CNVPhen.txt")
#' sapply(phen.info, class)

.loadPhen <- function(file.nam, pheno.path){
  
  all <- data.table::fread(file.path(pheno.path, file.nam), header=TRUE, sep="\t") ### Import phenotypes
  
  if(!all(colnames(all)[1:3]==c("sample.id", "fam", "sex"))){
    stop("Unexpected first three column names")}
  
  samplesPhen <- unique(all$sample.id) ## Greb sample names
  phenotypes <- names(all[, 4:ncol(all)]) ## Greb phenotype names
  
  phenotypesdf <- all[, 4:ncol(all), drop = FALSE]
  phenotypesdf <- as.data.frame(phenotypesdf) ### Greb phenotypes values
  phenotypesdf <- as.data.frame(phenotypesdf[!duplicated(phenotypesdf),,drop = FALSE])
  
  phenotypesSam <- data.frame(samplesPhen, phenotypesdf)
  phenotypesSam$samplesPhen <- as.character(phenotypesSam$samplesPhen)
  
  FamID <- data.frame(samplesPhen, "V2"=all$fam[as.numeric(rownames(phenotypesSam))])
  
  SexIds <- data.frame(samplesPhen,"V2"=all$sex[as.numeric(rownames(phenotypesSam))])
  
  phen.info <- list("samplesPhen" = samplesPhen, "phenotypes" = phenotypes, "phenotypesdf" = phenotypesdf,
                    "phenotypesSam" = phenotypesSam, "FamID" = FamID, "SexIds" = SexIds)
  return(phen.info)
}

#' HELPER - Check and convert CNV input
#' @param cnvs Data-frame with the CNVs to be analyzed from PennCNV, general format or sequencing
#' 

.CheckConvertCNVs <- function(cnvs){
  
  ## If PennCNV file
  if(unique(sub("^([[:alpha:]]*).*", "\\1", cnvs$V2))=="numsnp"){
    
    df1 <- reshape2::colsplit(cnvs$V1, ":", c("chr", "loc"))
    df2 <- reshape2::colsplit(df1$loc, "-", c("start", "end"))
    df1 <- cbind(df1[,-2, drop=FALSE], df2)
    
    cnvs <- cbind(df1, cnvs)
    cnvs$V1 <- gsub("chr", "", cnvs$V1)
    colnames(cnvs)[1:3] <- c("chr", "start", "end")
    CNVs <- cnvs
    
    ## If general format
  } else if(all(colnames(cnvs))==c("chr", "start", "end","sample.id", "state", "num.snps", "start.probe", "end.probe")){ 
    
    ## Convert to PennCNV format
    cnvs$V1 <- paste0(cnvs$chr, ":",cnvs$start, "-", cnvs$end)
    cnvs$length <- (cnvs$end - cnvs$start)+1
    cnvs$length <- paste0("length=", cnvs$length)
    
    cnvs$num.snps <- paste0("numsnp=", cnvs$length)
    
    cnvs$start.probe <- paste0("startSNP=", cnvs$start.probe)
    cnvs$end.probe <- paste0("endSNP=", cnvs$start.probe)
    
    cnvs <- subset(cnvs, select =c(chr, start, end, V1, num.snps, length, state, sample.id, start.probe, end.probe))
    colnames(cnvs) <- c("chr", "start", "end", "V1", "V2", "V3", "V4", "V5", "V6", "V7")
    CNVs <- cnvs
    
    ## If sequencing info
  } else if(all(colnames(cnvs))==c("chr", "start", "end", "sample.id", "state")){
    ## Convert to PennCNV 
    
    ## TODO 
    
    
  } else stop("Unexpected CNV input format - is it tab delimited?")
  
  return(CNVs)
}

#' Setup the folders and files to apply the CNV-GWAS analysis
#'
#' This function creates the (i) necessary folders in disk to perform downstream analysis on CNV genome-wide
#' association and (ii) import the necessary input files (i.e phenotypes, probe map and CNV list) from other 
#' locations in disk. 
#' 
#' @param name String with a project code or name
#' @param phen.loc Path to the tab separated text file cointaining phenotype and sample info
#' @param cnv.out.loc Path to the CNV analysis output (i.e. PennCNV output or tab delimited text file)
#' @param map.loc Path to the probe map (e.g. used in PennCNV analysis)
#' @param folder Choose manually the project folder (i.e. path as the root folder)
#' @return List phen.info with samplesPhen, phenotypes, phenotypesdf, phenotypesSam, FamID, SexIds and all.paths
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples
#' name <- "Project1"
#' phen.loc <- ".../Pheno.txt"
#' map.loc <- ".../Map.txt"
#' cnv.out.loc <- ".../CNVOutput.txt"
#' phen.info <- setupCnvGWAS(name, phen.loc, map.loc, cnv.out.loc)
#' @export

setupCnvGWAS <- function(name, phen.loc, cnv.out.loc, map.loc=NULL, folder=NULL){
  
  ## Create the folder structure for all subsequent analysis
  all.paths <- .createFolderTree(name, folder)
  
  ## Import the phenotype and sample info from external folder
  file.copy(phen.loc, file.path(all.paths[1], "/Pheno.txt"))
  ## Import the PennCNV output from external folder
  file.copy(cnv.out.loc, file.path(all.paths[1], "/CNVOut.txt"))
  
  if(is.null(map.loc)){
    ## Import the probe map from external folder
    cnvs <- read.table(file.path(all.paths[1], "CNVOut.txt"), sep="", header=F) ### CNV table 
    CNVs <- .CheckConvertCNVs(cnvs)
    CGr <- makeGRangesFromDataFrame(CNVs)
    startP <- cbind(CNVs$chr, start(CGr))
    endP <- cbind(CNVs$chr, end(CGr))
    middleP <- cbind(CNVs$chr, end(CGr)-(width(CGr)/2))
    probe.like.map <- rbind(startP, endP, middleP)
    probe.like.map <- as.data.frame(probe.like.map)
    plGr <- makeGRangesFromDataFrame(probe.like.map, start.field = "V2", end.field = "V2", seqnames.field = "V1")
    plGr <- BiocGenerics::unique(plGr)
    probe.like.map <- as.data.frame(plGr)
    probe.like.map$Name <- paste0("probe.like_", seq(from=1, to=nrow(probe.like.map)))
    probe.like.map <- subset(probe.like.map, select = c(Name, seqnames, start))
    colnames(probe.like.map) <- c("Name", "Chr", "Position")
    probe.like.map$Chr <- gsub("chr", "", probe.like.map$Chr)
    write.table(probe.like.map, file.path(all.paths[1], "/MapPenn.txt"), col.names = TRUE, row.names = FALSE,
                quote = FALSE, sep="\t")
    
  }else file.copy(map.loc, file.path(all.paths[1], "/MapPenn.txt"))
  
  phen.info <- .loadPhen("Pheno.txt", all.paths[1])
  
  PlinkVer <- .getPLINK(all.paths[2])
  
  if(PlinkVer != "TRUE"){
    stop("PLINK setup failed")}
  
  NewPLINK <- file.path(all.paths[2], dir(path = ".", pattern = "plink-1*"))
  all.paths[2] <- NewPLINK
  
  phen.info$all.paths <- all.paths
  return(phen.info)
}

#' HELPER - Write CNV-genotype in a probe-wise fashion
#' 

.writeProbesCNV <- function(lo, all.samples, genofile, CNVsGr, probes.cnv.gr, n){
  sampleX <- all.samples[[lo]]
  g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id=sampleX))
  CNVSamByState <- split(CNVsGr[values(CNVsGr)$V5 == all.samples[lo]], 
                         as.character(elementMetadata(CNVsGr[values(CNVsGr)$V5 == all.samples[lo]])$V4))
  for(slo in seq_len(length(CNVSamByState))){
    state <- unique(substr(x=as.character(values(CNVSamByState[[slo]])$V4),
                           start=nchar(as.character(values(CNVSamByState[[slo]])$V4)),
                           stop=nchar(as.character(values(CNVSamByState[[slo]])$V4))))
    SamCNVSNP <- values(subsetByOverlaps(probes.cnv.gr, CNVSamByState[[slo]]))$snp.id
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
#' Function to produce the GDS file in a probe-wise fashion for CNV genotypes. The GDS file which is produced
#' also incorporates one phenotype to be analyzed. If several phenotypes are enclosed in the \sQuote{phen.info}
#' object, the user may specify the phenotype to be analyzed with the \sQuote{lo.phe} parameter.  
#'
#'
#' @param phen.info Returned by \code{setupCnvAna}
#' @param freq.cn Minimum frequency where 1 = 100%. Default 0.01 = 1%
#' @param snp.matrix Only FALSE implemented - If TRUE B allele frequencies (BAF) would be used to reconstruct CNV-SNP genotypes
#' @param n.cor Number of cores to be used
#' @param lo.phe The phenotype to be analyzed in the PhenInfo$phenotypesSam data-frame 
#' @param chr.code.name A data-frame with the integer name in the first column and the original name for each chromosome 
#' @return probes.cnv.gr object with information about all probes to be used in the downstream CNV-GWAS
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @seealso (\code{\link{snpgdsCreateGeno}})
#' @examples
#' 
#'  Construct the data-frame with integer and original chromosome names 
#'df <- "27 1A
#' 25 25LG1
#' 29 4A
#' 30 LGE22
#' 31 25LG2
#' 32 LG22
#' 33 Z"
#' 
#' chr.code.name <- read.table(text=df, header=F)
#' 
#' head(chr.code.name)
#'
#' probes.cnv.gr <- prodGdsCnv(phen.info, chr.code.name=chr.code.name)
#' 
#'   @export

prodGdsCnv <- function(phen.info, freq.cn=0.01, snp.matrix=FALSE, n.cor=1, lo.phe=1, chr.code.name=NULL){  
  
  message("initializing ...", appendLF = FALSE)
  
  phenotypesSam <- phen.info$phenotypesSam
  samplesPhen <- phen.info$samplesPhen
  FamID <- phen.info$FamID
  SexIds <- phen.info$SexIds
  all.paths <- phen.info$all.paths
  
  phenotypesSamX <- phenotypesSam[,c(1,(lo.phe+1))]
  phenotypesSamX <- na.omit(phenotypesSamX) ### Exclude samples without phenotypes. i.e NA
  
  ######################  Import CNVs to data-frame
  
  cnvs <- read.table(file.path(all.paths[1], "CNVOut.txt"), sep="", header=F) ### CNV table 
  
  CNVs <- .CheckConvertCNVs(cnvs)
  
  
  ####################### Check if the chromosomes are numeric
  
  chr.names <- CNVs$chr
  chr.names <- gsub("chr", "", chr.names)
  
  if(any(is.na(chr.names))){
    stop("Your chromosome names should be integers. If they are not, make them integers and include the correspondent names 
         in chr.code.name parameter")
  }
  
  #############################
  
  CNVs$start <- as.integer(as.character(CNVs$start))
  CNVs$end <- as.integer(as.character(CNVs$end))
  CNVsGr <- makeGRangesFromDataFrame(CNVs, keep.extra.columns = TRUE)
  CNVsGr <- CNVsGr[values(CNVsGr)$V5 %in% samplesPhen] ### Subset CNVs in phenotyped samples
  
  ######################  Import SNP map to data-frame
  probes <- data.table::fread("MapPenn.txt", header=TRUE, sep="\t")
  probes <- as.data.frame(probes)
  probes$Position <- as.numeric(as.character(probes$Position))
  probes <- probes[complete.cases(probes), ]
  probesGr <- makeGRangesFromDataFrame(probes, seqnames.field = "Chr", start.field = "Position", end.field = "Position",
                                       keep.extra.columns = TRUE)
  
  all.samples <- unique(as.character(values(CNVsGr)$V5))
  
  ###################### Select probes  within CNVs
  probesCNV <- suppressWarnings(as.character(values(subsetByOverlaps(probesGr, 
                                                                     CNVsGr))$Name))
  probesCNV <- unique(unlist(probesCNV))
  probes.cnv.gr <- probesGr[values(probesGr)$Name %in% probesCNV]
  
  counts <- as.integer(countOverlaps(probes.cnv.gr, CNVsGr))
  values(probes.cnv.gr)$freq <- counts
  
  ##### Subset by frequency 
  NumSam <- freq.cn*length(all.samples)
  probes.cnv.gr <- probes.cnv.gr[values(probes.cnv.gr)$freq>=NumSam]
  
  ##### Order the probes
  probes.cnv.gr <- GenomeInfoDb::sortSeqlevels(probes.cnv.gr)
  probes.cnv.gr <- sort(probes.cnv.gr)
  
  values(probes.cnv.gr)$snp.id <- seq(1:length(probes.cnv.gr))
  
  ###################### SNP genotype matrix not available 
  if(snp.matrix == "FALSE"){
    CNVBiMa <- matrix(2, nrow = length(probes.cnv.gr), 
                      ncol = length(all.samples))}
  
  ###################### SNP genotype matrix available - TODO
  if(snp.matrix == "TRUE"){
    ...
  }
  
  ###################### Create a GDS with chr and SNP names to numeric
  SNPRelate::snpgdsCreateGeno("CNV.gds", 
                              genmat = CNVBiMa,
                              sample.id = all.samples, 
                              snp.id = as.character(values(probes.cnv.gr)$snp.id),
                              snp.rs.id = as.character(values(probes.cnv.gr)$Name),
                              snp.chromosome = as.character(seqnames(probes.cnv.gr)),
                              snp.position = start(probes.cnv.gr)
  )
  
  testit(15)
  
  ###################### Replace genotype matrix with CNV genotypes
  (genofile <- SNPRelate::snpgdsOpen("CNV.gds", allow.fork=TRUE, readonly=FALSE))
  
  CNVBiMaCN <- matrix(2, nrow = length(probes.cnv.gr), 
                      ncol = length(all.samples))
  n <- gdsfmt::add.gdsn(genofile, "CNVgenotype", CNVBiMaCN, replace=TRUE)
  read.gdsn(n)
  
  gback <- as.numeric(CNVBiMaCN[,1])
  
  if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
    multicoreParam <-  BiocParallel::MulticoreParam(workers = n.cor)
    suppressMessages(BiocParallel::bplapply(1:length(all.samples), .writeProbesCNV, BPPARAM = multicoreParam, all.samples=all.samples,
                                            genofile=genofile, CNVsGr=CNVsGr, probes.cnv.gr=probes.cnv.gr, n=n))
  }
  
  if(rappdirs:::get_os() == "win"){
    param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
    suppressMessages(BiocParallel::bplapply(1:length(all.samples), .writeProbesCNV, BPPARAM = param, all.samples=all.samples, 
                                            genofile=genofile, CNVsGr=CNVsGr, probes.cnv.gr=probes.cnv.gr, n=n))
  }
  
  testit(15) ## Wait to make sure that gds is writen
  
  message("GDS writen...", appendLF = FALSE)
  
  ###################### Include the phenotype 'lo' in the GDS
  phenotypesSamX <- phenotypesSamX[match(all.samples, phenotypesSamX$samplesPhen),] ### Check order correspondence
  gdsfmt::add.gdsn(genofile, name="phenotype", val=phenotypesSamX, replace=TRUE)
  FamID <- FamID[match(all.samples, FamID$samplesPhen),] ### Check order correspondence
  gdsfmt::add.gdsn(genofile, name="FamID", val=FamID[,2], replace=TRUE)
  FamID <- SexIds[match(all.samples, SexIds$samplesPhen),] ### Check order correspondence
  gdsfmt::add.gdsn(genofile, name="Sex", val=SexIds[,2], replace=TRUE)
  
  ## Include chr names if needed
  if(!is.null(chr.code.name)){
    gdsfmt::add.gdsn(genofile, name="Chr.names", val=chr.code.name, replace=TRUE)  
  }
  
  SNPRelate::snpgdsClose(genofile)
  ########################################### END #######################################
  message("done", appendLF = FALSE)
  return(probes.cnv.gr)}

#' HELPER - Internal function of parallel implementation to produce the .gvar file requested for the CNV-GWAS in PLINK
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

#' HELPER - Produce PLINK probe map in disk
#' 

.prodPLINKmap <- function(all.paths){
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  Pmap <- as.data.frame(cbind(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                              gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
                              0, gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position"))))
  
  colnames(Pmap) <- c("Chr", "NAME", "GD", "Position")
  write.table(Pmap, file.path(all.paths[2], "mydata.map"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  SNPRelate::snpgdsClose(genofile)
}

#' HELPER - Produce the PLINK .gvar file in disk
#' 

.prodPLINKgvar <- function(all.paths){
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
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
  gentype.all1 <- as.data.frame(gentype.all)
  
  write.table(gentype.all1, file.path(all.paths[2],"mydata.gvar"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  SNPRelate::snpgdsClose(genofile)
}


#' HELPER - Produce the PLINK fam file in disk
#' 

.prodPLINKfam <- function(all.paths){
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=T, readonly=FALSE)) ## read GDS
  SamGen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")) ## Extract samples
  SNPs <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")) ## Extract SNPs
  FAMID <- as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "FamID"))) ## Extract family
  Sex <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Sex")) ## Extract sex
  Phen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")) ## Extract phenotype
  
  famdf <- as.data.frame(cbind(FAMID, SamGen, SamGen, "NC", Sex, Phen[,2]))
  colnames(famdf) <- c("FAID", "IID", "PID", "MID", "Sex", "Phenotype")
  
  write.table(famdf, file.path(all.paths[2]), "mydata.fam", sep = "\t", col.names = TRUE, row.names=FALSE, quote=FALSE)
  SNPRelate::snpgdsClose(genofile)
}

#' HELPER - Run PLINK
#' 

.runPLINK <- function(all.paths){
  
  plinkPath <- paste0("\'", all.paths[2], "/plink\'")
  system(paste(plinkPath, "--gfile mydata --noweb"), wait=TRUE)
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=T, readonly=FALSE)) ## read GDS
  sumphen <- paste0("plink.gvar.summary.",names(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")))[[2]])
  
  file.rename(file.path(all.paths[2], "plink.gvar.summary"), file.path(all.paths[2],sumphen))
  SNPRelate::snpgdsClose(genofile)
  
}

#' HELPER - Produce CNV segments
#' 

.prodCNVseg <- function(all.paths, probes.cnv.gr){
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=T, readonly=FALSE)) ## read GDS
  
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
  
  values(probes.cnv.gr)$starts.seg <- unlist(all.segs)
  
  pstarts <- probes.cnv.gr[values(probes.cnv.gr)$starts.seg == "start"]
  pstarts.split = split(pstarts, GenomeInfoDb::seqnames(pstarts))
  
  all.segs.gr <- GRangesList()
  for(chx in 1:length(pstarts.split)){
    pstartsx <- pstarts.split[[chx]]
    if(length(pstartsx)==0){
      all.segs.gr[[chx]] <- pstartsx
      next
    }
    all.segs.grX <- GRangesList()
    
    for(st in 1:length(pstartsx)){
      prx <- pstartsx[st]
      if(st < length(pstartsx)){
        LimPrx <- BiocGenerics::start(pstartsx[st+1])-1}
      if(st == length(pstartsx)){
        LimPrx <-  max(BiocGenerics::start(probes.cnv.gr[GenomeInfoDb::seqnames(probes.cnv.gr) == GenomeInfoDb::seqnames(prx)]))
      }
      seqLar <- GRanges(as.character(GenomeInfoDb::seqnames(prx)), IRanges::IRanges(BiocGenerics::start(prx), LimPrx))
      seqPrbs <- subsetByOverlaps(probes.cnv.gr,seqLar)
      seqx <- GRanges(as.character(GenomeInfoDb::seqnames(prx)), IRanges::IRanges(BiocGenerics::start(prx), BiocGenerics::start(seqPrbs[length(seqPrbs)])))
      all.segs.grX[[st]] <- seqx
    }
    all.segs.gr[[chx]] <- unlist(all.segs.grX)
  }
  
  snpgdsClose(genofile)
  
  return(all.segs.gr)
  
}


#' HELPER - Associate probes with CNV segments and draw p-values
#' 

.assoPrCNV  <- function(all.paths, all.segs.gr){
  
  segs.pvalue.gr <- unlist(all.segs.gr)
  values(segs.pvalue.gr)$SegName <- seq(from=1, to=length(segs.pvalue.gr), by=1)
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  
  sum.file <- paste0("plink.gvar.summary.",names(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")))[[2]])
  results <- read.table(file.path(all.paths[2], sum.file), header=TRUE) 
  mydata.map <- read.table("mydata.map", header=TRUE)
  resultsp <- merge(results, mydata.map, by.x="NAME", by.y="NAME", all.x=TRUE, all.y=FALSE)
  
  #### Choose the method - TODO - Change if SNPMatrix implemented - Allow to get other p-values
  resultsp <- resultsp[ grep("P\\(CNP\\)", resultsp$FIELD),]
  resultsp$VALUE <- as.numeric(as.character(resultsp$VALUE))
  resultsp <- resultsp[order(resultsp$VALUE, decreasing=FALSE),]
  resultspGr <- makeGRangesFromDataFrame(resultsp, seqnames.field = "Chr",
                                         start.field = "Position", end.field = "Position",
                                         keep.extra.columns = TRUE)
  
  ######## Assign the lowest p-value to the CNV segment
  
  values.all <- vector("list", length(segs.pvalue.gr))
  names.all <- vector("list", length(segs.pvalue.gr))
  for(se in 1:length(segs.pvalue.gr)){
    values.all[[se]] <- min(values(subsetByOverlaps(resultspGr, segs.pvalue.gr[se]))$VALUE)
    Prbx <- subsetByOverlaps(resultspGr, segs.pvalue.gr[se])
    names.all[[se]] <- as.character(Prbx[values(Prbx)$VALUE %in% min(values(subsetByOverlaps(resultspGr, segs.pvalue.gr[se]))$VALUE)]$NAME[[1]])
  }
  
  values(segs.pvalue.gr)$MinPvalue <- unlist(values.all)
  values(segs.pvalue.gr)$NameProbe <- unlist(names.all)
  
  values(segs.pvalue.gr)$MinPvalueAdjusted <- stats::p.adjust(values(segs.pvalue.gr)$MinPvalue, method =method.m.test, n = length(segs.pvalue.gr))
  segs.pvalue.gr$Phenotype <- names(phenotypesSamX)[[2]]
  SNPRelate::snpgdsClose(genofile)
  
  return(segs.pvalue.gr)
  
}

#' Run the CNV-GWAS 
#' 
#' (i) Produces the requested inputs for a PLINK analysis and (ii) run a 
#' CNV-GWAS analysis using the GDS wrote by the last \code{\link{prodGdsCnv}}
#' run.
#'
#' @param phen.info Returned by \code{\link{setupCnvAna}}
#' @param n.cor Number of cores to be used
#' @param min.sim Minimum CNV genotype distribuition similarity among subsequent probes. Default is 0.95 = 95%
#' @param freq.cn Minimum frequency where 1 = 100%. Default 0.01 = 1%
#' @param snp.matrix Only FALSE implemented - If TRUE B allele frequencies (BAF) would be used to reconstruct CNV-SNP genotypes
#' @param method.m.test Correction for multiple tests to be used. FDR is default, see \code{\link(p.adjust)} for
#' other methods.
#' @param lo.phe The phenotype to be analyzed in the PhenInfo$phenotypesSam data-frame 
#' @param chr.code.name A data-frame with the integer name in the first column and the original name for each chromosome  
#' @return The CNV segments and the representative probes and their respective p-value
#' @examples
#' @export

cnvGWAS <- function(phen.info, n.cor, min.sim = 0.95, freq.cn = 0.01, snp.matrix=FALSE, method.m.test = "fdr", lo.phe=1, chr.code.name=NULL){
  phenotypesSam <- phen.info$phenotypesSam
  phenotypesSamX <- phenotypesSam[,c(1,(lo.phe+1))]
  phenotypesSamX <- na.omit(phenotypesSamX) ### Exclude samples without phenotypes. i.e NA
  all.paths <- phen.info$all.paths
  
  #################################### Produce the GDS for a given phenotype
  
  probes.cnv.gr <- prodGdsCnv(phen.info, freq.cn, snp.matrix, n.cor, lo.phe, chr.code.name)
  
  ##################################### Produce PLINK map ##################
  
  .prodPLINKmap(all.paths)
  
  ########################################### END #######################################
  
  ##################################### Produce gvar to use as PLINK input #########
  
  .prodPLINKgvar(all.paths)
  
  ########################################### END #######################################
  
  ##################################### Produce fam (phenotype) to use as PLINK input #############
  
  .prodPLINKfam(all.paths)
  
  ########################################### END #######################################
  
  ##################################### Run PLINK #############
  
  .runPLINK(all.paths)
  
  ########################################### END #######################################
  
  ##################################### Produce CNV segments #############
  
  all.segs.gr <- .prodCNVseg(all.paths, probes.cnv.gr)
  
  ########################################### END #######################################
  
  ##################################### Associate SNPs with CNV segments #############
  
  segs.pvalue.gr <- .assoPrCNV(all.paths, all.segs.gr)
  
  ######################################################## END ###########################################
  
  #################################### Plot the QQ-plot of the analysis
  
  pdf(file.path(all.paths[3], paste0(unique(values(segs.pvalue.gr)$Phenotype),"QQ-PLOT.pdf")))
  qqunif.plot(values(segs.pvalue.gr)$MinPvalue, auto.key=list(corner=c(.95,.05)))
  
  dev.off()
  ######################################################## END ###########################################
  
  #################################### Reconvert the chrs to original names if applicable
  
  if(!is.null(chr.code.name)){
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=T, readonly=FALSE)) ## read GDS
  chr.names.df < gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Chr.names"))
  
  segs.pvalue <- data.frame(segs.pvalue.gr)
  segs.pvalue <- segs.pvalue[ , -which(names(segs.pvalue) %in% c("strand"))]
  
  for(lopN in len_seq(nrow(chr.names.df))){
  segs.pvalue$seqnames <- gsub(chr.names.df[lopN,1],chr.names.df[lopN,2], segs.pvalue$seqnames)
  }
  
  segs.pvalue.gr <- makeGRangesFromDataFrame(segs.pvalue, keep.extra.columns = TRUE)
  
  SNPRelate::snpgdsClose(genofile)
  }
  
  ######################################################## END ########################################
  
  
  #################################### Plot the manhattan ######################
  
  #TODO
  
  ######################################################## END ###########################################
  

  return(segs.pvalue.gr)
}

