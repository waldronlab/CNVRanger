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
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree
#' @param version PLINK version. Only 1.07 implemented
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
#' @param pops.names Indicate the name of the populations. Only used if more than one population exists
#' @return List \sQuote{phen.info} with \sQuote{samplesPhen}, \sQuote{phenotypes}, \sQuote{phenotypesdf}, 
#' \sQuote{phenotypesSam}, \sQuote{FamID} and \sQuote{SexIds} and \sQuote{pops.names} (if more than one population)
#' @examples
#' phen.info <- loadPhen("/home/.../Phen.txt")
#' sapply(phen.info, class)

.loadPhen <- function(file.nam, pheno.path, pops.names=NULL){
  
  samplesPhen.all <- NULL
  phenotypes.all <- NULL
  phenotypesdf.all <- NULL
  phenotypesSam.all <- NULL
  FamID.all <- NULL
  SexIds.all <- NULL
  
  for(npop in seq_len(length(file.nam))){
    dfN <- data.table::fread(file.path(pheno.path, file.nam[npop]))
    if(all(names(dfN) == "INEXISTENT")){
      ## Extract the sample names from CNV file instead
      cnv.file <- paste0("CNVOutPop", npop,".txt")
      
      cnvs <- read.table(file.path(pheno.path, cnv.file), sep="", header=F) ### - adapt to multiple populations 
      CNVs <- .CheckConvertCNVs(cnvs)
      all.samples <- as.character(CNVs$V5)
      all.samples <- unique(all.samples)
      samplesPhen <- all.samples
      phenotypes <- phenotypes.all[[1]]
      phenotypesdf <- data.frame(rep("NA", length(all.samples)))
      colnames(phenotypesdf) <- phenotypes
      phenotypesSam <- cbind(all.samples, phenotypesdf)
      colnames(phenotypesSam)[1] <- "samplesPhen"
      FamID <- cbind(samplesPhen=as.character(phenotypesSam$samplesPhen), V2=rep("NA", length(all.samples)))
      FamID <- as.data.frame(FamID)
      SexIds <- cbind(samplesPhen=as.character(phenotypesSam$samplesPhen), V2=rep("NA", length(all.samples)))
      SexIds <- as.data.frame(SexIds)
      
    }else{
      all <- data.table::fread(file.path(pheno.path, file.nam[npop]), header=TRUE, sep="\t") ### Import phenotypes
      all <- as.data.frame(all)
      
      if(!all(colnames(all)[1:3]==c("sample.id", "fam", "sex"))){
        stop("Unexpected first three column names")}
      
      ### Include samples with CNV but without phenotypes
      cnv.file <- paste0("CNVOutPop", npop,".txt")
      
      cnvs <- read.table(file.path(pheno.path, cnv.file), sep="", header=F) 
      CNVs <- .CheckConvertCNVs(cnvs)
      all.samples <- as.character(CNVs$V5)
      all.samples <- unique(all.samples)
      
      pheno.na <- cbind(as.data.frame(all.samples), "fam"=NA, "sex"=NA) 
      colnames(pheno.na)[1] <- "sample.id"
      pheno.na$sample.id <- as.character(pheno.na$sample.id)
      pheno.na <- pheno.na[-which(pheno.na$sample.id %in% all$sample.id),]
      
      all <- plyr::rbind.fill(all, pheno.na)
      all[is.na(all)] <- -9
      
      ### Produce phen.info objects
      samplesPhen <- unique(all$sample.id) ## Greb sample names
      phenotypes <- names(all[, 4:ncol(all), drop=FALSE]) ## Greb phenotype names
      
      phenotypesdf <- all[, 4:ncol(all), drop = FALSE]
      phenotypesdf <- as.data.frame(phenotypesdf) ### Greb phenotypes values
      #phenotypesdf <- as.data.frame(phenotypesdf[!duplicated(phenotypesdf),,drop = FALSE])
      
      phenotypesSam <- data.frame(samplesPhen, phenotypesdf)
      phenotypesSam$samplesPhen <- as.character(phenotypesSam$samplesPhen)
      
      FamID <- data.frame(samplesPhen, "V2"=all$fam[as.numeric(rownames(phenotypesSam))])
      
      SexIds <- data.frame(samplesPhen,"V2"=all$sex[as.numeric(rownames(phenotypesSam))])
    }
    
    samplesPhen.all[[npop]] <- samplesPhen
    phenotypes.all[[npop]] <- phenotypes
    phenotypesdf.all[[npop]] <- phenotypesdf
    phenotypesSam.all[[npop]] <- phenotypesSam
    FamID.all[[npop]] <- FamID
    SexIds.all[[npop]] <- SexIds
  }
  
  samplesPhen <- unlist(samplesPhen.all)
  phenotypes <- unique(unlist(phenotypes.all))
  phenotypesdf <- data.table::rbindlist(phenotypesdf.all)
  phenotypesdf <- as.data.frame(phenotypesdf)
  phenotypesSam <- data.table::rbindlist(phenotypesSam.all)
  phenotypesSam <- as.data.frame(phenotypesSam)
  FamID <- data.table::rbindlist(FamID.all)
  FamID <- as.data.frame(FamID)
  SexIds <- data.table::rbindlist(SexIds.all)
  SexIds <- as.data.frame(SexIds)
  
  if(!is.null(pops.names)){
    
    all.pop.names <- NULL
    for(npop in seq_len(length(file.nam))){
      all.pop.names[[npop]] <- rep(pops.names[npop], length(samplesPhen.all[[npop]]))
    }
    
    pops.names <- unlist(all.pop.names)
    phen.info <- list("samplesPhen" = samplesPhen, "phenotypes" = phenotypes, "phenotypesdf" = phenotypesdf,
                      "phenotypesSam" = phenotypesSam, "FamID" = FamID, "SexIds" = SexIds, "pops.names" = pops.names)
  }
  else{
    phen.info <- list("samplesPhen" = samplesPhen, "phenotypes" = phenotypes, "phenotypesdf" = phenotypesdf,
                      "phenotypesSam" = phenotypesSam, "FamID" = FamID, "SexIds" = SexIds)}
  
  return(phen.info)
}

#' HELPER - Check and convert CNV input
#' @param cnvs Data-frame with the CNVs to be analyzed. From (i) PennCNV, (ii) SNP-chip general format 
#' or (iii) sequencing general format
#' 

.CheckConvertCNVs <- function(cnvs){
  
  ## If PennCNV file
  if(unique(sub("^([[:alpha:]]*).*", "\\1", cnvs$V2))[[1]]=="numsnp"){
    cnvs <- as.data.frame(cnvs)
    df1 <- reshape2::colsplit(cnvs$V1, ":", c("chr", "loc"))
    df2 <- reshape2::colsplit(df1$loc, "-", c("start", "end"))
    df1 <- cbind(df1[,-2, drop=FALSE], df2)
    
    cnvs <- cbind(df1, cnvs)
    cnvs$V1 <- gsub("chr", "", cnvs$V1)
    colnames(cnvs)[1:3] <- c("chr", "start", "end")
    CNVs <- cnvs
    
    ## If general format
  } 
  else if(all(as.character(as.matrix(cnvs[1,]))==c("chr", "start", "end","sample.id", "state", "num.snps", "start.probe", "end.probe"))){ 
    the.names <- as.character(as.matrix(cnvs[1,]))
    cnvs <- cnvs[-1,]
    colnames(cnvs) <- the.names
    cnvs <- as.data.frame(cnvs)
    cnvs$start <- as.numeric(cnvs$start)
    cnvs$end <- as.numeric(cnvs$end)
    ## Convert to PennCNV format
    cnvs$V1 <- paste0(cnvs$chr, ":",cnvs$start, "-", cnvs$end)
    cnvs$length <- (cnvs$end - cnvs$start)+1
    cnvs$length <- paste0("length=", cnvs$length)
    
    cnvs$num.snps <- paste0("numsnp=", cnvs$length)
    
    cnvs$start.probe <- paste0("startSNP=", cnvs$start.probe)
    cnvs$end.probe <- paste0("endSNP=", cnvs$end.probe)
    
    cnvs <- subset(cnvs, select =c(chr, start, end, V1, num.snps, length, state, sample.id, start.probe, end.probe))
    colnames(cnvs) <- c("chr", "start", "end", "V1", "V2", "V3", "V4", "V5", "V6", "V7")
    CNVs <- cnvs
    
    ## If sequencing info
  } 
  else if(all(as.character(as.matrix(cnvs[1,]))==c("chr", "start", "end", "sample.id", "state"))){
    cnvs <- as.data.frame(cnvs)
    ## Convert to PennCNV 
    
    ## TODO 
    
    
  } else stop("Unexpected CNV input format - is it tab delimited?")
  
  return(CNVs)
}

#' HELPER - Load multiple CNV inputs
#' @param cnv.file.all CNV file name
#' @param cnv.path CNV path name

.loadToMergeCNV <- function(cnv.file.all, cnv.path){
  ## Load and merge CNV files from multiple populations
  data.table::fread(file.path(cnv.path, cnv.file.all), header=TRUE)
}

#' Setup the folders and files to run CNV-GWAS analysis
#'
#' This function creates the (i) necessary folders in disk to perform downstream analysis on CNV genome-wide
#' association and (ii) import the necessary input files (i.e phenotypes, probe map and CNV list) from other 
#' locations in disk. 
#' 
#' The user can import several phenotypes at once. All information will be stored in the list returned by this function. 
#' User should be aware although several phenotypes can be imported, the \code{\link{cnvGWAS} or \code{\link{prodGdsCnv} functions
#' will handle only one phenotype per run. 
#' 
#' @param name String with a project code or name
#' @param phen.loc Path/paths to the tab separated text file cointaining phenotype and sample info. Using more than one population, if a given population
#' don't have phenotypes include the string "INEXISTENT" instead.
#' @param cnv.out.loc Path/paths to the CNV analysis output (i.e. PennCNV output, SNP-chip general format or sequencing general format)
#' @param map.loc Path to the probe map (e.g. used in PennCNV analysis). Column names containing probe name, 
#' chromosome and coordinate must be named as: Name, Chr and Position. Tab delimited.
#' @param folder Choose manually the project folder (i.e. path as the root folder). Otherwise, user-specific data dir
#' will be used automatically.  
#' @param pops.names Indicate the name of the populations, if using more than one
#' @return List \sQuote{phen.info} with \sQuote{samplesPhen}, \sQuote{phenotypes}, \sQuote{phenotypesdf}, 
#' \sQuote{phenotypesSam}, \sQuote{FamID}, \sQuote{SexIds}, \sQuote{pops.names} (if more than one population) and \sQuote{all.paths}
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples
#' 
#' folder.package <- system.file('extdata', package = 'cnvAnalyzeR')
#' 
#' ## One population
#' phen.info <- setupCnvGWAS(name="Example", phen.loc=file.path(folder.package, "Pheno.txt"), 
#' cnv.out.loc <- file.path(folder.package, "CNVOut.txt"), 
#' map.loc <- file.path(folder.package, "MapPenn.txt"))
#' 
#' ## Multiple populations
#' phen.loc <- c(".../PhenoPop1.txt", ".../PhenoPop2.txt")
#' map.loc <- ".../Map.txt"
#' cnv.out.loc <- c(".../CNVOutputPop1.txt", ".../CNVOutputPop2.txt")
#' phen.info <- setupCnvGWAS(name, phen.loc, cnv.out.loc, map.loc)
#'
#' ## Multiple populations where one is lacking phenotypes - always phenotyped populations first
#' phen.loc <- c(".../PhenoPop1.txt", "INEXISTENT")
#' map.loc <- ".../Map.txt"
#' cnv.out.loc <- c(".../CNVOutputPop1.txt", ".../CNVOutputPop2.txt")
#' phen.info <- setupCnvGWAS(name, phen.loc, cnv.out.loc, map.loc)
#' 
#' @export

setupCnvGWAS <- function(name, phen.loc, cnv.out.loc, map.loc=NULL, folder=NULL, pops.names=NULL){
  
  ## Create the folder structure for all subsequent analysis
  all.paths <- .createFolderTree(name, folder)
  
  ## Only one population
  if(length(phen.loc) == 1 && length(cnv.out.loc) == 1){
    ## Import the phenotype and sample info from external folder
    file.copy(phen.loc, file.path(all.paths[1], "/PhenoPop1.txt"), overwrite = TRUE)
    ## Import the PennCNV output from external folder
    file.copy(cnv.out.loc, file.path(all.paths[1], "/CNVOutPop1.txt"),  overwrite = TRUE)
    pheno.file <- "Pheno.txt"
    cnv.file <- "CNVOut.txt"
  }
  
  ## Multiple populations
  if(length(phen.loc) != length(cnv.out.loc)){
    stop("phen.loc and cnv.out.loc should have the same length. Use the string INEXISTENT if phenotypes are missing")}
  
  if(length(phen.loc) > 1 && length(cnv.out.loc) > 1){
    
    pheno.file.all <- NULL  
    cnv.file.all <- NULL
    
    for(npop in seq_len(length(phen.loc))){
      pheno.file <- paste0("/PhenoPop", npop,".txt")
      pheno.file.all[[npop]] <- pheno.file
      cnv.file <- paste0("/CNVOutPop", npop,".txt")
      cnv.file.all[[npop]] <- cnv.file
      
      if(phen.loc[npop] == "INEXISTENT"){
        write.table(c("INEXISTENT"), file.path(all.paths[1], pheno.file), quote=FALSE, 
                    col.names = FALSE, row.names=FALSE)
      }else{
        file.copy(phen.loc[npop], file.path(all.paths[1], pheno.file), overwrite = TRUE)}
      file.copy(cnv.out.loc[npop], file.path(all.paths[1], cnv.file), overwrite = TRUE)
    }
  }
  
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
    
  }
  else{
    if(length(map.loc)>1){
      stop("The map should be unique. If multiple populations with different probe map use map.loc=NULL")  
    }
    file.copy(map.loc, file.path(all.paths[1], "/MapPenn.txt"), overwrite = TRUE)}
  
  ## Write phenotype names and merge CNV file for multiple populations
  if(length(phen.loc) > 1 && length(cnv.out.loc) > 1){
    pheno.file<-unlist(pheno.file.all)
    all.cnvs <- lapply(cnv.file.all,.loadToMergeCNV, cnv.path=all.paths[1])
    all.cnvs <- data.table::rbindlist(all.cnvs)
    all.cnvs <- as.data.frame(all.cnvs)
    write.table(all.cnvs, file.path(all.paths[1], "CNVOut.txt"), sep="\t",
                col.names=TRUE, row.names = FALSE)
  }
  
  phen.info <- .loadPhen(pheno.file, all.paths[1], pops.names=pops.names)
  
  PlinkVer <- .getPLINK(all.paths[2])
  
  if(PlinkVer != "TRUE"){
    stop("PLINK setup failed")}
  
  NewPLINK <- file.path(all.paths[2], dir(path = all.paths[2], pattern = "plink-1*"))
  all.paths[2] <- NewPLINK
  
  phen.info$all.paths <- all.paths
  return(phen.info)
}

#' HELPER - Write CNV-genotype in a dosage-like fashion
#' 

.writeProbesCNV <- function(lo, all.samples, genofile, CNVsGr, probes.cnv.gr, n){
  sampleX <- all.samples[[lo]]
  g <- as.numeric(suppressMessages(SNPRelate::snpgdsGetGeno(genofile, sample.id=sampleX)))
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

#' HELPER - Write CNV-genotype in a SNP-like fashion
#' @param lo loop number
#' @param genofile loaded gds file
#' @param n is the object from \code{gdsfmt::add.gdsn}
#' @param ranges.gr is the CNVs of each sample

.replaceSNPtoCNV <- function(lo, all.samples, genofile, CNVsGr, probes.cnv.gr, n){
  sampleX <- all.samples[[lo]]
  g <- as.numeric(suppressMessages(SNPRelate::snpgdsGetGeno(genofile, sample.id=sampleX)))
  
  #2n 
  g[seq_len(length(g))] <- "1"
  
  # 0n 
  CNVsGr0n <- CNVsGr[which(values(CNVsGr)$sample == sampleX)]
  CNVsGr0n <- CNVsGr0n[which(values(CNVsGr0n)$type == "state1,cn=0")]
  probes.cnv.gr.0n <- subsetByOverlaps(probes.cnv.gr, CNVsGr0n)
  
  g[values(probes.cnv.gr.0n)$tag.snp] <- "0"
  
  # 1n
  CNVsGr1n <- CNVsGr[which(values(CNVsGr)$sample == sampleX)]
  CNVsGr1n <- CNVsGr1n[which(values(CNVsGr1n)$type == "state2,cn=1")]
  probes.cnv.gr.1n <- subsetByOverlaps(probes.cnv.gr, CNVsGr1n)
  
  g[values(probes.cnv.gr.1n)$tag.snp] <- "0"
  
  # 3n
  CNVsGr3n <- CNVsGr[which(values(CNVsGr)$sample == sampleX)]
  CNVsGr3n <- CNVsGr3n[which(values(CNVsGr3n)$type == "state5,cn=3")]
  probes.cnv.gr.3n <- subsetByOverlaps(probes.cnv.gr, CNVsGr3n)
  
  g[values(probes.cnv.gr.3n)$tag.snp] <- "2"
  
  # 4n
  CNVsGr4n <- CNVsGr[which(values(CNVsGr)$sample == sampleX)]
  CNVsGr4n <- CNVsGr4n[which(values(CNVsGr4n)$type == "state6,cn=4")]
  probes.cnv.gr.4n <- subsetByOverlaps(probes.cnv.gr, CNVsGr4n)
  
  g[values(probes.cnv.gr.4n)$tag.snp] <- "2"
  
  # 2n
  #g[c(-(values(probes.cnv.gr.0n)$tag.snp),-(values(probes.cnv.gr.1n)$tag.snp),-(values(probes.cnv.gr.3n)$tag.snp),-(values(probes.cnv.gr.4n)$tag.snp))] <- "1"
  
  g <- as.numeric(g)
  
  ### Replace the genotype in the gds file
  gdsfmt::write.gdsn(n, g, start=c(1,lo), count=c(length(g),1))
  
}

#' HELPER - Write CNV-genotype in SNP-like fashion in biallelic loci 
#' @param lo loop number, based on the number of SNPs per turn
#' @param genofile loaded gds file
#' @param n is the object from \code{gdsfmt::add.gdsn}
#' @param chunk Intervals of SNP to import and convert
#' @param coding.translate For 'CNVgenotypeSNPlike'. If NULL or unrecognized string use only biallelic CNVs. If 'all' code multiallelic CNVs as 0 
#' for loss; 1 for 2n and 2 for gain.  

.recodeCNVgenotype <- function(lo, genofile, all.samples, n, chunk, coding.translate){
  
  count.end <- (chunk[lo+1])-chunk[lo]
  
  if(length(chunk)-1 != lo){
    CNVgenoX <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(chunk[lo],1), count=c(count.end, length(all.samples))))
  }else{
    CNVgenoX <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(chunk[lo],1), count=c(count.end+1, length(all.samples))))
  }
  
  ### Put four copies as max
  CNVgenoX[CNVgenoX>4] <- 4
  CNVgenoX <- as.data.frame(CNVgenoX)
  
  CNVgenoX$zn <- 0
  CNVgenoX$on <- 1
  CNVgenoX$tn <- 2
  CNVgenoX$thn <- 3
  CNVgenoX$fn <- 4
  
  ### Count genotypes
  genos <- apply(CNVgenoX, 1,table)
  genos <- genos-1
  genos <- t(genos)
  
  #### Recode genotypes
  indexLoss <-apply(genos[,1:2], 1, sum)>0 & apply(genos[,4:5], 1, sum)==0
  indexGain <-apply(genos[,1:2], 1, sum)==0 & apply(genos[,4:5], 1, sum)>0
  indexExclude <-apply(genos[,1:2], 1, sum)>0 & apply(genos[,4:5], 1, sum)>0
  
  ## Exclude index in the matrix
  CNVgenoX  <- CNVgenoX[, -((ncol(CNVgenoX)-4):ncol(CNVgenoX))]
  
  ### To character
  CNVgenoX <- t(apply(CNVgenoX,1, paste0, "n"))
  
  ### Replace Gain
  CNVgenoX[indexGain,] <- gsub("2n","0", CNVgenoX[indexGain,])
  CNVgenoX[indexGain,] <- gsub("3n","1", CNVgenoX[indexGain,])
  CNVgenoX[indexGain,] <- gsub("4n","2", CNVgenoX[indexGain,])
  
  ### Replace Loss
  CNVgenoX[indexLoss,] <- gsub("2n","0", CNVgenoX[indexLoss,])
  CNVgenoX[indexLoss,] <- gsub("1n","1", CNVgenoX[indexLoss,])
  CNVgenoX[indexLoss,] <- gsub("0n","2", CNVgenoX[indexLoss,])
  
  
  ## Non bi-allelic CNVs = no genotype
  if(is.null(coding.translate)){
    coding.translate <- "biallelic"
  }
  if(coding.translate == "all"){
    CNVgenoX[indexExclude,] <- gsub("4n","2", CNVgenoX[indexExclude,])
    CNVgenoX[indexExclude,] <- gsub("3n","2", CNVgenoX[indexExclude,])
    CNVgenoX[indexExclude,] <- gsub("2n","-1", CNVgenoX[indexExclude,])
    CNVgenoX[indexExclude,] <- gsub("1n","0", CNVgenoX[indexExclude,])
    CNVgenoX[indexExclude,] <- gsub("0n","0", CNVgenoX[indexExclude,])
  }else{
    CNVgenoX[indexExclude,] <- -1  
  }
  
  CNVgenoX <- apply(CNVgenoX, 2, as.numeric)
  
  ### Replace the genotype in the gds file
  if(length(chunk)-1 != lo){
    gdsfmt::write.gdsn(n, CNVgenoX, start=c(chunk[lo],1), count=c(count.end,length(all.samples)))
  }else{
    gdsfmt::write.gdsn(n, CNVgenoX, start=c(chunk[lo],1), count=c(count.end+1,length(all.samples)))
  }
}

#' HELPER Wait x seconds
#' # from https://stat.ethz.ch/R-manual/R-devel/library/base/html/Sys.sleep.html

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
#' @param genotype.nodes Expression data type. Nodes with CNV genotypes to be produced in the gds file. 
#' Use 'CNVGenotype' for dosage-like genotypes (i.e. from 0 to Inf). Use 'CNVgenotypeSNPlike' alonside for SNP-like 
#' CNV genotype in a separated node (i.e.  '0, 1, 2, 3, 4' as '0/0, 0/1, 1/1, 1/2, 2/2').
#' @param coding.translate For 'CNVgenotypeSNPlike'. If NULL or unrecognized string use only biallelic CNVs. If 'all' code multiallelic CNVs as 0 
#' for loss; 1 for 2n and 2 for gain. 
#' @param lrr logical. Import LRR/BAF values from external files
#' @param path.files Folder containing the input CNV files used for the CNV calling (i.e. one text file with 5 collumns 
#' for each sample). Columns should contain (i) probe name, (ii) Chromosome, (iii) Position, (iv) LRR and (v) BAF. Only if lrr = TRUE
#' @param list.of.files Data-frame with two columns where the (i) is the file name with signals and (ii) is the 
#' correspondent name of the sample in the gds file. Only if lrr = TRUE
#' @return probes.cnv.gr object with information about all probes to be used in the downstream CNV-GWAS
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples
#' 
#' # Load phenotype-CNV information
#' 
#' folder.package <- system.file('extdata', package = 'cnvAnalyzeR')
#' 
#' phen.info <- setupCnvGWAS(name="Example", phen.loc=file.path(folder.package, "Pheno.txt"), 
#' cnv.out.loc <- file.path(folder.package, "CNVOut.txt"), 
#' map.loc <- file.path(folder.package, "MapPenn.txt"))
#' 
#' # Construct the data-frame with integer and original chromosome names 
#'  
#' # Define chr correspondence to numeric, if necessary
#' df <- "16 1A
#' 25 4A
#' 29 25LG1
#' 30 25LG2
#' 31 LGE22"
#' 
#' chr.code.name <- read.table(text=df, header=F)
#' 
#'
#' probes.cnv.gr <- prodGdsCnv(phen.info, chr.code.name=chr.code.name)
#' 
#'@export

prodGdsCnv <- function(phen.info, freq.cn=0.01, snp.matrix=FALSE, n.cor=1, lo.phe=1, 
                       #chr.code.name=NULL, genotype.nodes="CNVGenotype", chunk.len =1000000,
                       chr.code.name=NULL, genotype.nodes="CNVGenotype", coding.translate=NULL, 
                       #lrr=FALSE, path.files=NULL, list.of.files=NULL){  
                       path.files=NULL, list.of.files=NULL){  
  
  message("initializing ...", appendLF = FALSE)
  
  phenotypesSam <- phen.info$phenotypesSam
  samplesPhen <- phen.info$samplesPhen
  FamID <- phen.info$FamID
  SexIds <- phen.info$SexIds
  all.paths <- phen.info$all.paths
  
  phenotypesSamX <- phenotypesSam[,c(1,(lo.phe+1))]
  #phenotypesSamX <- na.omit(phenotypesSamX) ### Exclude samples without phenotypes. i.e NA
  
  ######################  Import CNVs to data-frame
  
  cnvs <- data.table::fread(file.path(all.paths[1], "CNVOut.txt"), sep="\t", header=FALSE) ### CNV table 
  CNVs <- .CheckConvertCNVs(cnvs)
  
  ####################### Check if the chromosomes are numeric
  
  chr.names <- CNVs$chr
  chr.names <- gsub("chr", "", chr.names)
  
  if(any(is.na(chr.names))){
    stop("Chromosome names should be integers. If they are not, make them integers and include the correspondent names 
         in chr.code.name parameter")
  }
  
  #############################
  
  CNVs$start <- as.integer(as.character(CNVs$start))
  CNVs$end <- as.integer(as.character(CNVs$end))
  CNVsGr <- makeGRangesFromDataFrame(CNVs, keep.extra.columns = TRUE)
  CNVsGr <- CNVsGr[values(CNVsGr)$V5 %in% samplesPhen] ### Subset CNVs in phenotyped samples
  
  ######################  Import SNP map to data-frame
  probes <- data.table::fread(file.path(all.paths[1], "MapPenn.txt"), header=TRUE, sep="\t")
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
    ... ## CHECK IF WE HAVE THE SAME SAMPLES - INCLUDE NA FOR SAMPLES WITH ONLY ONE GENOTYPE TYPE? (i.e CNV or SNP)
  }
  
  ###################### Create a GDS with chr and SNP names to numeric
  SNPRelate::snpgdsCreateGeno(file.path(all.paths[1], "CNV.gds"), 
                              genmat = CNVBiMa,
                              sample.id = all.samples, 
                              snp.id = as.character(values(probes.cnv.gr)$snp.id),
                              snp.rs.id = as.character(values(probes.cnv.gr)$Name),
                              snp.chromosome = as.character(seqnames(probes.cnv.gr)),
                              snp.position = start(probes.cnv.gr)
  )
  
  
  if("CNVGenotype" %in% genotype.nodes){
    ###################### Replace genotype matrix with CNV genotypes
    (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
    
    CNVBiMaCN <- matrix(2, nrow = length(probes.cnv.gr), 
                        ncol = length(all.samples))
    n <- gdsfmt::add.gdsn(genofile, "CNVgenotype", CNVBiMaCN, replace=TRUE)
    gdsfmt::read.gdsn(n)
    
    if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
      multicoreParam <-  BiocParallel::MulticoreParam(workers = n.cor)
      BiocParallel::bplapply(1:length(all.samples), .writeProbesCNV, BPPARAM = multicoreParam, all.samples=all.samples,
                             genofile=genofile, CNVsGr=CNVsGr, probes.cnv.gr=probes.cnv.gr, n=n)
    }
    
    if(rappdirs:::get_os() == "win"){
      param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
      BiocParallel::bplapply(1:length(all.samples), .writeProbesCNV, BPPARAM = param, all.samples=all.samples, 
                             genofile=genofile, CNVsGr=CNVsGr, probes.cnv.gr=probes.cnv.gr, n=n)
    }
    
    testit(15) ## Wait to make sure that gds is writen in parallel
    SNPRelate::snpgdsClose(genofile)
    
  }
  
  if("CNVgenotypeSNPlike" %in% genotype.nodes){
    
    (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
    
    CNVBiMaCN <- matrix(2, nrow = length(probes.cnv.gr), 
                        ncol = length(all.samples))
    
    n <- gdsfmt::add.gdsn(genofile, "CNVgenotypeSNPlike", CNVBiMaCN, replace=TRUE)
    gdsfmt::read.gdsn(n)
    
    ### Define the chunks of SNPs to import 
    #chunk <- seq(from=1, to=length(probes.cnv.gr), by=chunk.len)
    chunk <- seq(from=1, to=length(probes.cnv.gr), by=length(probes.cnv.gr))
    if(chunk[length(chunk)]!=length(probes.cnv.gr)){
      chunk[length(chunk)+1] <- length(probes.cnv.gr)}
    
    if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
      multicoreParam <-  BiocParallel::MulticoreParam(workers = n.cor)
      BiocParallel::bplapply(1:(length(chunk)-1), .recodeCNVgenotype, BPPARAM = multicoreParam,
                             genofile=genofile, all.samples=all.samples, n=n, chunk=chunk, coding.translate=coding.translate)
    }
    
    if(rappdirs:::get_os() == "win"){
      param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
      BiocParallel::bplapply(1:(length(chunk)-1), .recodeCNVgenotype, BPPARAM = param,
                             genofile=genofile, all.samples=all.samples, n=n, chunk=chunk, coding.translate=coding.translate)
    }
    testit(15) ## Wait to make sure that gds is writen in parallel
    gdsfmt::closefn.gds(genofile)
    
    
  }
  
  if(!all("CNVGenotype" %in% genotype.nodes | "CNVgenotypeSNPlike" %in% genotype.nodes)){
    stop("Undefined method in CNV genotype")
  }
  
  ### Include LRR/BAF info
  # if(lrr == TRUE){
  #importLRR_BAF(all.paths=all.paths, path.files=path.files, list.of.files=list.of.files, n.cor=1) ## TODO: Change n.cor option once parallel BUG fixed
  # }
  
  
  ###################### Include the phenotype 'lo' in the GDS
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  phenotypesSamX <- phenotypesSamX[match(all.samples, phenotypesSamX$samplesPhen),] ### Check order correspondence
  #phenotypesSamX
  gdsfmt::add.gdsn(genofile, name="phenotype", val=phenotypesSamX, replace=TRUE)
  FamID <- FamID[match(all.samples, FamID$samplesPhen),] ### Check order correspondence
  gdsfmt::add.gdsn(genofile, name="FamID", val=as.character(FamID[,2]), replace=TRUE, storage="string")
  FamID <- SexIds[match(all.samples, SexIds$samplesPhen),] ### Check order correspondence
  gdsfmt::add.gdsn(genofile, name="Sex", val=as.character(SexIds[,2]), replace=TRUE, storage="string")
  
  ## Include chr names if needed
  if(!is.null(chr.code.name)){
    gdsfmt::add.gdsn(genofile, name="Chr.names", val=chr.code.name, replace=TRUE)  
  }
  
  if(!is.null(phen.info$pops.names)){
    gdsfmt::add.gdsn(genofile, name="Pops.names", val=phen.info$pops.names, replace=TRUE)
  }
  
  SNPRelate::snpgdsClose(genofile)
  message("GDS writen...", appendLF = FALSE)
  ########################################### END #######################################
  testit(15)
  message("done", appendLF = FALSE)
  return(probes.cnv.gr)} ## BUG - small value in chuck mixing genotypes

#' HELPER - Produce the probes.cnv.gr
#' 
#' 

.prodProbes <- function(phen.info, lo.phe=1, freq.cn=0.01){
  message("initializing ...", appendLF = FALSE)
  
  phenotypesSam <- phen.info$phenotypesSam
  samplesPhen <- phen.info$samplesPhen
  FamID <- phen.info$FamID
  SexIds <- phen.info$SexIds
  all.paths <- phen.info$all.paths
  
  phenotypesSamX <- phenotypesSam[,c(1,(lo.phe+1))]
  #phenotypesSamX <- na.omit(phenotypesSamX) ### Exclude samples without phenotypes. i.e NA
  
  ######################  Import CNVs to data-frame
  
  cnvs <- data.table::fread(file.path(all.paths[1], "CNVOut.txt"), sep="\t", header=FALSE) ### CNV table 
  CNVs <- .CheckConvertCNVs(cnvs)
  
  ####################### Check if the chromosomes are numeric
  
  chr.names <- CNVs$chr
  chr.names <- gsub("chr", "", chr.names)
  
  if(any(is.na(chr.names))){
    stop("Chromosome names should be integers. If they are not, make them integers and include the correspondent names 
         in chr.code.name parameter")
  }
  
  #############################
  
  CNVs$start <- as.integer(as.character(CNVs$start))
  CNVs$end <- as.integer(as.character(CNVs$end))
  CNVsGr <- makeGRangesFromDataFrame(CNVs, keep.extra.columns = TRUE)
  CNVsGr <- CNVsGr[values(CNVsGr)$V5 %in% samplesPhen] ### Subset CNVs in phenotyped samples
  
  ######################  Import SNP map to data-frame
  probes <- data.table::fread(file.path(all.paths[1], "MapPenn.txt"), header=TRUE, sep="\t")
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
  
  return(probes.cnv.gr)}


#' HELPER - Internal function of parallel implementation to produce the .gvar file (absolute copy number) requested for the CNV-GWAS in PLINK
#' @examples
#' Gens <- .prodGvar(lo, genofile, sam.gen, snps, fam.id)

.prodGvar <- function(lo, genofile, sam.gen, snps, fam.id, snp.matrix){
  g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id=sam.gen[[lo]]))
  A <- rep("A", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
  B <- rep("B", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
  CNVg <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))
  CNVg <- CNVg-1
  Dose1 <- rep(1, length(CNVg))
  Dose1[CNVg==-1] <- 0
  CNVg[CNVg==-1] <- 0
  B[CNVg==1] <- "A"
  
  if(snp.matrix=="FALSE"){
    df <- as.data.frame(cbind(as.character(fam.id[[lo]]), sam.gen[[lo]], snps, 
                              A, CNVg, B, Dose1))
    colnames(df) <- c("FID", "IID", "NAME", "ALLELE1", "DOSAGE1", "ALLELE2", "DOSAGE2")}
  
  if(snp.matrix=="TRUE"){ #### TODO - To implement
    ...
  }
  return(df)
}


#' HELPER - Internal function of parallel implementation to produce the .gvar file (LRR) requested for the CNV-GWAS in PLINK
#' @examples
#' Gens <- .prodGvarLRR(lo, genofile, sam.gen, snps, fam.id)

.prodGvarLRR <- function(lo, genofile, sam.gen, snps, fam.id, snp.matrix){
  
  ## Transform LRR calclulations into genotypes
  g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id=sam.gen[[lo]]))
  A <- rep("A", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
  B <- rep("B", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
  LRRg <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "LRR"), start=c(1,lo), count=c(length(g),1))
  LRRg[is.na(LRRg)] <- 0
  
  CNVg <- LRRg
  Dose1 <- rep(0, length(CNVg))
  
  
  if(snp.matrix=="FALSE"){
    df <- as.data.frame(cbind(as.character(fam.id[[lo]]), sam.gen[[lo]], snps, 
                              A, CNVg, B, Dose1))
    colnames(df) <- c("FID", "IID", "NAME", "ALLELE1", "DOSAGE1", "ALLELE2", "DOSAGE2")}
  
  if(snp.matrix=="TRUE"){ #### TODO - To implement
    ...
  }
  return(df)
}


#' HELPER - Produce PLINK probe map in disk
#' 

.prodPLINKmap <- function(all.paths){
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE))
  Pmap <- as.data.frame(cbind(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                              gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
                              0, gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position"))))
  
  colnames(Pmap) <- c("Chr", "NAME", "GD", "Position")
  write.table(Pmap, file.path(all.paths[2], "mydata.map"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  SNPRelate::snpgdsClose(genofile)
}

#' HELPER - Produce the PLINK .gvar file in disk
#' 

.prodPLINKgvar <- function(all.paths, n.cor, snp.matrix, run.lrr=FALSE){
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
  sam.gen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")) ## Extract samples
  snps <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")) ## Extract SNPs
  fam.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "FamID")) ## Extract family
  
  ## Produce gvar file in parallel processing
  # Identify the SO
  if(run.lrr == TRUE){
    if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
      multicoreParam <- BiocParallel::MulticoreParam(workers = n.cor)
      Gens <- suppressMessages(BiocParallel::bplapply(1:length(sam.gen), .prodGvarLRR, BPPARAM = multicoreParam, genofile=genofile, sam.gen=sam.gen,
                                                      snps=snps, fam.id=fam.id, snp.matrix=snp.matrix))
    }
    
    if(rappdirs:::get_os() == "win"){
      param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
      Gens <- suppressMessages(BiocParallel::bplapply(1:length(sam.gen), .prodGvarLRR, BPPARAM = param, genofile=genofile, sam.gen=sam.gen,
                                                      snps=snps, fam.id=fam.id, snp.matrix=snp.matrix))
    }
    
  }else{
    if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
      multicoreParam <- BiocParallel::MulticoreParam(workers = n.cor)
      Gens <- suppressMessages(BiocParallel::bplapply(1:length(sam.gen), .prodGvar, BPPARAM = multicoreParam, genofile=genofile, sam.gen=sam.gen,
                                                      snps=snps, fam.id=fam.id, snp.matrix=snp.matrix))
    }
    
    if(rappdirs:::get_os() == "win"){
      param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
      Gens <- suppressMessages(BiocParallel::bplapply(1:length(sam.gen), .prodGvar, BPPARAM = param, genofile=genofile, sam.gen=sam.gen,
                                                      snps=snps, fam.id=fam.id, snp.matrix=snp.matrix))
    }
    
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
  
  famdf <- as.data.frame(cbind(FAMID, SamGen, SamGen, "NC", Sex, as.numeric(as.character(Phen[,2]))))
  colnames(famdf) <- c("FAID", "IID", "PID", "MID", "Sex", "Phenotype")
  
  write.table(famdf, file.path(all.paths[2], "mydata.fam"), sep = "\t", col.names = TRUE, row.names=FALSE, quote=FALSE)
  SNPRelate::snpgdsClose(genofile)
}

#' HELPER - Run PLINK
#' 

.runPLINK <- function(all.paths){
  
  plinkPath <- paste0("\'", all.paths[2], "/plink\'")
  system(paste(plinkPath, "--gfile", file.path(all.paths[2], "mydata"), paste("--out", file.path(all.paths[2], "plink")), "--noweb"), wait=TRUE)
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=T, readonly=FALSE)) ## read GDS
  sumphen <- paste0("plink.gvar.summary.",names(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")))[[2]])
  
  file.rename(file.path(all.paths[2], "plink.gvar.summary"), file.path(all.paths[2],sumphen))
  SNPRelate::snpgdsClose(genofile)
  
}

#' HELPER - Produce CNV segments
#' 

.prodCNVseg <- function(all.paths, probes.cnv.gr, min.sim){
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
  
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
    
    if(length(snpindex)==1){
      SimiLar.all[[lo]] <- "UNIQUE"
      next
    }
    
    if(length(snpindex)==2){
      g1x <- rbind(g1x, rep(2, ncol(g1x))) # If only two SNPs input a 2n row
    }
    
    g1xM <- g1x[-(nrow(g1x)),]
    
    
    g1CN <- apply(g1xM, 1, function (x) sum(table(x)[names(table(x))!=2]))
    g1x <- as.data.frame(t(sapply(1:nrow(g1x), function(i) replace(g1x[i,], g1x[i,] == 2, paste0(i, 'R')))))
    
    g2x <- g2[snpindex,]
    
    if(length(snpindex)==2){
      g2x <- rbind(g2x, rep(2, ncol(g2x))) # If only two SNPs input a 2n row
    }
    
    g2x <- as.data.frame(t(sapply(1:nrow(g2x), function(i) replace(g2x[i,], g2x[i,] == 2, paste0(i, 'R')))))
    g1x <- g1x[-(nrow(g1x)),]
    g2x <- g2x[-1,]
    ftma <- g1x==g2x
    SimiLar <- apply(ftma, 1, function(x) as.numeric(table(x)["TRUE"]))
    SimiLar[is.na(SimiLar)] <- 0
    
    simX <- as.numeric(SimiLar)/g1CN
    if(length(snpindex)==2){
      simX <- simX[-2]
    }
    
    SimiLar.all[[lo]] <- simX ### CNV genotype similarity between subsequent SNPs
  }
  
  all.segs <- NULL
  for(ch in 1:length(SimiLar.all)){
    SimiLarX <- SimiLar.all[[ch]]
    if(all(SimiLarX=="UNIQUE")){
      all.segs[[ch]] <- "start"
      next
    }
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
          all.segsch[[se+1]] <- "cont"}
        if(SimiLarX[[se]]<min.sim){
          all.segsch[[se+1]] <- "start"}
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
  
  SNPRelate::snpgdsClose(genofile)
  
  return(all.segs.gr)
  
}


#' HELPER - Associate probes with CNV segments and draw p-values
#' 

.assoPrCNV  <- function(all.paths, all.segs.gr, phenotypesSamX, method.m.test, probes.cnv.gr, assign.probe = "min.pvalue"){
  
  segs.pvalue.gr <- unlist(all.segs.gr)
  values(segs.pvalue.gr)$SegName <- seq(from=1, to=length(segs.pvalue.gr), by=1)
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  
  sum.file <- paste0("plink.gvar.summary.",names(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")))[[2]])
  results <- read.table(file.path(all.paths[2], sum.file), header=TRUE) 
  mydata.map <- read.table(file.path(all.paths[2],"mydata.map"), header=TRUE)
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
  freqs.all <- vector("list", length(segs.pvalue.gr))
  
  if(assign.probe == "min.pvalue"){
    for(se in 1:length(segs.pvalue.gr)){
      if(length(subsetByOverlaps(resultspGr, segs.pvalue.gr[se]))>0){
        values.all[[se]] <- min(values(subsetByOverlaps(resultspGr, segs.pvalue.gr[se]))$VALUE)
        Prbx <- subsetByOverlaps(resultspGr, segs.pvalue.gr[se])
        names.all[[se]] <- as.character(Prbx[values(Prbx)$VALUE %in% min(values(subsetByOverlaps(resultspGr, segs.pvalue.gr[se]))$VALUE)]$NAME[[1]])
        probe.sel <- subsetByOverlaps(resultspGr, segs.pvalue.gr[se])
        probe.sel <- probe.sel[order(values(probe.sel)$VALUE)][1] 
        probe.sel <- subsetByOverlaps(probes.cnv.gr, probe.sel)
        freqs.all[[se]] <- as.character(values(probe.sel)$freq[1])}
      else{
        probe.sel <- subsetByOverlaps(probes.cnv.gr, segs.pvalue.gr[se])[1]
        values.all[[se]] <- 1
        names.all[[se]] <- as.character(values(probe.sel)$Name)
        freqs.all[[se]] <- as.character(values(probe.sel)$freq[1])}
    }
  }
  else if(assign.probe == "high.freq"){
    for(se in 1:length(segs.pvalue.gr)){
      probe.sel <- subsetByOverlaps(probes.cnv.gr, segs.pvalue.gr[se])
      probe.sel <- probe.sel[rev(order(values(probe.sel)$freq))][1] 
      resultX <- subsetByOverlaps(resultspGr, probe.sel)[1] ## Take the first probe if the position is duplicated
      values.all[[se]] <- values(resultX)$VALUE
      names.all[[se]] <- as.character(values(probe.sel)$Name)
      freqs.all[[se]] <- as.character(values(probe.sel)$freq[1])
    }
  }
  
  values(segs.pvalue.gr)$MinPvalue <- unlist(values.all)
  values(segs.pvalue.gr)$NameProbe <- unlist(names.all)
  values(segs.pvalue.gr)$Frequency <- unlist(freqs.all)
  
  values(segs.pvalue.gr)$MinPvalueAdjusted <- stats::p.adjust(values(segs.pvalue.gr)$MinPvalue, method = method.m.test, n = length(segs.pvalue.gr))
  segs.pvalue.gr$Phenotype <- names(phenotypesSamX)[[2]]
  SNPRelate::snpgdsClose(genofile)
  
  return(segs.pvalue.gr)
  
}


#' Run the CNV-GWAS 
#'
#' Wraps all the necessary functions to run a CNV-GWAS using the output of \code{\link{setupCnvGWAS}} function
#'  
#' (i) Produces the GDS file containing the genotype information  (if produce.gds == TRUE),
#' (ii) Produces the requested inputs for a PLINK analysis, (iii) run a 
#' CNV-GWAS analysis using linear model implemented in PLINK (http://zzz.bwh.harvard.edu/plink/gvar.shtml) and 
#' (iv) export a QQ-plot displaying the adjusted p-values. In this release only the p-value for the copy number is available  (i.e. 'P(CNP)'). 
#'
#'   
#' @param phen.info Returned by \code{\link{setupCnvGWAS}}
#' @param n.cor Number of cores to be used
#' @param min.sim Minimum CNV genotype distribuition similarity among subsequent probes. Default is 0.95 = 95%
#' @param freq.cn Minimum CNV frequency where 1 = 100%, or all samples deviating from diploid state. Default 0.01 = 1%
#' @param snp.matrix Only FALSE implemented - If TRUE B allele frequencies (BAF) would be used to reconstruct CNV-SNP genotypes
#' @param method.m.test Correction for multiple tests to be used. FDR is default, see \code{\link(p.adjust)} for
#' other methods.
#' @param lo.phe The phenotype to be analyzed in the PhenInfo$phenotypesSam data-frame 
#' @param chr.code.name A data-frame with the integer name in the first column and the original name for each chromosome  
#' @param genotype.nodes Expression data type. Nodes with CNV genotypes to be produced in the gds file. 
#' @param coding.translate For 'CNVgenotypeSNPlike'. If NULL or unrecognized string use only biallelic CNVs. If 'all' code multiallelic CNVs as 0 
#' for loss; 1 for 2n and 2 for gain.
#' @param lrr logical. Import LRR/BAF values from external files
#' @param path.files Folder containing the input CNV files used for the CNV calling (i.e. one text file with 5 collumns 
#' for each sample). Columns should contain (i) probe name, (ii) Chromosome, (iii) Position, (iv) LRR and (v) BAF. Only if \sQuote{lrr} == TRUE
#' @param list.of.files Data-frame with two columns where the (i) is the file name with signals and (ii) is the 
#' correspondent name of the sample in the gds file. Only if \sQuote{lrr} == TRUE
#' @param produce.gds logical. If TRUE produce a new gds, if FALSE use gds previously created   
#' @param run.lrr If TRUE use LRR values instead absolute copy numbers in the association
#' @param assign.probe \sQuote{min.pvalue} or \sQuote{high.freq} to represent the CNV segment
#' @param select.population Select a subset of populations to be used in the GWAS. Use the population name or names
#' @return The CNV segments and the representative probes and their respective p-value
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples
#' 
#' # Load phenotype-CNV information
#' 
#' folder.package <- system.file('extdata', package = 'cnvAnalyzeR')
#' 
#' phen.info <- setupCnvGWAS(name="Example", phen.loc=file.path(folder.package, "Pheno.txt"), 
#' cnv.out.loc = file.path(folder.package, "CNVOut.txt"), 
#' map.loc = file.path(folder.package, "MapPenn.txt"))
#' 
#' # Define chr correspondence to numeric, if necessary
#' df <- "16 1A
#' 25 4A
#' 29 25LG1
#' 30 25LG2
#' 31 LGE22"
#' 
#' chr.code.name <- read.table(text=df, header=F)
#' 
#' segs.pvalue.gr <- cnvGWAS(phen.info, chr.code.name=chr.code.name)
#'  
#' @export

cnvGWAS <- function(phen.info, n.cor=1, min.sim = 0.95, freq.cn = 0.01, snp.matrix=FALSE, 
                    method.m.test = "fdr", lo.phe=1, chr.code.name=NULL, genotype.nodes="CNVGenotype",
                    #coding.translate="all", lrr = FALSE, path.files=NULL, list.of.files=NULL, 
                    coding.translate="all", path.files=NULL, list.of.files=NULL, 
                    produce.gds=TRUE, run.lrr=FALSE, assign.probe="min.pvalue",
                    select.population = NULL){
  
  if(!is.null(select.population)){
    phen.info$samplesPhen <- phen.info$samplesPhen[which(phen.info$pops.names %in% select.population)]
    phen.info$phenotypesdf <- phen.info$phenotypesdf[which(phen.info$pops.names %in% select.population), ]
    phen.info$phenotypesSam <- phen.info$phenotypesSam[which(phen.info$pops.names %in% select.population), ]
    phen.info$FamID <- phen.info$FamID[which(phen.info$pops.names %in% select.population), ]
    phen.info$SexIds <- phen.info$SexIds[which(phen.info$pops.names %in% select.population), ]
    phen.info$pops.names <- phen.info$pops.names[which(phen.info$pops.names %in% select.population), ]
  }
  
  phenotypesSam <- phen.info$phenotypesSam
  phenotypesSamX <- phenotypesSam[,c(1,(lo.phe+1))]
  phenotypesSamX <- na.omit(phenotypesSamX) ### Exclude samples without phenotypes. i.e NA
  all.paths <- phen.info$all.paths
  
  #################################### Produce the GDS for a given phenotype
  
  if(produce.gds == TRUE){
    message("Produce the GDS for a given phenotype")
    
    probes.cnv.gr <- prodGdsCnv(phen.info=phen.info, freq.cn=freq.cn, snp.matrix=snp.matrix, n.cor=n.cor, 
                                lo.phe=lo.phe, chr.code.name=chr.code.name, genotype.nodes=genotype.nodes,
                                #coding.translate=coding.translate, lrr=lrr)
                                coding.translate=coding.translate)
  }else if(produce.gds == FALSE){
    message("Using existent gds file") 
    probes.cnv.gr <- .prodProbes(phen.info, lo.phe, freq.cn)
  }else{
    stop("TRUE/FALSE needed")  
  }
  
  ##################################### Produce PLINK map ##################
  message("Produce PLINK map")
  
  .prodPLINKmap(all.paths)
  
  ########################################### END #######################################
  
  ##################################### Produce gvar to use as PLINK input #########
  message("Produce gvar to use as PLINK input")
  
  .prodPLINKgvar(all.paths, n.cor, snp.matrix, run.lrr)
  
  ########################################### END #######################################
  
  ##################################### Produce fam (phenotype) to use as PLINK input #############
  message("Produce fam (phenotype) to use as PLINK input")
  
  .prodPLINKfam(all.paths)
  
  ########################################### END #######################################
  
  ##################################### Run PLINK #############
  message("Run PLINK")
  .runPLINK(all.paths)
  
  ########################################### END #######################################
  
  ##################################### Produce CNV segments #############
  message("Produce CNV segments")
  all.segs.gr <- .prodCNVseg(all.paths, probes.cnv.gr, min.sim)
  
  ########################################### END #######################################
  
  ##################################### Associate SNPs with CNV segments #############
  message("Associate SNPs with CNV segments")
  segs.pvalue.gr <- .assoPrCNV(all.paths, all.segs.gr, phenotypesSamX, method.m.test, probes.cnv.gr, assign.probe)
  
  ######################################################## END ###########################################
  
  #################################### Plot the QQ-plot of the analysis
  message("Plot the QQ-plot of the analysis")
  
  pdf(file.path(all.paths[3], paste0(unique(values(segs.pvalue.gr)$Phenotype), "-LRR-", run.lrr, "QQ-PLOT.pdf")))
  qq.plot.pdf <- qqunif.plot(values(segs.pvalue.gr)$MinPvalue, auto.key=list(corner=c(.95,.05)))
  print(qq.plot.pdf)
  dev.off()
  
  ######################################################## END ###########################################
  
  #################################### Reconvert the chrs to original names if applicable
  
  if(!is.null(chr.code.name)){
    (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
    chr.code.name < gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Chr.names"))
    
    segs.pvalue <- data.frame(segs.pvalue.gr)
    segs.pvalue <- segs.pvalue[ , -which(names(segs.pvalue) %in% c("strand"))]
    
    for(lopN in seq_len(nrow(chr.code.name))){
      segs.pvalue$seqnames <- gsub(chr.code.name[lopN,1],chr.code.name[lopN,2], segs.pvalue$seqnames)
    }
    
    segs.pvalue.gr <- makeGRangesFromDataFrame(segs.pvalue, keep.extra.columns = TRUE)
    
    SNPRelate::snpgdsClose(genofile)
  }
  
  ######################################################## END ########################################
  
  ######################################################## END ###########################################
  segs.pvalue.gr <- segs.pvalue.gr[order(values(segs.pvalue.gr)$MinPvalue)]  
  
  return(segs.pvalue.gr)
}


#' HELPER Extract the CNV genotypes in dosage format from GDS
#' @param genofile is the loaded gds file
#' @param snp.id is the SNP id to retrieve 
#' @param node.to.extract  \sQuote{CNVgenotype} or  \sQuote{CNVgenotypeSNPlike}. Default is  \sQuote{CNVgenotype}

.snpgdsGetGenoCNV <- function(genofile, snp.id, node.to.extract="CNVgenotype"){
  
  map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")), 
                             chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                             position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")), 
                             probes = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
  ))
  
  all.samples <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
  
  if(node.to.extract == "CNVgenotype"){
    gens <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(snp.id,1), count=c(1,length(all.samples))))
  }else if(node.to.extract == "CNVgenotypeSNPlike"){
    gens <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotypeSNPlike"), start=c(snp.id,1), count=c(1,length(all.samples))))
  }else{
    stop("Undefined node to extract CNV genotype")
  }
  
  return(gens)
  
}


#' Import LRR and BAF from text files used in the CNV analysis
#' 
#' This function imports the LRR/BAF values and create a node for each one in the GDS file at the working folder 'Inputs'
#' created by the \code{\link{setupCnvGWAS}} function. Once imported, the LRR values can be used to perform a GWAS directly as
#' an alternative to copy number dosage
#'
#' @param all.paths Object returned from \code{CreateFolderTree} function with the working folder tree 
#' @param path.files Folder containing the input CNV files used for the CNV calling (i.e. one text file with 5 collumns 
#' for each sample). Columns should contain (i) probe name, (ii) Chromosome, (iii) Position, (iv) LRR and (v) BAF
#' @param list.of.files Data-frame with two columns where the (i) is the file name with signals and (ii) is the 
#' correspondent name of the sample in the gds file
#' @param n.cor Number of cores to be used. Only one core implemented in this version. 
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples
#' 
#' # Load phenotype-CNV information
#' 
#' folder.package <- system.file('extdata', package = 'cnvAnalyzeR')
#' 
#' phen.info <- setupCnvGWAS(name="Example", phen.loc=file.path(folder.package, "Pheno.txt"), 
#' cnv.out.loc = file.path(folder.package, "CNVOut.txt"), 
#' map.loc = file.path(folder.package, "MapPenn.txt"))
#' 
#' # Extract path names
#' all.paths <- phen.info$all.paths
#' 
#' # List files to import LRR/BAF 
#' list.of.files <- list.files(path = folder.package, pattern = "\\.cnv.txt.adjusted$")
#' list.of.files <- as.data.frame(list.of.files)
#' colnames(list.of.files)[1] <- "file.names"
#' list.of.files$sample.names <- gsub(".cnv.txt.adjusted","", list.of.files$file.names)
#' 
#' # All missing samples will have LRR = '0' and BAF = '0.5' in all SNPs listed in the GDS file
#' importLRR_BAF(all.paths, folder.package, list.of.files)
#' 
#' # Read the GDS to check if the LRR/BAF nodes were added
#' (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
#' 
#' @export

importLRR_BAF <- function(all.paths, path.files, list.of.files, n.cor=1){
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
  
}


#' HELPER Create the pairwise genotypes for epistasis analysis - TODO
#' 
#' 
#' 
