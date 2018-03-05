##################
# Author: Vinicius Henrique da Silva
# Script description: Produce PLINK inputs for CNV segment association - Modules CNV analyzeR
# Date: October 9, 2017
# Code: CNVASSOPACK002
###################

library(SNPRelate)
library(data.table)
library(GenomicRanges)
library(doMC)
library(GWASTools)
library(doBy)

##### Parameters ####
Ncor <- 10 ## Number of cores to be used
Inputset <- "/home/vhsilva/Projects/CNV_chicken/Inputs" ### Folder with inputs
Plinkset <- "/home/vhsilva/Projects/CNV_chicken/PLINK/plink-1.07-x86_64" ### PLINK folder
ResSet <- "/home/vhsilva/Projects/CNV_chicken/Results" ### Results folder
SNPmatrix <- "FALSE" ### Only FALSE implemented - If TRUE the SNP genotypes would be used
FreqCN <- 0.01 ### Minimum frequency - 1 = 100%
MinSim <- 0.95 ### Minimum CNV genotype distribuition similarity among subsequent SNPs - 1 = 100%
#########


##################################### (i) Load the CNVs and phenotypes  ######
setwd(Inputset)
all <- fread("Allinformation_GWAS_popTT.txt", header=TRUE, sep="\t") ### All info - Input4
all$sample_id <- gsub(".cnv.txt-edit.txt.adjusted", "",all$sample_id)
samplesPhen <- unique(all$sample_id) ## Greb sample names
#length(samplesPhen)
phenotypes <- names(all[, 26:41]) ## Greb phenotype names
Samphenotypes <- paste0("samples_",names(all[, 26:41])) ## Greb sample names with a given phenotype
phenotypesdf <- as.data.frame(all[, 26:41]) ### Greb phenotypes values
phenotypesdf <- phenotypesdf[!duplicated(phenotypesdf),]
phenotypesSam <- as.data.frame(cbind(samplesPhen, phenotypesdf))
head(phenotypesSam)
FamID <- as.data.frame(cbind(samplesPhen, all$NLABfather[as.numeric(rownames(phenotypesSam))]))
SexIds <- as.data.frame(cbind(samplesPhen,all$sex[as.numeric(rownames(phenotypesSam))]))

########################################### END #######################################
SegsGr.all <- GRangesList()

for(loPhe in 1:ncol(phenotypesdf)){ ### Subsequent items need to run independly for each phenotype - i.e. different number of phenotyped samples
  #loPhe =1
  
  phenotypesSamX <- phenotypesSam[,c(1,(loPhe+1))]
  phenotypesSamX <- na.omit(phenotypesSamX) ### Exclude samples without phenotypes. i.e NA
  
  ##################################### (ii) Produce CNV GDS for phenotyped samples #####
  ######################  Import CNVs to data-frame
  setwd(Inputset)
  cnvs <- read.table("inputGWAS_CNV_popTT.txt", sep="", header=F) ### PennCNV output - Input1 - OK
  cnvs$V5 <- gsub(".cnv.txt-edit.txt.adjusted", "",cnvs$V5)
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
  CNVsGr <- makeGRangesFromDataFrame(CNVs, keep.extra.columns = T)
  CNVsGr <- CNVsGr[values(CNVsGr)$V5 %in% samplesPhen] ### Subset CNVs in phenotyped samples
  
  ######################  Import SNP map to data-frame
  setwd(Inputset)
  probes <- fread("mapa_chicken.txt", header=T, sep="\t") #### Map of probes - Input2 - OK
  probes <- as.data.frame(probes)
  probes$Position <- as.numeric(as.character(probes$Position))
  probes <- probes[complete.cases(probes), ]
  probesGr <- makeGRangesFromDataFrame(probes, seqnames.field = "Chr", start.field = "Position", end.field = "Position",
                                       keep.extra.columns = T)
  
  Samples <- unique(as.character(values(CNVsGr)$V5))
  
  ###################### Select probes  within CNVs
  probesCNV <- suppressWarnings(as.character(values(subsetByOverlaps(probesGr, 
                                                                     CNVsGr))$Name))
  probesCNV <- unique(unlist(probesCNV))
  length(probesCNV)
  probesCNVGr <- probesGr[values(probesGr)$Name %in% probesCNV]
  
  counts <- as.integer(countOverlaps(probesCNVGr, CNVsGr))
  #subsetByOverlaps(probesCNVGr[[1]], CNVsGr)
  #subsetByOverlaps(CNVsGr, probesCNVGr[1])
  values(probesCNVGr)$freq <- counts
  
  ##### Subset by frequency 
  NumSam <- FreqCN*length(Samples)
  probesCNVGr <- probesCNVGr[values(probesCNVGr)$freq>=NumSam]
  
  ##### Order the probes
  probesCNVGr <- sortSeqlevels(probesCNVGr)
  probesCNVGr <- sort(probesCNVGr)
  probesCNVGr
  
  values(probesCNVGr)$snp.id <- seq(1:length(probesCNVGr))
  
  ###################### SNP genotype matrix not available 
  if(SNPmatrix == "FALSE"){
    CNVBiMa <- matrix(2, nrow = length(probesCNVGr), 
                      ncol = length(Samples))}
  
  ###################### SNP genotype matrix available - TODO
  if(SNPmatrix == "TRUE"){
    ...
  }
  
  ###################### Create a GDS with chr and SNP names to numeric
  snpgdsCreateGeno("CNV.gds", 
                   genmat = CNVBiMa,
                   sample.id = Samples, 
                   #sample.id = SamplesS,
                   snp.id = as.character(values(probesCNVGr)$snp.id),
                   snp.rs.id = as.character(values(probesCNVGr)$Name),
                   snp.chromosome = as.character(seqnames(probesCNVGr)),
                   snp.position = start(probesCNVGr)
  )
  
  testit <- function(x)
  {
    p1 <- proc.time()
    Sys.sleep(x)
    proc.time() - p1 
  }
  testit(15)
  
  ###################### Replace genotype matrix with CNV genotypes
  setwd(Inputset)
  (genofile <- snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE))
  
  CNVBiMaCN <- matrix(2, nrow = length(probesCNVGr), 
                      ncol = length(Samples))
  n <- add.gdsn(genofile, "CNVgenotype", CNVBiMaCN, replace=TRUE)
  read.gdsn(n)
  
  gback <- as.numeric(CNVBiMaCN[,1])
  registerDoMC(Ncor)
  probesCNV  <- foreach(lo=1:length(Samples), .verbose=FALSE, .errorhandling=c('pass')) %dopar% {
    print(paste(lo,"of",length(Samples)))
    #sampleX <- SamplesS[[lo]]
    sampleX <- Samples[[lo]]
    #g <- gback
    g <- as.numeric(snpgdsGetGeno(genofile, sample.id=sampleX))
    CNVSamByState <- split(CNVsGr[values(CNVsGr)$V5 == Samples[lo]], 
                           as.character(elementMetadata(CNVsGr[values(CNVsGr)$V5 == Samples[lo]])$V4))
    for(slo in 1:length(CNVSamByState)){
      state <- unique(substr(x=as.character(values(CNVSamByState[[slo]])$V4),
                             start=nchar(as.character(values(CNVSamByState[[slo]])$V4)),
                             stop=nchar(as.character(values(CNVSamByState[[slo]])$V4))))
      SamCNVSNP <- values(subsetByOverlaps(probesCNVGr, CNVSamByState[[slo]]))$snp.id
      g[SamCNVSNP] <- state
    }
    g <- as.integer(g)
    write.gdsn(n, g, start=c(1,lo), count=c(length(g),1))
  }
  
  
  testit(15)
  
  ###################### Include the phenotype 'lo' in the GDS
  phenotypesSamX <- phenotypesSamX[match(Samples, phenotypesSamX$samplesPhen),] ### Check order correspondence
  #phenotypesSamX$samplesPhen <- Samples
  #phenotypesSamX
  add.gdsn(genofile, name="phenotype", val=phenotypesSamX, replace=TRUE)
  FamID <- FamID[match(Samples, FamID$samplesPhen),] ### Check order correspondence
  add.gdsn(genofile, name="FamID", val=FamID[,2], replace=TRUE)
  FamID <- SexIds[match(Samples, SexIds$samplesPhen),] ### Check order correspondence
  add.gdsn(genofile, name="Sex", val=SexIds[,2], replace=TRUE)
  snpgdsClose(genofile)
  ########################################### END #######################################
  
  ##################################### (iii) Produce PLINK map ##################
  setwd(Inputset)
  (genofile <- snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE))
  Pmap <- as.data.frame(cbind(read.gdsn(index.gdsn(genofile, "snp.chromosome")), 
                              read.gdsn(index.gdsn(genofile, "snp.rs.id")),
                              0, read.gdsn(index.gdsn(genofile, "snp.position"))))
  
  colnames(Pmap) <- c("Chr", "NAME", "GD", "Position")
  setwd(Plinkset)
  write.table(Pmap, "mydata.map", sep="\t", col.names=T, row.names=F, quote=F)
  snpgdsClose(genofile)
  ########################################### END #######################################
  
  ##################################### (iv) Produce gvar to use as PLINK input #############
  setwd(Inputset)
  (genofile <- snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE)) ## read GDS
  SamGen <- read.gdsn(index.gdsn(genofile, "sample.id")) ## Extract samples
  SNPs <- read.gdsn(index.gdsn(genofile, "snp.rs.id")) ## Extract SNPs
  FAMID <- read.gdsn(index.gdsn(genofile, "FamID")) ## Extract family
  #A <- "A"
  #B <- "B"
  #Dose1 <- 1
  
  #for(lo in 1:SamGen){
  Gens  <- foreach(lo=1:length(SamGen), .verbose=FALSE, .errorhandling=c('pass')) %dopar% {
    print(paste0("Sample ",lo," of ", length(SamGen)))
    g <- as.numeric(snpgdsGetGeno(genofile, sample.id=SamGen[[lo]]))
    A <- rep("A", length(read.gdsn(index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
    B <- rep("B", length(read.gdsn(index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
    CNVg <- read.gdsn(index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))
    CNVg <- CNVg-1
    Dose1 <- rep(1, length(CNVg))
    Dose1[CNVg==-1] <- 0
    #A[CNVg==-1] <- "0"
    #B[CNVg==-1] <- "0"
    CNVg[CNVg==-1] <- 0
    B[CNVg==1] <- "A"
    
    if(SNPmatrix=="FALSE"){
      df <- as.data.frame(cbind(as.character(FAMID[[lo]]), SamGen[[lo]], SNPs, 
                                A, CNVg, B, Dose1))
      colnames(df) <- c("FID", "IID", "NAME", "ALLELE1", "DOSAGE1", "ALLELE2", "DOSAGE2")}
    
    if(SNPmatrix=="TRUE"){ #### TODO - To implement
      ...
    }
    df
  }
  
  gentype.all <- rbindlist(Gens)
  nrow(gentype.all)
  gentype.all1 <- as.data.frame(gentype.all)
  #table(gentype.all$DOSAGE1)
  
  setwd(Plinkset)
  write.table(gentype.all1, "mydata.gvar", sep="\t", col.names=T, row.names=F, quote=F)
  #fwrite(gentype.all, "mydata.gvar", sep = "\t", col.names = T, row.names=F, quote=F, nThread=Ncor)
  snpgdsClose(genofile)
  ########################################### END #######################################
  
  ##################################### (v) Produce fam (phenotype) to use as PLINK input #############
  setwd(Inputset)
  (genofile <- snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE)) ## read GDS
  SamGen <- read.gdsn(index.gdsn(genofile, "sample.id")) ## Extract samples
  SNPs <- read.gdsn(index.gdsn(genofile, "snp.rs.id")) ## Extract SNPs
  FAMID <- as.character(read.gdsn(index.gdsn(genofile, "FamID"))) ## Extract family
  #FAMID <- rep(NA, length(as.character(read.gdsn(index.gdsn(genofile, "FamID"))))) ## Extract family
  Sex <- read.gdsn(index.gdsn(genofile, "Sex")) ## Extract sex
  Phen <- read.gdsn(index.gdsn(genofile, "phenotype")) ## Extract phenotype
  
  setwd(Plinkset)
  famdf <- as.data.frame(cbind(FAMID, SamGen, SamGen, "NC", Sex, Phen[,2]))
  #famdf <- as.data.frame(cbind(FAMID, SamGen, FAMID, "NC", Sex, Phen[,2]))
  colnames(famdf) <- c("FAID", "IID", "PID", "MID", "Sex", "Phenotype")
  
  
  setwd(Plinkset)
  write.table(famdf, "mydata.fam", sep = "\t", col.names = T, row.names=F, quote=F)
  #fwrite(famdf, "mydata.fam", sep = "\t", col.names = T, row.names=F, quote=F, nThread=Ncor)
  snpgdsClose(genofile)
  ########################################### END #######################################
  
  ##################################### (vi) Run PLINK #############
  setwd(Plinkset)
  system("./plink --gfile mydata --noweb", wait=TRUE)
  
  setwd(Inputset)
  (genofile <- snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE)) ## read GDS
  sumphen <- paste0("plink.gvar.summary.",names(read.gdsn(index.gdsn(genofile, "phenotype")))[[2]])
  cod <- c("mv plink.gvar.summary NEW")
  cod <- gsub("NEW", sumphen, cod)
  setwd(Plinkset)
  system(cod)
  snpgdsClose(genofile)
  ########################################### END #######################################
  
  ##################################### (vii) Produce CNV segments #############
  setwd(Inputset)
  (genofile <- snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE)) ## read GDS
  
  chrx <- as.data.frame(cbind(read.gdsn(index.gdsn(genofile, "snp.chromosome")), 
                              seq(from=1, to=length(read.gdsn(index.gdsn(genofile, "snp.chromosome"))), by=1) ))
  g1 <- read.gdsn(index.gdsn(genofile, "CNVgenotype"))
  g2 <- read.gdsn(index.gdsn(genofile, "CNVgenotype"))
  
  chrx$V1 <- factor(chrx$V1, levels=unique(chrx$V1))
  chrx.split <- split(chrx, chrx$V1)
  #sapply(chrx.split[[lo]], class)
  
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
    print(paste0("Chr",ch))
    SimiLarX <- SimiLar.all[[ch]]
    nextP <- NULL
    all.segsch <- NULL
    for(se in 1:length(SimiLarX)){
      print(paste("SNP",se,"of",length(SimiLarX)))
      if(SimiLarX[[se]]>=MinSim){
        nextP[[se]] <- "TRUE"
      }
      if(SimiLarX[[se]]<MinSim){
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
        if(SimiLarX[[se]]>=MinSim){
          all.segsch[[se+1]] <- "start"}
        if(SimiLarX[[se]]<MinSim){
          all.segsch[[se+1]] <- "cont"}
      }
      
    }
    all.segs[[ch]] <- all.segsch
  }
  
  length(unlist(all.segs))
  
  values(probesCNVGr)$starts.seg <- unlist(all.segs)
  
  pstarts <- probesCNVGr[values(probesCNVGr)$starts.seg == "start"]
  pstarts.split = split(pstarts, seqnames(pstarts))
  
  all.seqsGr <- GRangesList()
  #all.seqPrbs <- GRangesList()
  for(chx in 1:length(pstarts.split)){
    print(paste("chr", chx, "segs"))
    pstartsx <- pstarts.split[[chx]]
    if(length(pstartsx)==0){
      all.seqsGr[[chx]] <- pstartsx
      next
    }
    all.seqsGrX <- GRangesList()
    #all.seqPrbsX <- GRangesList()
    for(st in 1:length(pstartsx)){
      print(st)
      prx <- pstartsx[st]
      if(st < length(pstartsx)){
        LimPrx <- start(pstartsx[st+1])-1}
      if(st == length(pstartsx)){
        LimPrx <-  max(start(probesCNVGr[seqnames(probesCNVGr) == seqnames(prx)]))
      }
      seqLar <- GRanges(as.character(seqnames(prx)), IRanges(start(prx), LimPrx))
      seqPrbs <- subsetByOverlaps(probesCNVGr,seqLar)
      seqx <- GRanges(as.character(seqnames(prx)), IRanges(start(prx), start(seqPrbs[length(seqPrbs)])))
      all.seqsGrX[[st]] <- seqx
      #all.seqPrbsX[[st]] <- seqPrbs
    }
    all.seqsGr[[chx]] <- unlist(all.seqsGrX)
    #all.seqPrbs[[chrx]] <- unlist(all.seqPrbsX)
  }
  
  setwd("/home/vhsilva/Projects/CNV_chicken/Results")
  save(all.seqsGr, file="segs.rda")
  snpgdsClose(genofile)
  
  ########################################### END #######################################
  
  ##################################### (viii) Associate SNPs with CNV segments #############
  
  setwd("/home/vhsilva/Projects/CNV_chicken/Results")
  load(file="segs.rda")
  SegsGr <- unlist(all.seqsGr)
  #PrbGr <- unlist(all.seqPrbs)
  values(SegsGr)$SegName <- seq(from=1, to=length(SegsGr), by=1)
  
  setwd(Inputset)
  (genofile <- snpgdsOpen("CNV.gds", allow.fork=T, readonly=FALSE))
  
  setwd("/home/vhsilva/Projects/CNV_chicken/PLINK/plink-1.07-x86_64")
  results <- read.table(paste0("plink.gvar.summary.",names(read.gdsn(index.gdsn(genofile, "phenotype")))[[2]]), header=T) ### FIX WITH GDS NAME
  mydata.map <- read.table("mydata.map", header=T)
  resultsp <- merge(results, mydata.map, by.x="NAME", by.y="NAME", all.x=T, all.y=F)
  tail(resultsp)
  #resultsp <- resultsp[ grep("P\\(", resultsp$FIELD),] ### some probes disapear n 
  #resultsp <- resultsp[ grep("CNV", resultsp$FIELD),] 
  
  #### Choose the method
  resultsp <- resultsp[ grep("P\\(CNP\\)", resultsp$FIELD),]
  resultsp$VALUE <- as.numeric(as.character(resultsp$VALUE))
  resultsp <- resultsp[order(resultsp$VALUE, decreasing=F),]
  resultspGr <- makeGRangesFromDataFrame(resultsp, seqnames.field = "Chr",
                                         start.field = "Position", end.field = "Position",
                                         keep.extra.columns = T)
  
  ######## Assign the lowest p-value to the CNV segment
  
  values.all <- vector("list", length(SegsGr))
  names.all <- vector("list", length(SegsGr))
  for(se in 1:length(SegsGr)){
    print(se)
    values.all[[se]] <- min(values(subsetByOverlaps(resultspGr, SegsGr[se]))$VALUE)
    Prbx <- subsetByOverlaps(resultspGr, SegsGr[se])
    names.all[[se]] <- as.character(Prbx[values(Prbx)$VALUE %in% min(values(subsetByOverlaps(resultspGr, SegsGr[se]))$VALUE)]$NAME[[1]])
  }
  
  values(SegsGr)$MinPvalue <- unlist(values.all)
  values(SegsGr)$NameProbe <- unlist(names.all)
  values(SegsGr)$MinPvalueAdjusted <- p.adjust(values(SegsGr)$MinPvalue, method ="fdr", n = length(SegsGr))
  SegsGr$Phenotype <- names(phenotypesSamX)[[2]]
  
  SegsGr.all[[loPhe]] <- SegsGr
  snpgdsClose(genofile)
  
  ######################################################## END ###########################################
}

setwd(ResSet)
save(SegsGr.all, file="CNVSegAssociation.rda")
