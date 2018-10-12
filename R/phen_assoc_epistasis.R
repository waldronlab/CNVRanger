
#' Run the CNV-GWAS 
#'
#' Wraps all the necessary functions to run a CNV-GWAS using the output of the \code{\link{setupCnvGWAS}} function
#'  
#' (i) Produces the GDS file containing the genotype information  (if produce.gds == TRUE),
#' (ii) Produces the requested inputs for a linear (PLINK, http://zzz.bwh.harvard.edu/plink/gvar.shtml) or linear mixed 
#' model analysis, (iii) run a CNV-GWAS analysis and (iv) export a QQ-plot displaying the adjusted p-values. In this release, 
#' SNP genotypes are not considered and only the p-value for the copy number is available  (i.e. 'P(CNP)' in PLINK). 
#'   
#' @param phen.info Returned by \code{\link{setupCnvGWAS}}
#' @param n.cor Number of cores to be used
#' @param min.sim Minimum CNV genotype distribuition similarity among subsequent probes. Default is 0.95 (i.e. 95\%)
#' @param freq.cn Minimum CNV frequency where 1 (i.e. 100\%), or all samples deviating from diploid state. Default 0.01 (i.e. 1\%)
#' @param snp.matrix Only FALSE implemented - If TRUE B allele frequencies (BAF) would be used to reconstruct CNV-SNP genotypes
#' @param method.m.test Correction for multiple tests to be used. FDR is default, see \code{\link{(p.adjust)}} for
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
#' @param epistasis.mode logical. In addition, perform the GWAS analysis using pairwise combination of the genotypes of all CNV segments to
#' understand epistatic effect among CNVs. Default is FALSE.
#' @param both.up.down logical. Probe genotype similarity considering the up- and downstream vicinity probes. If FALSE will only
#' downstream probes are considered
#' @param cnvr.specific Compare probes inside the same CNVR. Recommended to run with \sQuote{min.sim}>1
#' @param exclude.outliers logical. Use quadrants to exclude samples with extreme phenotypes (i.e. outliers)
#' @param minq Minimum quadrant to exclude. Only if \sQuote{exclude.outliers} is TRUE
#' @param minq Maximum quadrant to exclude. Only if \sQuote{exclude.outliers} is TRUE
#' @param log.pheno logical. log2 transform the phenotypes
#' @param verbose Show progress in the analysis
#' @return The CNV segments and the representative probes with their respective assigned p-value
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
#' segs.pvalue.gr <- cnvGWAS.V2(phen.info, chr.code.name=chr.code.name)
#'  
#' @export

cnvGWAS.V2 <- function(phen.info, n.cor=1, min.sim = 0.95, freq.cn = 0.01, snp.matrix=FALSE, 
                    method.m.test = "fdr", lo.phe=1, chr.code.name=NULL, genotype.nodes="CNVGenotype",
                    coding.translate="all", path.files=NULL, list.of.files=NULL, 
                    produce.gds=TRUE, run.lrr=FALSE, assign.probe="min.pvalue",
                    select.population = NULL, epistasis.mode =FALSE, correct.inflation=FALSE,
                    both.up.down=FALSE, cnvr.specific=FALSE, exclude.outliers=FALSE,
                    minq=.01, maxq=.99, log.pheno=FALSE, verbose=FALSE){
  
  if(!is.null(select.population)){
    phen.info$samplesPhen <- phen.info$samplesPhen[which(phen.info$pops.names %in% select.population)]
    phen.info$phenotypesdf <- phen.info$phenotypesdf[which(phen.info$pops.names %in% select.population), ]
    phen.info$phenotypesSam <- phen.info$phenotypesSam[which(phen.info$pops.names %in% select.population), ]
    phen.info$FamID <- phen.info$FamID[which(phen.info$pops.names %in% select.population), ]
    phen.info$SexIds <- phen.info$SexIds[which(phen.info$pops.names %in% select.population), ]
    phen.info$pops.names <- phen.info$pops.names[which(phen.info$pops.names %in% select.population), ]
  }
  
  ### Exclude outliers
  if(exclude.outliers){
  fun <- function(x, minq, maxq){
    quantiles <- quantile( x, c(minq, maxq) )
    x[ x < quantiles[1] ] <- quantiles[1]
    x[ x > quantiles[2] ] <- quantiles[2]
    x
  }
  
  phenos <- suppressWarnings(as.numeric(as.character(phen.info$phenotypesdf[,1])))
  phenos[which(phenos==-9)] <- NA
  if(log.pheno){
    phenos <- log(2*min(phenos, na.rm=TRUE)*-1+phenos) ## TODO - Re-check
  }
  phenos.na <- phenos
  phenos <- as.numeric(na.omit(phenos))
  
  phenos.na[which(phenos.na<min(fun(phenos, minq, maxq)))] <- NA
  phenos.na[which(phenos.na>max(fun(phenos,minq, maxq)))] <- NA
  
  #phenos.na[is.na(phenos.na)] <- -9 ## Trying to avoid -9 in the phenotype dataset
  phen.info$phenotypesSam[,(lo.phe+1)] <- phenos.na 
  phen.info$phenotypesdf[,1] <- phenos.na
  }
  ### END 
  
  phenotypesSam <- phen.info$phenotypesSam
  phenotypesSamX <- phenotypesSam[,c(1,(lo.phe+1))]
  phenotypesSamX <- na.omit(phenotypesSamX) ### Exclude samples without phenotypes. i.e NA
  all.paths <- phen.info$all.paths
  
  #################################### Produce the GDS for a given phenotype
  
  if(produce.gds == TRUE){
    if(verbose==TRUE){message("Produce the GDS for a given phenotype")}
    
    probes.cnv.gr <- prodGdsCnv(phen.info=phen.info, freq.cn=freq.cn, snp.matrix=snp.matrix, 
                                lo.phe=lo.phe, chr.code.name=chr.code.name, genotype.nodes=genotype.nodes,
                                coding.translate=coding.translate)
  }else if(produce.gds == FALSE){
    if(verbose==TRUE){
      message("Using existent gds file") }
    probes.cnv.gr <- .prodProbes(phen.info, lo.phe, freq.cn)
  }else{
    stop("TRUE/FALSE needed")  
  }
  
  if(PLINK){
  ##################################### Produce PLINK map ##################
  if(verbose==TRUE){
    message("Produce PLINK map")}
  .prodPLINKmap(all.paths)
  
  ########################################### END #######################################
  
  ##################################### Produce gvar to use as PLINK input #########
  if(verbose==TRUE){message("Produce gvar to use as PLINK input")}
  .prodPLINKgvar(all.paths, n.cor, snp.matrix, run.lrr)
  
  ########################################### END #######################################
  
  ##################################### Produce fam (phenotype) to use as PLINK input #############
  if(verbose==TRUE){message("Produce fam (phenotype) to use as PLINK input")}
  .prodPLINKfam(all.paths)
  
  ########################################### END #######################################
  
  ##################################### Run PLINK #############
  if(verbose==TRUE){message("Run PLINK")}
  .runPLINK(all.paths)
  
  ########################################### END #######################################
  
  }else{
  
  ### Produce inputs to the lme4qtl analysis
    
  ### Run the LMM analysis for all probes  
      
  }
  
  ##################################### Produce CNV segments #############
  if(verbose==TRUE){message("Produce CNV segments")}
  all.segs.gr <- .prodCNVseg(all.paths, probes.cnv.gr, min.sim, both.up.down)
  
  ########################################### END #######################################
  
  ##################################### Associate SNPs with CNV segments #############
  if(verbose==TRUE){message("Associate SNPs with CNV segments")}
  segs.pvalue.gr <- .assoPrCNV(all.paths, all.segs.gr, phenotypesSamX, method.m.test, probes.cnv.gr, assign.probe, 
                               correct.inflation=correct.inflation)
  
  ######################################################## END ###########################################
  
  ################################### Run bayesian  INCOMPLETE ######################################################
  ## TODO - Might use to choose the better modelling for each CNV segment (i.e. linear, quadratic etc...)
  
  if(FALSE){
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
  
  map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")), 
                             chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                             position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")), 
                             probes = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
  ))
  
  map.segs <- subset(map, probes %in% values(segs.pvalue.gr)$NameProbe)
  
  #all.samples <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
  gens <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"))) ## Select CNV genotypes
  gens <- t(gens) ## Translate the genotypes
  gens <- gens[, as.numeric(as.character(map.segs$snp.id))] ## Subset segments
  phen <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")))
  gens <- cbind(phen, gens)
  gens[gens[,2]==-9, 2] <- NA
  gens <- gens[,-1]
  colnames(gens)[1] <- "trait"
  gens <- gens[which(gens$trait != "NA"),]
  gens$trait <- as.numeric(as.character(gens$trait))
  #.completeFun <- function(data, desiredCols){
    #completeVec <- complete.cases(data[, desiredCols])
    #return(data[completeVec, ])
  #}
  #gens <- .completeFun(gens, "trait")
   #gens <- subset(gens, !is.na(trait))  
   #dim(gens)
  options(mc.cores = 40)
  mcmc <- list(nChains = 8, burnin = 500, chainLength = 20000,
               thin = 10, reduceRet = FALSE)
  
  colnames(gens)[2:length(gens)] <- paste0("var", seq_len(ncol(gens)-1))
  formB <- formula(paste0("trait~", as.character(paste0(colnames(gens)[2:length(colnames(gens))], collapse="+"))))

  number.of.v0 <- 30
  v0.default <- 2.5e-04
  v0.all <- NULL
  for(lo in seq_len(number.of.v0)){
    v0.all[[lo]] <- v0.default
    v0.default <- v0.default*0.5
  }
  #v0.all
  
  baymod.all <- NULL
  for(lo in seq_len(length(v0.all))){
  hyper <- list(gamma = c(v0 = v0.all[[lo]]), tau2 = c(5,25), w = c(1, 1))
  baymod.all[[lo]] <- spikeSlabGAM::spikeSlabGAM(formB, data=gens, mcmc = mcmc, hyperparameters=hyper)
  }
  
  ## Prepare data to plot
  
  map.segs.all <- list()
  for(lo in seq_len(length(baymod.all))){
  baymod <- baymod.all[[lo]]
  #model.prob <- summary(baymod)$modelTable[[1]]
  #all.p.model <- as.numeric(summary(baymod)$trmSummary[,2])[-1]
  
  map.segs$PBayes <- as.numeric(summary(baymod)$trmSummary[,2])[-1] 
  map.segs$model.P <- summary(baymod)$modelTable[[1]]
  map.segs$v0 <- v0.all[[lo]]
  map.segs$Reg <- "other"
  map.segs$Reg[c(4,5,6)] <- "chr2" 
  map.segs$Reg[c(37,38,39,40)] <- "chr27"
  map.segs.all[[lo]] <- map.segs
  
  }
  
  
  #summary(baymod) 
  
  
  #### Interaction Bayes
  
  interact.bayes <- expand.grid(colnames(gens[,-1]), colnames(gens[,-1]))
  formIntB <- formula(paste0("trait~", as.character(paste0(colnames(gens)[2:length(colnames(gens))], collapse="+")),
                             "+",paste0(interact.bayes$Var1, ":", interact.bayes$Var2, collapse="+")))
  baymod <- spikeSlabGAM::spikeSlabGAM(formIntB, data=gens, mcmc = mcmc)
  
  SNPRelate::snpgdsClose(genofile)
  }
  
  ######################################################## END ###########################################
  
  #################################### Plot the QQ-plot of the analysis
  if(verbose==TRUE){message("Plot the QQ-plot of the analysis")}
  
  pdf(file.path(all.paths[3], paste0(unique(values(segs.pvalue.gr)$Phenotype), "-LRR-", run.lrr, "QQ-PLOT.pdf")))
  qq.plot.pdf <- qqunif.plot(values(segs.pvalue.gr)$MinPvalueAdjusted, auto.key=list(corner=c(.95,.05)))
  print(qq.plot.pdf)
  invisible(dev.off())
  
  ######################################################## END ###########################################
  
  #################################### Reconvert the chrs to original names if applicable
  
  if(!is.null(chr.code.name)){
    (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
    chr.code.name <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Chr.names"))
    
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
  
  #################   Return results or run the epistasis analysis ####################################
  
  segs.pvalue.gr <- segs.pvalue.gr[order(values(segs.pvalue.gr)$MinPvalue)]  
  
  if(epistasis.mode == FALSE){
  return(segs.pvalue.gr)}
  else if(epistasis.mode == TRUE){
  
  if(verbose==TRUE){message("Analysis of CNV pairwise epistasis")}
    
  ## Produce the interaction map
  probes.cnv.gr.gwas <- probes.cnv.gr[which(values(probes.cnv.gr)$Name %in% values(segs.pvalue.gr)$NameProbe)]  # Select probes assigned in the canonical GWAS
  
  
  #index.epistasis <- expand.grid(seq_len(length(probes.cnv.gr.gwas)), seq_len(length(probes.cnv.gr.gwas)))
  index.epistasis <- as.data.frame(t(utils::combn(seq_len(length(probes.cnv.gr.gwas)), 2)), stringsAsFactors=FALSE)
  colnames(index.epistasis) <- c("Var1", "Var2") 
  index.epistasis <- index.epistasis[-which(index.epistasis$Var1 == index.epistasis$Var2),]
  
  ## TODO - CNVR specific
  #if(FALSE){
  if(cnvr.specific){ ## TODO
  ## Make CNVRs
  cnvs <- data.table::fread(file.path(all.paths[1], "CNVOut.txt"), sep = "\t", 
                              header = FALSE)  ### CNV table 
  CNVs <- .checkConvertCNVs(cnvs)
  ####################### Check if the chromosomes are numeric
  chr.names <- CNVs$chr
  chr.names <- gsub("chr", "", chr.names)
  if (any(is.na(chr.names))) {
      stop("Chromosome names should be integers. If they are not, make them integers and include the correspondent names 
           in chr.code.name parameter")
    }
  ############################# 
  CNVs$start <- as.integer(as.character(CNVs$start))
  CNVs$end <- as.integer(as.character(CNVs$end))
  CNVsGr <- makeGRangesFromDataFrame(CNVs, keep.extra.columns = TRUE)
  cnvrs.mom <- reduce(CNVsGr)
  values(cnvrs.mom)$tag <- paste0("cnvr", seq_len(length(cnvrs.mom)))
  
  .subsetCNVR <- function(lo, probes.cnv.gr.gwas, cnvrs.mom){
    return(values(subsetByOverlaps(cnvrs.mom, probes.cnv.gr.gwas[lo]))$tag)
  }
  
  if(rappdirs:::get_os() == "win"){
    param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
    cnvr.tags <- suppressMessages(BiocParallel::bplapply(1:length(probes.cnv.gr.gwas), 
                                                    .subsetCNVR, BPPARAM = param, probes.cnv.gr.gwas=probes.cnv.gr.gwas,
                                                    cnvrs.mom=cnvrs.mom))
  }
  if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
      multicoreParam <- BiocParallel::MulticoreParam(workers = n.cor)
      cnvr.tags <- suppressMessages(BiocParallel::bplapply(1:length(probes.cnv.gr.gwas), 
                                                           .subsetCNVR, BPPARAM = multicoreParam, probes.cnv.gr.gwas=probes.cnv.gr.gwas,
                                                           cnvrs.mom=cnvrs.mom))
    }
  
  
  ## Subset based on the CNVR location   
  values(probes.cnv.gr.gwas)$cnvr <- unlist(cnvr.tags)
  #index.epistasis$cnvr1
  #index.epistasis$cnvr2
  #index.epistasis[index.epistasis$cnvr1==index.epistasis$cnvr2,]
  index.epistasis <- index.epistasis[values(probes.cnv.gr.gwas[index.epistasis$Var1])$cnvr==values(probes.cnv.gr.gwas[index.epistasis$Var2])$cnvr,]
  }
  
  gr.interactions <- InteractionSet::GInteractions(index.epistasis$Var1, index.epistasis$Var2, probes.cnv.gr.gwas)
  values(gr.interactions)$InteractName <- seq_len(length(gr.interactions))
  
  n <- round(length(gr.interactions)/2000)
  #n <- 1
  gr.interactions.split <- split(gr.interactions, 
                                 factor(sort(rank(seq_len(length(gr.interactions)))%%n))) ## Split interactions
  
  pb <- txtProgressBar(min = 0, max = length(gr.interactions.split), style = 3)
  
  #results.all <- list()
  
  for(runP in seq_len(length(gr.interactions.split))){
  if(verbose==TRUE){setTxtProgressBar(pb, runP, title="Running all GWAS analysis for parwise CNV genotypes")}
  
  gr.interactions.splitX <- gr.interactions.split[[runP]]
    
  #if(verbose==TRUE){message("Produce gvar to use as PLINK input (epistasis)")}
  suppressWarnings(.prodPLINKgvar.epistasis(all.paths, n.cor, snp.matrix, run.lrr, gr.interactions=gr.interactions.splitX, 
                           probes.cnv.gr.gwas=probes.cnv.gr.gwas)) ## TODO: Warnings only running the function - why?
    
  #if(verbose==TRUE){message("Produce PLINK map (epistasis)")}
  .prodPLINKmap.epistasis(all.paths, gr.interactions=gr.interactions.splitX)  
    
  #if(verbose==TRUE){message("Run PLINK (epistasis)")}
  .runPLINK(all.paths)
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  sum.file <- paste0("plink.gvar.summary.",names(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")))[[2]])
  results <- read.table(file.path(all.paths[2], sum.file), header=TRUE) 
  #results.all[[runP]] <- results
  SNPRelate::snpgdsClose(genofile)
  
  if(runP==1){
    write.table(results, file.path(all.paths[2],"resultsEpi.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)}
  if(runP>1){
    write.table(results, file.path(all.paths[2],"resultsEpi.txt"), append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
  }
  
  file.rename(file.path(all.paths[2],"resultsEpi.txt"), file.path(all.paths[2], sum.file))
  
  #gentype.all <- data.table::rbindlist(results.all)
  #gentype.all1 <- as.data.frame(gentype.all)
  #write.table(gentype.all1, file.path(all.paths[2],"mydata.gvar"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  
  if(verbose==TRUE){message("Assign probe-interactions and draw p-values for each one")}
  segs.pvalue.gr.epistasis <- .assoPrCNV.epistasis(all.paths, phenotypesSamX, method.m.test, gr.interactions)
  
  if(verbose==TRUE){message("Numeric to character chromosomes, if any")}
  if(!is.null(chr.code.name)){
    (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
    chr.code.name <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Chr.names"))
    
    #segs.pvalue <- data.frame(segs.pvalue.gr)
    #segs.pvalue <- segs.pvalue[ , -which(names(segs.pvalue) %in% c("strand"))]
    #segs.pvalue.gr.epistasis <- renameSeqlevels(segs.pvalue.gr.epistasis, c("16"="1A"))
    
    #seqNew <- seqlevels(segs.pvalue.gr.epistasis)
    for(lopN in seq_len(nrow(chr.code.name))){
      ##seqx <- seqlevels(segs.pvalue.gr.epistasis)[which(match(seqlevels(segs.pvalue.gr.epistasis), 
                                        #chr.code.name[lopN,1], nomatch=0) == 1)] 
      
      #seqlevels(segs.pvalue.gr.epistasis)<- sub(chr.code.name[lopN,1],chr.code.name[lopN,2],seqlevels(segs.pvalue.gr.epistasis))
      #seqNew[which(match(seqNew, 
      #chr.code.name[lopN,1], nomatch=0) == 1)] <-  as.character(chr.code.name[lopN,2])
      #seqx <- as.character(chr.code.name[lopN,2])
      seqlevels(InteractionSet::regions(segs.pvalue.gr.epistasis))[which(match(seqlevels(InteractionSet::regions(segs.pvalue.gr.epistasis)), 
                                                      chr.code.name[lopN,1], nomatch=0) == 1)] <- as.character(chr.code.name[lopN,2])
      #segs.pvalue$seqnames <- gsub(chr.code.name[lopN,1],chr.code.name[lopN,2], segs.pvalue$seqnames)
      
     #c(seqlevels(segs.pvalue.gr.epistasis) = gsub(paste0("\\<",chr.code.name[lopN,1], "\\>"),as.character(chr.code.name[lopN,2]), seqlevels(segs.pvalue.gr.epistasis)))
      
     # seqlevels(segs.pvalue.gr.epistasis) <- gsub(paste0("\\<",chr.code.name[lopN,1], "\\>"), 
      #                                            as.character(chr.code.name[lopN,2]), seqlevels(segs.pvalue.gr.epistasis))
                                               #seqlevels(segs.pvalue.gr), fixed=TRUE)
      #seqlevels(segs.pvalue.gr)
      #seqlevels(segs.pvalue.gr.epistasis, force=TRUE)[which(match(seqlevels(segs.pvalue.gr.epistasis), 
       #                                               chr.code.name[lopN,1], nomatch=0) == 1)] <- gsub("1", "X", seqlevels(segs.pvalue.gr), fixed=TRUE)
    }
    #seqlevels(segs.pvalue.gr.epistasis) <- as.factor(seqNew)
    #segs.pvalue.gr <- makeGRangesFromDataFrame(segs.pvalue, keep.extra.columns = TRUE)
    SNPRelate::snpgdsClose(genofile)
  }
  
  segs.and.interaction <- list()
  
  segs.and.interaction[[1]] <- segs.pvalue.gr
  segs.and.interaction[[2]] <- segs.pvalue.gr.epistasis
  
  return(segs.and.interaction)
  
  }
  ######################################################## END ########################################
}


# HELPER - Produce the PLINK .gvar file in disk - epistasis
# 

.prodPLINKgvar.epistasis <- function(all.paths, n.cor, snp.matrix, run.lrr=FALSE, gr.interactions, probes.cnv.gr.gwas){
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE)) ## read GDS
  sam.gen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")) ## Extract samples
  #snps <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")) ## Extract SNPs
  snps <- as.numeric(as.character(values(probes.cnv.gr.gwas)$snp.id))
  fam.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "FamID")) ## Extract family
  
  ## Produce gvar file in parallel processing
  # Identify the SO
  if(run.lrr == TRUE){
    if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
      multicoreParam <- BiocParallel::MulticoreParam(workers = n.cor)
      Gens <- suppressMessages(BiocParallel::bplapply(1:length(sam.gen), .prodGvarLRR.epistasis, BPPARAM = multicoreParam, genofile=genofile, sam.gen=sam.gen,
                                                      snps=snps, fam.id=fam.id, snp.matrix=snp.matrix, gr.interactions=gr.interactions))
    }
    
    if(rappdirs:::get_os() == "win"){
      param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
      Gens <- suppressMessages(BiocParallel::bplapply(1:length(sam.gen), .prodGvarLRR.epistasis, BPPARAM = param, genofile=genofile, sam.gen=sam.gen,
                                                      snps=snps, fam.id=fam.id, snp.matrix=snp.matrix, gr.interactions=gr.interactions))
    }
    
  }else{
    if(rappdirs:::get_os() == "unix" | rappdirs:::get_os() == "mac"){
      multicoreParam <- BiocParallel::MulticoreParam(workers = n.cor)
      Gens <- suppressMessages(BiocParallel::bplapply(1:length(sam.gen), .prodGvar.epistasis, BPPARAM = multicoreParam, genofile=genofile, sam.gen=sam.gen,
                                                      snps=snps, fam.id=fam.id, snp.matrix=snp.matrix, gr.interactions=gr.interactions))
    }
    
    if(rappdirs:::get_os() == "win"){
      param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
      Gens <- suppressMessages(BiocParallel::bplapply(1:length(sam.gen), .prodGvar.epistasis, BPPARAM = param, genofile=genofile, sam.gen=sam.gen,
                                                      snps=snps, fam.id=fam.id, snp.matrix=snp.matrix, gr.interactions=gr.interactions))
    }
    
  }
  
  #gentype.all <- data.table::rbindlist(Gens)
  #gentype.all1 <- as.data.frame(gentype.all)
  #write.table(gentype.all1, file.path(all.paths[2],"mydata.gvar"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  
  
  #fileConn<-file(file.path(all.paths[2],"mydata.gvar"), open="r")
  
  #pb <- txtProgressBar(min = 0, max = length(Gens), style = 3)
  
  for(wr in seq_len(length(Gens))){
  #message(wr)
  #lapply(Gens, write.table, fileConn, append=TRUE, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  GenX <- as.data.frame(Gens[[wr]])
  if(wr==1){
  write.table(GenX, file.path(all.paths[2],"mydata.gvar"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)}
  else{
  write.table(GenX, file.path(all.paths[2],"mydata.gvar"), append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
  #setTxtProgressBar(pb, wr, title="Writing gvar file for all interactions")
  }
  
  
  #close(fileConn)
  SNPRelate::snpgdsClose(genofile)
} ## TODO - modify - LRR TODO

# HELPER - Produce PLINK probe map in disk - epistasis
# 

.prodPLINKmap.epistasis <- function(all.paths, gr.interactions){
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE))
  #Pmap <- as.data.frame(cbind(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
  #                            gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
  #                            0, gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position"))))
  
  chr.x <- as.character(seqnames(InteractionSet::anchors(gr.interactions)$first))
  chr.y <- as.character(seqnames(InteractionSet::anchors(gr.interactions)$second))
  chr.both <- paste0(chr.x,"-", chr.y)
  #head(chr.both)
  #snps.interactions <- seq_len(length(gr.interactions)) ## BUG!
  snps.interactions <- values(gr.interactions)$InteractName
  positions.both <- paste0(start(InteractionSet::anchors(gr.interactions)$first),"-",
                           start(InteractionSet::anchors(gr.interactions)$second))
  #positions.both
  
  Pmap <- as.data.frame(cbind(chr.both, snps.interactions, 0, positions.both))
  
  colnames(Pmap) <- c("Chr", "NAME", "GD", "Position")
  write.table(Pmap, file.path(all.paths[2], "mydata.map"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  SNPRelate::snpgdsClose(genofile)
}

# HELPER - Internal function of parallel implementation to produce the .gvar file (absolute copy number) 
# requested for the CNV-GWAS in PLINK - epistasis
# @examples
# Gens <- .prodGvar(lo, genofile, sam.gen, snps, fam.id)

.prodGvar.epistasis <- function(lo, genofile, sam.gen, snps, fam.id, snp.matrix, gr.interactions){
  g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id=sam.gen[[lo]], verbose=FALSE))
  #A <- rep("A", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
  #B <- rep("B", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))))
  CNVg <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), start=c(1,lo), count=c(length(g),1))
  
  ###
  CNVg <- CNVg[snps] #Subset to have only SNPs used in the canonical GWAS
  index.x <- as.numeric(InteractionSet::anchors(gr.interactions, id=TRUE)$first) ## Use anchor to make the interaction gvar object
  index.y <- as.numeric(InteractionSet::anchors(gr.interactions, id=TRUE)$second)
  Dose1 <- CNVg[index.x]
  Dose2 <- CNVg[index.y]
  
  #snps.interactions <- seq_len(length(gr.interactions))
  snps.interactions <- values(gr.interactions)$InteractName
  A <- rep("A", length(snps.interactions))
  B <- rep("B", length(snps.interactions))
  
  #CNVg <- CNVg-1
  #Dose1 <- rep(1, length(CNVg))
  #Dose1[CNVg==-1] <- 0
  #CNVg[CNVg==-1] <- 0
  #B[CNVg==1] <- "A"
  ###
  
  if(snp.matrix=="FALSE"){
    df <- as.data.frame(cbind(as.character(fam.id[[lo]]), sam.gen[[lo]], snps.interactions, 
                              #A, CNVg, B, Dose1))
                              A, Dose1, B, Dose2))
    colnames(df) <- c("FID", "IID", "NAME", "ALLELE1", "DOSAGE1", "ALLELE2", "DOSAGE2")}
  
  if(snp.matrix=="TRUE"){ #### TODO - To implement
    ...
  }
  return(df)
} 

# HELPER - Internal function of parallel implementation to produce the .gvar file (LRR) requested for the CNV-GWAS in PLINK - epistasis
# @examples
# Gens <- .prodGvarLRR(lo, genofile, sam.gen, snps, fam.id)

.prodGvarLRR.epistasis <- function(lo, genofile, sam.gen, snps, fam.id, snp.matrix){
  
  ## Transform LRR calclulations into genotypes
  g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id=sam.gen[[lo]], verbose=FALSE))
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
} ## TODO - modify


# HELPER - Associate tested probe-interactions and draw p-values/adjusted p-values for each interaction
# 

.assoPrCNV.epistasis  <- function(all.paths, phenotypesSamX, method.m.test, gr.interactions){
  
  #segs.pvalue.gr <- unlist(all.segs.gr)
  #values(segs.pvalue.gr)$SegName <- seq(from=1, to=length(segs.pvalue.gr), by=1)
  
  (genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE, readonly=FALSE))
  
  sum.file <- paste0("plink.gvar.summary.",names(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype")))[[2]])
  resultsp <- read.table(file.path(all.paths[2], sum.file), header=TRUE) 
  resultsp.back <- resultsp
  
  resultsp <- resultsp[ grep("P\\(CNP\\)", resultsp$FIELD),]
  resultsp$VALUE <- as.numeric(as.character(resultsp$VALUE))
  
  resultsp <- resultsp[match(values(gr.interactions)$InteractName, resultsp$NAME),]
  #resultspGr <- makeGRangesFromDataFrame(resultsp, seqnames.field = "Chr",
  #                                       start.field = "Position", end.field = "Position",
  #                                       keep.extra.columns = TRUE)
  
  values(gr.interactions)$CNP <- resultsp$VALUE
  values(gr.interactions)$CNP_Padjusted <- stats::p.adjust(values(gr.interactions)$CNP, 
                                                           method = method.m.test, n = length(gr.interactions))
  
  resultsp <- resultsp.back
  
  ### Include "P(CNP|SNP)" p.values
  resultsp.other.p <- resultsp[ which(resultsp$FIELD %in% c("P(CNP|SNP)")),] 
  resultsp.other.p <- resultsp.other.p[match(values(gr.interactions)$InteractName, resultsp.other.p$NAME),]
  resultsp.other.p <- as.numeric(as.character(resultsp.other.p$VALUE))
  resultsp.other.p[is.na(resultsp.other.p)] <- "no_test"
  values(gr.interactions)$CNP_SNP <- resultsp.other.p
  
  #p.values.other <- rep("no_test", length(gr.interactions))
  #p.values.other[which(values(gr.interactions)$InteractName %in% resultsp.other.p$NAME)] <- as.numeric(as.character(resultsp.other.p$VALUE))
  #values(resultspGr)$CNP_SNP <-  p.values.other
  
  ### Include "P(SNP|CNP)" p.values 
  #resultsp.other.p <- resultsp[ which(resultsp$FIELD %in% c("P(SNP|CNP)")),] 
  resultsp.other.p <- resultsp[ which(resultsp$FIELD %in% c("P(SNP|CNP)")),] 
  resultsp.other.p <- resultsp.other.p[match(values(gr.interactions)$InteractName, resultsp.other.p$NAME),]
  resultsp.other.p <- as.numeric(as.character(resultsp.other.p$VALUE))
  resultsp.other.p[is.na(resultsp.other.p)] <- "no_test"
  values(gr.interactions)$SNP_CNP <- resultsp.other.p
  
  ### Include "P(SNP&CNP)" p.values 
  #resultsp.other.p <- resultsp[ which(resultsp$FIELD %in% c("P(SNP&CNP)")),] 
  resultsp.other.p <- resultsp[ which(resultsp$FIELD %in% c("P(SNP&CNP)")),] 
  resultsp.other.p <- resultsp.other.p[match(values(gr.interactions)$InteractName, resultsp.other.p$NAME),]
  resultsp.other.p <- as.numeric(as.character(resultsp.other.p$VALUE))
  resultsp.other.p[is.na(resultsp.other.p)] <- "no_test"
  values(gr.interactions)$SNPeCNP <- resultsp.other.p
  
  
  #resultspGr <- resultspGr[which(values(resultspGr)$NAME %in% values(segs.pvalue.gr)$NameProbe)] # EXCLUDE?
  #resultspGr <- resultspGr[match(values(segs.pvalue.gr)$NameProbe, values(resultspGr)$NAME)] # EXCLUDE?
  #resultspGr[mcols(resultspGr)$CNP_SNP=="no_test"] # EXCLUDE?
  
  ### Include "P(CNP|SNP)" adjusted p.values 
  num.of.tests <- length(values(gr.interactions)$CNP_SNP[values(gr.interactions)$CNP_SNP != "no_test"])
  corr.pvalues <- suppressWarnings(stats::p.adjust(values(gr.interactions)$CNP_SNP, 
                                                   method = method.m.test, n = num.of.tests))
  corr.pvalues <- round(corr.pvalues, 4)
  corr.pvalues[is.na(corr.pvalues)] <- "no_test"
  
  values(gr.interactions)$CNP_SNP_Padjusted <- corr.pvalues
  
  
  ### Include "P(SNP|CNP)" adjusted p.values 
  num.of.tests <- length(values(gr.interactions)$SNP_CNP[values(gr.interactions)$SNP_CNP != "no_test"])
  corr.pvalues <- suppressWarnings(stats::p.adjust(values(gr.interactions)$SNP_CNP, 
                                                   method = method.m.test, n = num.of.tests))
  corr.pvalues <- round(corr.pvalues, 4)
  corr.pvalues[is.na(corr.pvalues)] <- "no_test"
  
  values(gr.interactions)$SNP_CNP_Padjusted <- corr.pvalues
  
  
  ### Include "P(SNP&CNP)" adjusted p.values 
  num.of.tests <- length(values(gr.interactions)$SNPeCNP[values(gr.interactions)$SNPeCNP != "no_test"])
  corr.pvalues <- suppressWarnings(stats::p.adjust(values(gr.interactions)$SNPeCNP, 
                                                   method = method.m.test, n = num.of.tests))
  corr.pvalues <- round(corr.pvalues, 4)
  corr.pvalues[is.na(corr.pvalues)] <- "no_test"
  
  values(gr.interactions)$CNPeSNP_Padjusted <- corr.pvalues
  gr.interactions
  
  SNPRelate::snpgdsClose(genofile)
  
  return(gr.interactions)
} 


# HELPER Create the pairwise genotypes for epistasis analysis - TODO
# 


# HELPER - Debug snippet - Dont't Run
# 

run.tests <- function(){
index.epistasis <- base::expand.grid(seq_len(length(probes.cnv.gr)), seq_len(length(probes.cnv.gr)))
old.row <- nrow(index.epistasis)
#index.epistasis$tag <- paste0(index.epistasis$Var1, "-", index.epistasis$Var2)
#index.epistasis <- index.epistasis[-which(duplicated(index.epistasis$tag)), ]
index.epistasis <- index.epistasis[-which(index.epistasis$Var1 == index.epistasis$Var2),]
old.row
nrow(index.epistasis)

gi <- GInteractions(index.epistasis$Var1, index.epistasis$Var2, probes.cnv.gr)
#gi.new <- swapAnchors(gi)




segs.pvalue.gr <- cnvGWAS.V2(phen.info, n.cor=1, min.sim = 0.9, freq.cn = 0.05, snp.matrix=FALSE, 
                          method.m.test = "fdr", lo.phe=1, chr.code.name=chr.code.name,
                          genotype.nodes=c("CNVGenotype"),
                          produce.gds=TRUE, verbose=TRUE, assign.probe = "min.pvalue", #"min.pvalue", "high.freq"
                          coding.translate= NULL, epistasis.mode=FALSE, exclude.outliers=FALSE,
                          minq=.01, maxq=.99, log.pheno=FALSE, both.up.down=FALSE)


### CNVR specific
segs.pvalue.gr <- cnvGWAS.V2(phen.info, n.cor=1, min.sim = 0.95, freq.cn = 0.05, snp.matrix=FALSE, 
                             method.m.test = "fdr", lo.phe=1, chr.code.name=chr.code.name,
                             genotype.nodes=c("CNVGenotype"),
                             produce.gds=FALSE, verbose=TRUE, assign.probe = "min.pvalue", #"min.pvalue", "high.freq"
                             coding.translate= NULL, epistasis.mode =TRUE, cnvr.specific=TRUE, BothUpDown=TRUE)

(genofile <- SNPRelate::snpgdsOpen(file.path(all.paths[1], "CNV.gds"), allow.fork=TRUE))

}