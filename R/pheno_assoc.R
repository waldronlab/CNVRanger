################## 
# Author: Vinicius Henrique da Silva 
# Script description: Functions related to CNV-GWAS
# Date: July 20, 2018 
# Code: CNVASSOPACK002
################## 

#' Run the CNV-GWAS 
#'
#' Wraps all the necessary functions to run a CNV-GWAS using the output of 
#' \code{\link{setupCnvGWAS}} function
#'  
#' (i) Produces the GDS file containing the genotype information  (if produce.gds == TRUE),
#' (ii) Produces the requested inputs for a PLINK analysis, 
#' (iii) run a CNV-GWAS analysis using linear model implemented in PLINK 
#' (http://zzz.bwh.harvard.edu/plink/gvar.shtml), and 
#' (iv) export a QQ-plot displaying the adjusted p-values. 
#' In this release only the p-value for the copy number is available  (i.e. 'P(CNP)'). 
#'   
#' @param phen.info Returned by \code{\link{setupCnvGWAS}}
#' @param n.cor Number of cores to be used
#' @param min.sim Minimum CNV genotype distribution similarity among subsequent
#' probes. Default is 0.95 (i.e. 95\%)
#' @param freq.cn Minimum CNV frequency where 1 (i.e. 100\%), or all samples 
#' deviating from diploid state. Default 0.01 (i.e. 1\%)
#' @param snp.matrix Only FALSE implemented - If TRUE B allele frequencies (BAF)
#' would be used to reconstruct CNV-SNP genotypes
#' @param method.m.test Correction for multiple tests to be used. FDR is default, 
#' see \code{\link{p.adjust}} for other methods.
#' @param lo.phe The phenotype to be analyzed in the PhenInfo$phenotypesSam data-frame 
#' @param chr.code.name A data-frame with the integer name in the first column 
#' and the original name for each chromosome  
#' @param genotype.nodes Expression data type. Nodes with CNV genotypes to be 
#' produced in the gds file. 
#' @param coding.translate For 'CNVgenotypeSNPlike'. If NULL or unrecognized 
#' string use only biallelic CNVs. If 'all' code multiallelic CNVs as 0 
#' for loss; 1 for 2n and 2 for gain.
#' @param path.files Folder containing the input CNV files used for the CNV 
#' calling (i.e. one text file with 5 collumns for each sample). Columns should 
#' contain (i) probe name, (ii) Chromosome, (iii) Position, (iv) LRR, and (v) BAF.
#' @param list.of.files Data-frame with two columns where the (i) is the file 
#' name with signals and (ii) is the correspondent name of the sample in the gds file
#' @param produce.gds logical. If TRUE produce a new gds, if FALSE use gds previously created   
#' @param run.lrr If TRUE use LRR values instead absolute copy numbers in the association
#' @param assign.probe \sQuote{min.pvalue} or \sQuote{high.freq} to represent the CNV segment
#' @param correct.inflation logical. Estimate lambda from raw p-values and correct for genomic inflation.
#' Use with argument \code{method.m.test} to generate strict p-values. 
#' @param both.up.down Check for CNV genotype similarity in both directions. 
#' Default is FALSE (i.e. only downstream)
#' @param verbose Show progress in the analysis
#' @return The CNV segments and the representative probes and their respective p-value
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @seealso \code{link{setupCnvGWAS}} to setup files needed for the CNV-GWAS.
#' @references da Silva et al. (2016) Genome-wide detection of CNVs and their 
#' association with meat tenderness in Nelore cattle. PLoS One, 11(6):e0157711.
#'
#' @examples
#' 
#' # Load phenotype-CNV information
#' data.dir <- system.file("extdata", package="CNVRanger")
#' 
#' phen.loc <- file.path(data.dir, "Pheno.txt")
#' cnv.out.loc <- file.path(data.dir, "CNVOut.txt")
#' map.loc <- file.path(data.dir, "MapPenn.txt")
#'
#' phen.info <- setupCnvGWAS('Example', phen.loc, cnv.out.loc, map.loc)
#' 
#' # Define chr correspondence to numeric, if necessary
#' df <- '16 1A
#' 25 4A
#' 29 25LG1
#' 30 25LG2
#' 31 LGE22'
#' 
#' chr.code.name <- read.table(text=df, header=FALSE)
#' segs.pvalue.gr <- cnvGWAS(phen.info, chr.code.name=chr.code.name)
#'  
#' @export

cnvGWAS <- function(phen.info, n.cor = 1, min.sim = 0.95, freq.cn = 0.01, snp.matrix = FALSE, 
    method.m.test = "fdr", lo.phe = 1, chr.code.name = NULL, genotype.nodes = "CNVGenotype", 
    coding.translate = "all", path.files = NULL, list.of.files = NULL, produce.gds = TRUE, 
    run.lrr = FALSE, assign.probe = "min.pvalue", correct.inflation = FALSE, both.up.down = FALSE, 
    verbose = FALSE) {
    
    phenotypesSam <- phen.info$phenotypesSam
    phenotypesSamX <- phenotypesSam[, c(1, (lo.phe + 1))]
    phenotypesSamX <- stats::na.omit(phenotypesSamX)
    all.paths <- phen.info$all.paths
    
    # Produce the GDS for a given phenotype
    if (produce.gds) {
        if (verbose) 
            message("Produce the GDS for a given phenotype")
        probes.cnv.gr <- generateGDS(phen.info = phen.info, freq.cn = freq.cn, snp.matrix = snp.matrix, 
            lo.phe = lo.phe, chr.code.name = chr.code.name, genotype.nodes = genotype.nodes, 
            coding.translate = coding.translate)
    } else {
        if (verbose) 
            message("Using existent gds file")
        probes.cnv.gr <- .prodProbes(phen.info, lo.phe, freq.cn)
    }
    
    # Produce PLINK map 
    if (verbose) 
        message("Produce PLINK map")
    .prodPLINKmap(all.paths)
    
    # Produce gvar to use as PLINK input
    if (verbose) 
        message("Produce gvar to use as PLINK input")
    .prodPLINKgvar(all.paths, n.cor, snp.matrix, run.lrr)
    
    # Produce fam (phenotype) to use as PLINK input
    if (verbose) 
        message("Produce fam (phenotype) to use as PLINK input")
    .prodPLINKfam(all.paths)
    
    # Run PLINK
    if (verbose) 
        message("Run PLINK")
    suppressMessages(.runPLINK(all.paths))
    
    # Produce CNV segments
    if (verbose) 
        message("Produce CNV segments")
    all.segs.gr <- .prodCNVseg(all.paths, probes.cnv.gr, min.sim, both.up.down)
    
    # Associate SNPs with CNV segments
    if (verbose) 
        message("Associate SNPs with CNV segments")
    segs.pvalue.gr <- .assoPrCNV(all.paths, all.segs.gr, phenotypesSamX, method.m.test, 
        probes.cnv.gr, assign.probe, correct.inflation = correct.inflation)
    
    # Plot the QQ-plot of the analysis
    if (verbose) 
        message("Plot the QQ-plot of the analysis")
   
    uphen <- unique(segs.pvalue.gr$Phenotype)
    qq.plot.pdf <- paste0(uphen, "-LRR-", run.lrr, "QQ-PLOT.pdf")
    qq.plot.pdf <- file.path(all.paths["Results"], qq.plot.pdf)
    pdf(qq.plot.pdf)
    print(.qqunifPlot(segs.pvalue.gr$MinPvalueAdjusted, 
                        auto.key = list(corner = c(0.95, 0.05))))
    invisible(dev.off())
    
    # Reconvert the chrs to original names if applicable
    if (!is.null(chr.code.name)) {
        cnv.gds <- file.path(all.paths["Inputs"], "CNV.gds")
        genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork = TRUE, readonly = FALSE)
        chr.code.name <- gdsfmt::index.gdsn(genofile, "Chr.names")
        chr.code.name <- gdsfmt::read.gdsn(chr.code.name)
        
        segs.pvalue <- data.frame(segs.pvalue.gr)
        segs.pvalue <- segs.pvalue[, !(names(segs.pvalue) %in% "strand")]

        segs.pvalue$seqnames <- as.vector(segs.pvalue$seqnames)
        cmap <- as.vector(chr.code.name[, 2])
        names(cmap) <- as.vector(chr.code.name[, 1])
        ind <- segs.pvalue$seqnames %in% names(cmap)
        segs.pvalue$seqnames[ind] <- cmap[segs.pvalue$seqnames[ind]]       

        segs.pvalue.gr <- GenomicRanges::makeGRangesFromDataFrame(segs.pvalue, keep.extra.columns = TRUE)
        SNPRelate::snpgdsClose(genofile)
    }
    
    segs.pvalue.gr <- segs.pvalue.gr[order(segs.pvalue.gr$MinPvalue)]
    return(segs.pvalue.gr)
}


#' Setup the folders and files to run CNV-GWAS analysis
#'
#' This function creates the (i) necessary folders in disk to perform downstream 
#' analysis on CNV genome-wide association and (ii) import the necessary input 
#' files (i.e. phenotypes, probe map and CNV list) from other locations in disk.
#' 
#' The user can import several phenotypes at once. All information will be 
#' stored in the list returned by this function. 
#' The user should be aware although several phenotypes can be imported, the 
#' \code{\link{cnvGWAS}} or \code{\link{generateGDS}} functions will handle only 
#' one phenotype per run. 
#' 
#' @param name String with a project code or name (e.g. 'Project1')
#' @param phen.loc Path/paths to the tab separated text file containing phenotype 
#' and sample info. When using more than one population, for populations without
#' phenotypes include the string 'INEXISTENT' instead the path for a file.
#' @param cnv.out.loc Path(s) to the CNV analysis output (i.e. PennCNV output, 
#' SNP-chip general format or sequencing general format). It is also possible to
#' use a \code{\linkS4class{RaggedExperiment}} or a \code{\linkS4class{GRangesList}} 
#' object instead if the run includes only one population.
#' @param map.loc Path to the probe map (e.g. used in PennCNV analysis). Column 
#' names containing probe name, chromosome and coordinate must be named as: Name, 
#' Chr and Position. Tab delimited. If NULL, artificial probes will be generated 
#' based on the CNV breakpoints.
#' @param folder Choose manually the project folder (i.e. path as the root folder). 
#' Otherwise, user-specific data dir will be used automatically.  
#' @param pops.names Indicate the name of the populations, if using more than one.
#' @param n.cor Number of cores
#' @return List \sQuote{phen.info} with \sQuote{samplesPhen}, \sQuote{phenotypes}, 
#' \sQuote{phenotypesdf}, \sQuote{phenotypesSam}, \sQuote{FamID}, \sQuote{SexIds}, 
#' \sQuote{pops.names} (if more than one population) and \sQuote{all.paths}
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples
#' 
#' data.dir <- system.file("extdata", package="CNVRanger")
#' 
#' phen.loc <- file.path(data.dir, "Pheno.txt")
#' cnv.out.loc <- file.path(data.dir, "CNVOut.txt")
#' map.loc <- file.path(data.dir, "MapPenn.txt")
#'
#' phen.info <- setupCnvGWAS('Example', phen.loc, cnv.out.loc, map.loc)
#' 
#' 
#' @export
setupCnvGWAS <- function(name, phen.loc, cnv.out.loc, 
    map.loc = NULL, folder = NULL, pops.names = NULL, n.cor = 1) 
{  
    ## Create the folder structure for all subsequent analysis
    all.paths <- .createFolderTree(name, folder)  
  
    if(is(cnv.out.loc, "RaggedExperiment"))
    {
        phen.info <- colData(cnv.out.loc)
        stopifnot(all(names(phen.info)[1:3] == c("sample.id", "fam", "sex")))
        cnv.out.loc <- as(cnv.out.loc, "GRangesList")
      
        phen.loc <- file.path(all.paths["Inputs"], "PhenN.txt")
        write.table(as.data.frame(phen.info), 
            file=phen.loc, sep="\t", row.names=FALSE, quote=FALSE)
    }
    
    ## Check if the CNVs are in GRanges format   
    if(is(cnv.out.loc, "GRangesList"))
    {
        df.cnv <- as.data.frame(cnv.out.loc)
        if(all(c("num.snps", "start.probe", "end.probe") %in% colnames(df.cnv)))
        { 
            rel.cols <- c("seqnames", "start", "end", 
                "group_name", "state", "num.snps", "start.probe", "end.probe")
            df.cnv <- df.cnv[,rel.cols]
            colnames(df.cnv)[c(1,4)] <- c("chr", "sample.id")
        }
        else df.cnv <- df.cnv[,c("seqnames", "start", "end", "group_name", "state")]
        colnames(df.cnv)[c(1,4)] <- c("chr", "sample.id")
        
        cnv.out.loc <- file.path(all.paths["Inputs"], "CNVOutN.txt")
        write.table(df.cnv, file=cnv.out.loc, sep="\t", row.names=FALSE, quote=FALSE)
    }   
    
    ## Only one population
    if (length(phen.loc) == 1 && length(cnv.out.loc) == 1) 
    {
        ## Import the phenotype and sample info from external folder
        file.copy(phen.loc, file.path(all.paths["Inputs"], "/PhenoPop1.txt"), overwrite = TRUE)
        ## Import the PennCNV output from external folder
        file.copy(cnv.out.loc, file.path(all.paths["Inputs"], "/CNVOut.txt"), overwrite = TRUE)
        file.copy(cnv.out.loc, file.path(all.paths["Inputs"], "/CNVOutPop1.txt"), overwrite = TRUE)
        pheno.file <- "PhenoPop1.txt"
        cnv.file <- "CNVOut.txt"
    }
    
    ## Multiple populations
    else if (length(phen.loc) != length(cnv.out.loc)) 
        stop("phen.loc and cnv.out.loc should have the same length. Use the string INEXISTENT if phenotypes are missing")
      
    else if (length(phen.loc) > 1 && length(cnv.out.loc) > 1)
    {
        grid <- seq_along(phen.loc)
        pheno.files <- paste0("/PhenoPop", grid, ".txt")
        cnv.files <- paste0("/CNVOutPop", grid, ".txt")
        
        abs.pheno.files <- file.path(all.paths["Inputs"], pheno.files)
        for (npop in grid) 
        {
            if (phen.loc[npop] == "INEXISTENT") cat("INEXISTENT", file=abs.pheno.files[npop])
            else file.copy(phen.loc[npop], abs.pheno.files[npop], overwrite = TRUE)
        }
        file.copy(cnv.out.loc, file.path(all.paths["Inputs"], cnv.files), overwrite = TRUE)
       
        ## Write phenotype names and merge CNV file for multiple populations
        pheno.file <- pheno.files
        all.cnvs <- lapply(cnv.files, .loadToMergeCNV, cnv.path = all.paths["Inputs"])
        all.cnvs <- data.table::rbindlist(all.cnvs)
        all.cnvs <- as.data.frame(all.cnvs)
        write.table(all.cnvs, file.path(all.paths["Inputs"], "CNVOut.txt"), sep = "\t", 
            col.names = TRUE, row.names = FALSE)
    }
    
    if (is.null(map.loc)) {
        ## Import the probe map from external folder
        cnvs <- read.table(file.path(all.paths["Inputs"], "CNVOut.txt"), sep = "", header = FALSE)  ### CNV table 
        CNVs <- .checkConvertCNVs(cnvs, all.paths, n.cor)
        CGr <- GenomicRanges::makeGRangesFromDataFrame(CNVs)
        
    } else {
        if (length(map.loc) > 1) 
            stop("The map should be unique. If multiple populations with different probe map use map.loc=NULL")
        file.copy(map.loc, file.path(all.paths["Inputs"], "MapPenn.txt"), overwrite = TRUE)
    }
    
    phen.info <- .loadPhen(pheno.file, all.paths, pops.names = pops.names, n.cor)
    plink.dir <- dir(all.paths["PLINK"], pattern = "plink-1*")
    
    if (!length(plink.dir)) {
        got.plink <- .getPLINK(all.paths["PLINK"])
        if (!got.plink) stop("PLINK setup failed")
        plink.dir <- dir(all.paths["PLINK"], pattern = "plink-1*")
    }
    plink.bin <- file.path(all.paths["PLINK"], plink.dir, "plink")
    all.paths["PLINK"] <- dirname(plink.bin)
    
    phen.info$all.paths <- all.paths
    phen.info$phenotypesSam$samplesPhen <- as.character(phen.info$phenotypesSam$samplesPhen)
    
    ## Delete temporary files
    cnv.out.loc <- file.path(all.paths["Inputs"], "CNVOutN.txt")
    if(file.exists(cnv.out.loc)) file.remove(cnv.out.loc)
    
    phen.loc <- file.path(all.paths["Inputs"], "PhenN.txt")
    if(file.exists(phen.loc)) file.remove(phen.loc)
    
    return(phen.info)
}


#' Produce CNV-GDS for the phenotyped samples
#'
#' Function to produce the GDS file in a probe-wise fashion for CNV genotypes. 
#' The GDS file which is produced also incorporates one phenotype to be analyzed. 
#' If several phenotypes are enclosed in the \sQuote{phen.info} object, the user 
#' may specify the phenotype to be analyzed with the \sQuote{lo.phe} parameter. 
#' Only diploid chromosomes should be included.  
#'
#' @param phen.info Returned by \code{setupCnvGWAS}
#' @param freq.cn Minimum frequency. Default is 0.01 (i.e. 1\%)
#' @param snp.matrix Only FALSE implemented. If TRUE, B allele frequencies (BAF) 
#' and SNP genotypes would be used to reconstruct CNV-SNP genotypes - under development
#' @param lo.phe The phenotype to be analyzed in the PhenInfo$phenotypesSam dataframe 
#' @param chr.code.name A data-frame with the integer name in the first column 
#' and the original name in the second for each chromosome previously converted to numeric
#' @param genotype.nodes Nodes with CNV genotypes to be produced in the gds file. 
#' Use 'CNVGenotype' for dosage-like genotypes (i.e. from 0 to Inf). 
#' Use 'CNVgenotypeSNPlike' alongside for SNP-like CNV genotype in a separated 
#' node (i.e.  '0, 1, 2, 3, 4' as '0/0, 0/1, 1/1, 1/2, 2/2').
#' @param coding.translate For 'CNVgenotypeSNPlike'. If NULL or unrecognized 
#' string use only biallelic CNVs. If 'all' code multiallelic CNVs as 0 for loss; 
#' 1 for 2n and 2 for gain. 
#' @param n.cor Number of cores
#' @return probes.cnv.gr Object with information about all probes to be used in 
#' the downstream CNV-GWAS. Only numeric chromosomes
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples
#' 
#' # Load phenotype-CNV information
#' data.dir <- system.file("extdata", package="CNVRanger")
#' 
#' phen.loc <- file.path(data.dir, "Pheno.txt")
#' cnv.out.loc <- file.path(data.dir, "CNVOut.txt")
#' map.loc <- file.path(data.dir, "MapPenn.txt")
#'
#' phen.info <- setupCnvGWAS('Example', phen.loc, cnv.out.loc, map.loc)
#'
#' # Construct the data-frame with integer and original chromosome names 
#'  
#' # Define chr correspondence to numeric, if necessary
#' df <- '16 1A
#' 25 4A
#' 29 25LG1
#' 30 25LG2
#' 31 LGE22'
#' 
#' chr.code.name <- read.table(text=df, header=FALSE)
#' probes.cnv.gr <- generateGDS(phen.info, chr.code.name=chr.code.name)
#' 
#'@export

generateGDS <- function(phen.info, freq.cn = 0.01, snp.matrix = FALSE, lo.phe = 1,
    chr.code.name = NULL, genotype.nodes = c("CNVGenotype", "CNVgenotypeSNPlike"), 
    coding.translate = NULL, n.cor = 1) 
{
    phenotypesSam <- phen.info$phenotypesSam
    samplesPhen <- phen.info$samplesPhen
    FamID <- phen.info$FamID
    SexIds <- phen.info$SexIds
    all.paths <- phen.info$all.paths
    phenotypesSamX <- phenotypesSam[, c(1, (lo.phe + 1))]
    
    # Import CNVs to data-frame
    cnv.table <- file.path(all.paths["Inputs"], "CNVOut.txt")
    cnvs <- data.table::fread(cnv.table, sep = "\t", header = FALSE)
    CNVs <- .checkConvertCNVs(cnvs, all.paths, n.cor)
    
    # Check if the chromosomes are numeric
    chr.names <- CNVs$chr
    chr.names <- gsub("chr", "", chr.names)
    
    stopifnot(all(!is.na(chr.names)))
    
    CNVsGr <- GenomicRanges::makeGRangesFromDataFrame(CNVs, keep.extra.columns = TRUE)
    # Subset CNVs in phenotyped samples
    CNVsGr <- CNVsGr[CNVsGr$V5 %in% samplesPhen]
    
    # Import SNP map to data-frame
    map.file <- file.path(all.paths["Inputs"], "MapPenn.txt")
    probes <- data.table::fread(map.file, header = TRUE, sep = "\t")
    probes <- as.data.frame(probes)
    probes <- probes[stats::complete.cases(probes), ]
    probesGr <- GenomicRanges::makeGRangesFromDataFrame(probes, seqnames.field = "Chr", 
        start.field = "Position", end.field = "Position", keep.extra.columns = TRUE)
    
    all.samples <- unique(CNVsGr$V5)
    
    # Select probes within CNVs
    probesCNV <- suppressWarnings(IRanges::subsetByOverlaps(probesGr, CNVsGr)$Name)
    probesCNV <- unique(unlist(probesCNV))
    probes.cnv.gr <- probesGr[probesGr$Name %in% probesCNV]
    
    counts <- GenomicRanges::countOverlaps(probes.cnv.gr, CNVsGr)
    probes.cnv.gr$freq <- unname(counts)
    
    # Subset by frequency
    NumSam <- freq.cn * length(all.samples)
    probes.cnv.gr <- probes.cnv.gr[probes.cnv.gr$freq >= NumSam]
    
    # Order the probes
    probes.cnv.gr <- GenomeInfoDb::sortSeqlevels(probes.cnv.gr)
    probes.cnv.gr <- GenomicRanges::sort(probes.cnv.gr)
    
    probes.cnv.gr$snp.id <- seq_along(probes.cnv.gr)
    
    # SNP genotype matrix not available
    if (!snp.matrix) CNVBiMa <- matrix(2, nrow = length(probes.cnv.gr), ncol = length(all.samples)) 
    else stop("Option to consider SNP matrix is not implemented yet")
    # TODO: SNP genotype matrix available 
    # CHECK IF WE HAVE THE SAME SAMPLES -
    # INCLUDE NA FOR SAMPLES WITH ONLY ONE GENOTYPE TYPE? (i.e CNV or SNP)
    
    # Create a GDS with chr and SNP names to numeric
    cnv.gds <- file.path(all.paths["Inputs"], "CNV.gds")
    SNPRelate::snpgdsCreateGeno(cnv.gds, genmat = CNVBiMa, sample.id = all.samples, 
        snp.id = as.character(probes.cnv.gr$snp.id), snp.rs.id = probes.cnv.gr$Name, 
        snp.chromosome = as.character(GenomicRanges::seqnames(probes.cnv.gr)), 
        snp.position = GenomicRanges::start(probes.cnv.gr))
    

    genotype.nodes <- match.arg(genotype.nodes)    
    genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork = TRUE, readonly = FALSE)
    if (genotype.nodes == "CNVGenotype") {
        # Replace genotype matrix with CNV genotypes
        n <- gdsfmt::add.gdsn(genofile, "CNVgenotype", CNVBiMa, replace = TRUE)
        
        param <- BiocParallel::SnowParam(workers = 1, type = "SOCK")
        BiocParallel::bplapply(seq_along(all.samples), .writeProbesCNV, BPPARAM = param, 
                all.samples = all.samples, genofile = genofile, CNVsGr = CNVsGr, 
                probes.cnv.gr = probes.cnv.gr, n = n)
    }
    else {
        # CNVgenotypeSNPlike
        CNVBiMaCN <- matrix(2, nrow = length(probes.cnv.gr), ncol = length(all.samples))
        n <- gdsfmt::add.gdsn(genofile, "CNVgenotypeSNPlike", CNVBiMaCN, replace = TRUE)
        
        #Define the chunks of SNPs to import chunk
        chunk <- seq(from = 1, to = length(probes.cnv.gr), by = length(probes.cnv.gr))  ## No chunk
        
        if (chunk[length(chunk)] != length(probes.cnv.gr)) {
            chunk[length(chunk) + 1] <- length(probes.cnv.gr)
        }
        
        os <- rappdirs:::get_os()
        if (os %in% c("unix", "mac")) param <- BiocParallel::MulticoreParam(workers = 1)
        else param <- BiocParallel::SnowParam(workers = 1, type = "SOCK")
        BiocParallel::bplapply(seq_len(length(chunk) - 1), .recodeCNVgenotype, BPPARAM = param, 
                genofile = genofile, all.samples = all.samples, n = n, chunk = chunk, 
                coding.translate = coding.translate)
    }
    
    # Include the phenotype 'lo' in the GDS
    phenotypesSamX <- phenotypesSamX[match(all.samples, phenotypesSamX$samplesPhen),]  
    gdsfmt::add.gdsn(genofile, name = "phenotype", val = phenotypesSamX, replace = TRUE)
    FamID <- FamID[match(all.samples, FamID$samplesPhen), ]
    gdsfmt::add.gdsn(genofile, name = "FamID", val = as.character(FamID[, 2]), replace = TRUE, 
        storage = "string")
    FamID <- SexIds[match(all.samples, SexIds$samplesPhen), ]
    gdsfmt::add.gdsn(genofile, name = "Sex", val = as.character(SexIds[, 2]), replace = TRUE, 
        storage = "string")
    
    # Include chr names if needed
    if (!is.null(chr.code.name))
        gdsfmt::add.gdsn(genofile, name = "Chr.names", val = chr.code.name, replace = TRUE)
    
    if (!is.null(phen.info$pops.names))
        gdsfmt::add.gdsn(genofile, name = "Pops.names", val = phen.info$pops.names, 
            replace = TRUE)
    
    SNPRelate::snpgdsClose(genofile)
    return(probes.cnv.gr)
}


#' Import LRR and BAF from text files used in the CNV analysis
#' 
#' This function imports the LRR/BAF values and create a node for each one in 
#' the GDS file at the working folder 'Inputs' created by the 
#' \code{\link{setupCnvGWAS}} function. Once imported, the LRR values can be 
#' used to perform a GWAS directly as an alternative to copy number dosage
#'
#' @param all.paths Object returned from \code{CreateFolderTree} function with 
#' the working folder tree 
#' @param path.files Folder containing the input CNV files used for the CNV 
#' calling (i.e. one text file with 5 collumns for each sample). Columns should 
#' contain (i) probe name, (ii) Chromosome, (iii) Position, (iv) LRR, and (v) BAF.
#' @param list.of.files Data-frame with two columns where the (i) is the file 
#' name with signals and (ii) is the correspondent name of the sample in the gds file
#' @param gds.file Path to the GDS file which contains nodes harboring respective LRR and 
#' BAF values. The \sQuote{snp.rs.id}, \sQuote{sample.id}, \sQuote{LRR} and \sQuote{BAF}
#' nodes are mandatory. Both the SNPs and samples should follow the order and length in the CNV.gds
#' (located at all.paths["Inputs"] folder). \sQuote{path.files} and \sQuote{list.of.files}
#' will be ignored if \sQuote{gds.file} is not NULL
#' @param verbose Print the samples while importing
#' @return Writes to the specified GDS file by side effect.
#' @author Vinicius Henrique da Silva <vinicius.dasilva@@wur.nl>
#' @examples
#' 
#' # Load phenotype-CNV information
#' data.dir <- system.file("extdata", package="CNVRanger")
#' 
#' phen.loc <- file.path(data.dir, "Pheno.txt")
#' cnv.out.loc <- file.path(data.dir, "CNVOut.txt")
#' map.loc <- file.path(data.dir, "MapPenn.txt")
#'
#' phen.info <- setupCnvGWAS('Example', phen.loc, cnv.out.loc, map.loc)
#'
#' # Extract path names
#' all.paths <- phen.info$all.paths
#' 
#' # List files to import LRR/BAF 
#' list.of.files <- list.files(path=data.dir, pattern="cnv.txt.adjusted$")
#' list.of.files <- as.data.frame(list.of.files)
#' colnames(list.of.files)[1] <- "file.names"
#' list.of.files$sample.names <- sub(".cnv.txt.adjusted$", "", list.of.files$file.names)
#' 
#' # All missing samples will have LRR = '0' and BAF = '0.5' in all SNPs listed in the GDS file
#' importLrrBaf(all.paths, data.dir, list.of.files)
#' 
#' # Read the GDS to check if the LRR/BAF nodes were added
#' cnv.gds <- file.path(all.paths["Inputs"], 'CNV.gds')    
#' genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork=TRUE, readonly=FALSE)
#' SNPRelate::snpgdsClose(genofile)
#' 
#' @export

importLrrBaf <- function(all.paths, path.files, list.of.files, gds.file=NULL, verbose=TRUE){
    
    cnv.gds <- file.path(all.paths["Inputs"], "CNV.gds")
    genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork=TRUE, readonly=FALSE)
    
    # Create file path for all text files
    file.paths <- file.path(path.files, list.of.files[, 1])
    gds.names <- list.of.files[, 2]
    list.filesLo <- data.frame(add = file.paths, gds = gds.names)
    
    snps.included <- gdsfmt::index.gdsn(genofile, "snp.rs.id")
    snps.included <- gdsfmt::read.gdsn(snps.included)
    snps.included <- as.character(snps.included)
   
    all.samples <- gdsfmt::index.gdsn(genofile, "sample.id")
    all.samples <- gdsfmt::read.gdsn(all.samples)
    all.samples <- as.character(all.samples)
    
    if(!is.null(gds.file)){
        lrr.baf.gds <- SNPRelate::snpgdsOpen(gds.file, allow.fork=TRUE, readonly=FALSE)
        nodes.ava <- GDSArray::gdsnodes(lrr.baf.gds)
        
        stopifnot(all(c("snp.rs.id", "sample.id", "LRR", "BAF") %in% nodes.ava))
        
        snps.included.signals <- gdsfmt::index.gdsn(lrr.baf.gds, "snp.rs.id")
        snps.included.signals <- gdsfmt::read.gdsn(snps.included.signals)
        snps.included.signals <- as.character(snps.included.signals)
        
        stopifnot(length(snps.included.signals)==length(snps.included))
        stopifnot(all(snps.included.signals==snps.included))
        
        all.samples.signals <- gdsfmt::index.gdsn(lrr.baf.gds, "sample.id")
        all.samples.signals <- gdsfmt::read.gdsn(all.samples.signals)
        all.samples.signals <- as.character(all.samples.signals)
        
        stopifnot(length(all.samples.signals)==length(all.samples))
        stopifnot(all(all.samples.signals==all.samples))
        
        LRR.matrix <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(lrr.baf.gds, "LRR"))
        BAF.matrix <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(lrr.baf.gds, "BAF"))
        
        gdsfmt::add.gdsn(lrr.baf.gds, "LRR", LRR.matrix, replace = TRUE)
        gdsfmt::add.gdsn(lrr.baf.gds, "BAF", BAF.matrix, replace = TRUE)
        
        SNPRelate::snpgdsClose(lrr.baf.gds)
    }
    else
    {
        # Check if all files
        list.filesLo.back <- list.filesLo
        list.filesLo <- list.filesLo[list.filesLo$gds %in% all.samples, ]
        if(verbose)
        {
            if (nrow(list.filesLo) != nrow(list.filesLo.back))
            warning("list.of.files has different length of the list of samples from gds")
        }
        
        LRR.matrix <- matrix(0, nrow = length(snps.included), ncol = length(all.samples))
        nLRR <- gdsfmt::add.gdsn(genofile, "LRR", LRR.matrix, replace = TRUE)
        gdsfmt::read.gdsn(nLRR)
        
        BAF.matrix <- matrix(0.5, nrow = length(snps.included), ncol = length(all.samples))
        nBAF <- gdsfmt::add.gdsn(genofile, "BAF", BAF.matrix, replace = TRUE)
        gdsfmt::read.gdsn(nBAF)
        
        if (verbose) message("Start parallel import of LRR/BAF values")
        param <- BiocParallel::SnowParam(workers = 1, type = "SOCK")
        BiocParallel::bplapply(seq_len(nrow(list.filesLo)), .freadImport, BPPARAM = param, 
            list.filesLo = list.filesLo, genofile = genofile, all.samples = all.samples, 
            nLRR = nLRR, nBAF = nBAF, snps.included = snps.included, verbose = verbose)
    }
    
    SNPRelate::snpgdsClose(genofile)
}

# Extract the CNV genotypes from the GDS file @param genofile is the loaded gds
# file @param snp.id is the SNP probe name to retrieve @param node.to.extract
# \sQuote{CNVgenotype} or \sQuote{CNVgenotypeSNPlike}. Default is
# \sQuote{CNVgenotype} @return CNV genotypes in a specific probe @author
# Vinicius Henrique da Silva <vinicius.dasilva@wur.nl>
.snpgdsGetGenoCNV <- function(genofile, snp.id, node.to.extract = "CNVgenotype") 
{
    map <- data.frame(
        snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")), 
        chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
        position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")), 
        probes = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")), 
        stringsAsFactors = FALSE)
    
    all.samples <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
    snp.id <- as.numeric(as.character(map[map$probes == snp.id]))
    gens <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, node.to.extract), 
        start = c(snp.id, 1), count = c(1, length(all.samples)))
    
    return(gens)
}

# HELPER - Create the folder tree to keep necessary files during and after the
# analysis @param name String with a project code or name (e.g. 'Project1')
# @param folder Choose manually the project folder (i.e. path as the root
# folder). If NULL, a standard program folder will be chosen.  @return List with
# paths placed to store the files produced by subsequent analysis @examples
# all.paths <- createFolderTree('Project_name')
.createFolderTree <- function(name, folder = NULL) 
{
    # data dir
    if (is.null(folder)) 
        folder <- rappdirs::user_data_dir("CNVRanger")
    if (!file.exists(folder)) 
        dir.create(folder, recursive = TRUE)
    
    # project directory
    proj.dir <- file.path(folder, name)
    if (!file.exists(proj.dir)) dir.create(proj.dir)
    
    # subdirs
    dir.names <- c("Inputs", "PLINK", "Results")
    dirs <- dir.names
    all.paths <- file.path(proj.dir, dirs)
    for (d in all.paths) if (!file.exists(d)) dir.create(d)
    names(all.paths) <- dir.names
    
    return(all.paths)
}

# HELPER - Download and test PLINK 1.07 @param all.paths Object returned from
# \code{CreateFolderTree} function with the working folder tree @param version
# PLINK version. Only 1.07 implemented @return boolean. Success (TRUE) or fail
# (FALSE) in running PLINK
.getPLINK <- function(plink.path, version = "1.07") 
{
    plink.url <- "http://zzz.bwh.harvard.edu/plink/dist/plink-"
    plink.url <- paste0(plink.url, version, "-")
    plink.file <- file.path(plink.path, "PLINK.zip")
    
    # get os-specific binary
    os <- rappdirs:::get_os()
    os <- switch(os, unix = "x86_64", mac = "mac-intel", win = "dos", os)
    plink.url <- paste0(plink.url, os, ".zip")
    
    # download & unzip
    download.file(plink.url, plink.file)
    unzip(plink.file, exdir = plink.path)
    
    # remove zip
    file.remove(plink.file)
    
    # path to plink binary
    plink.file <- dir(plink.path, pattern = "plink*")
    plink.path <- file.path(plink.path, plink.file, "plink")
    Sys.chmod(plink.path, mode = "755")
    
    suppressWarnings(res <- system2(plink.path, args = "--noweb", stdout = TRUE, 
        stderr = FALSE))
    
    res <- any(grepl(version, res))
    return(res)
}

# HELPER - Load phenotypes for each of the samples to be analyzed @param file.nam
# Name of the file to be imported @param pheno.path First item of the list
# returned by \code{CreateFolderTree} function @param pops.names Indicate the
# name of the populations. Only used if more than one population exists @param
# n.cor Number of cores @return List \sQuote{phen.info} with
# \sQuote{samplesPhen}, \sQuote{phenotypes}, \sQuote{phenotypesdf},
# \sQuote{phenotypesSam}, \sQuote{FamID} and \sQuote{SexIds} and
# \sQuote{pops.names} (if more than one population) @examples phen.info <-
# loadPhen('/home/.../Phen.txt') sapply(phen.info, class)
.loadPhen <- function(file.nam, all.paths, pops.names = NULL, n.cor = 1) 
{
    pheno.path <- all.paths["Inputs"]
    samplesPhen.all <- NULL
    phenotypes.all <- NULL
    phenotypesdf.all <- NULL
    phenotypesSam.all <- NULL
    FamID.all <- NULL
    SexIds.all <- NULL
    
    for (npop in seq_along(file.nam)) {
        dfN <- data.table::fread(file.path(pheno.path, file.nam[npop]))
        cnv.file <- paste0("CNVOutPop", npop, ".txt")
        cnvs <- read.table(file.path(pheno.path, cnv.file), sep = "", header = FALSE)
        CNVs <- .checkConvertCNVs(cnvs, all.paths, n.cor)
        all.samples <- as.character(CNVs$V5)
        all.samples <- unique(all.samples)

        if (all(names(dfN) == "INEXISTENT")) {
            ## Extract the sample names from CNV file instead
            samplesPhen <- all.samples
            phenotypesSam <- data.frame(all.samples, "NA", stringsAsFactors=FALSE)
            colnames(phenotypesSam) <- c("samplesPhen", phenotypes.all[[1]])
            FamID <- SexIds <- phenotypesSam
            colnames(FamID)[2] <- colnames(SexIds)[2] <- "V2"
        } 
        else {
            all <- data.table::fread(file.path(pheno.path, file.nam[npop]), header = TRUE, 
                sep = "\t")  ### Import phenotypes
            all <- as.data.frame(all)
            stopifnot(all(colnames(all)[1:3] == c("sample.id", "fam", "sex")))
            
            ### Include samples with CNV but without phenotypes
            pheno.na <- data.frame(sample.id=all.samples, fam = NA, sex = NA, stringsAsFactors=FALSE)
            pheno.na <- pheno.na[!(pheno.na$sample.id %in% all$sample.id), ]
            
            all <- plyr::rbind.fill(all, pheno.na)
            all[is.na(all)] <- -9
            
            ### Produce phen.info objects
            samplesPhen <- unique(all$sample.id)  ## Greb sample names
            phenotypesdf <- all[, 4:ncol(all), drop = FALSE]
            phenotypes <- names(phenotypesdf)  ## Greb phenotype names
            phenotypesdf <- as.data.frame(phenotypesdf)  ### Greb phenotypes values
            
            phenotypesSam <- data.frame(samplesPhen, phenotypesdf, stringsAsFactors=FALSE)
            ind <- as.numeric(rownames(phenotypesSam))
            FamID <- data.frame(samplesPhen, V2 = all$fam[ind])
            SexIds <- data.frame(samplesPhen, V2 = all$sex[ind])
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
    ind <- colnames(phenotypesSam) == phenotypes
    phenotypesSam[, ind] <- suppressWarnings(as.numeric(as.character(phenotypesSam[,ind])))
    phenotypesSam[is.na(phenotypesSam)] <- -9
    FamID <- data.table::rbindlist(FamID.all)
    FamID <- as.data.frame(FamID)
    SexIds <- data.table::rbindlist(SexIds.all)
    SexIds <- as.data.frame(SexIds)
    
    phen.info <- list(samplesPhen = samplesPhen, phenotypes = phenotypes, phenotypesdf = phenotypesdf, 
            phenotypesSam = phenotypesSam, FamID = FamID, SexIds = SexIds)

    if (!is.null(pops.names)) 
    {
        grid <- seq_along(file.nam)
        phen.info$pops.names <- rep(pops.names[grid], each=lengths(samplesPhen.all[grid]))    
    }
    return(phen.info)
}

# HELPER - Check and convert CNV input @param cnvs Data-frame with the CNVs to be
# analyzed. From (i) PennCNV, (ii) SNP-chip general format or (iii) sequencing
# general format
.checkConvertCNVs <- function(cnvs, all.paths, n.cor = 1) 
{
    .convertPenn <- function(cnvs)
    { 
        cnvs <- as.data.frame(cnvs)
        cnvs$start <- as.integer(cnvs$start)
        cnvs$end <- as.integer(cnvs$end)
        cnvs$V1 <- paste0(cnvs$chr, ":", cnvs$start, "-", cnvs$end)
        cnvs$length <- (cnvs$end - cnvs$start) + 1
        cnvs$length <- paste0("length=", cnvs$length)
        cnvs$num.snps <- paste0("numsnp=", cnvs$num.snps)
        cnvs$start.probe <- paste0("startSNP=", cnvs$start.probe)
        cnvs$end.probe <- paste0("endSNP=", cnvs$end.probe)
        
        cnvs <- cnvs[, c(stand.names[1:3], "V1", "num.snps", "length", "state", "sample.id", 
            "start.probe", "end.probe")]
        colnames(cnvs)[5:10] <- paste0("V", 2:7)
        return(cnvs)
    }

    the.names <- as.character(as.matrix(cnvs[1, ]))
    stand.names <- c("chr", "start", "end", "sample.id", "state")
    stand.names.ext <- c(stand.names, "num.snps", "start.probe", "end.probe")
    ## If PennCNV file
    if (unique(sub("^([[:alpha:]]*).*", "\\1", cnvs$V2))[[1]] == "numsnp") 
    {
        cnvs <- as.data.frame(cnvs)
        df1 <- reshape2::colsplit(cnvs$V1, ":", c("chr", "loc"))
        df2 <- reshape2::colsplit(df1$loc, "-", c("start", "end"))
        df1 <- cbind(df1[, -2, drop = FALSE], df2)
        
        cnvs <- cbind(df1, cnvs)
        cnvs$V1 <- gsub("chr", "", cnvs$V1)
        colnames(cnvs)[1:3] <- c("chr", "start", "end")
    } 
    ## If general format
    else if (suppressWarnings(all(the.names == stand.names.ext))) 
    {
        cnvs <- cnvs[-1, ]
        colnames(cnvs) <- the.names
        cnvs <- .convertPenn(cnvs)  
    } 
    ## If sequencing info
    else if (all(the.names == stand.names)) 
    {
        cnvs <- as.data.frame(cnvs)
        ### Include artificial probe tags
        cnvs.seq <- cnvs
        cnvs.seq <- cnvs.seq[-1, ]
        colnames(cnvs.seq) <- the.names
        cnv.seq.gr <- GenomicRanges::makeGRangesFromDataFrame(cnvs.seq, keep.extra.columns = TRUE)
        chr.seq <- as.character(GenomicRanges::seqnames(cnv.seq.gr))
        start.seq <- GenomicRanges::start(cnv.seq.gr)
        start.seq <- data.frame(chr = chr.seq, position = start.seq, stringsAsFactors = FALSE)
        end.seq <- GenomicRanges::end(cnv.seq.gr)
        end.seq <- data.frame(chr = chr.seq, position = end.seq, stringsAsFactors = FALSE)
        dif <- GenomicRanges::end(cnv.seq.gr) - GenomicRanges::start(cnv.seq.gr)
        middle.seq <- GenomicRanges::end(cnv.seq.gr) - dif
        middle.seq <- data.frame(chr = chr.seq, position = middle.seq, stringsAsFactors = FALSE)
        arti.pr <- rbind(start.seq, end.seq, middle.seq)  # Bind all artifical probes
        arti.pr <- GenomicRanges::makeGRangesFromDataFrame(arti.pr, seqnames.field = "chr", 
            start.field = "position", end.field = "position")
        arti.pr <- GenomicRanges::reduce(arti.pr)  ## Exclude duplicated positions
        arti.pr <- GenomeInfoDb::sortSeqlevels(arti.pr)
        arti.pr <- GenomicRanges::sort(arti.pr)
        
        probe.like.map <- as.data.frame(arti.pr)
        probe.like.map$Name <- paste0("probe.like_", seq_len(nrow(probe.like.map)))
        probe.like.map <- probe.like.map[, c("Name", "seqnames", "start")]
        colnames(probe.like.map) <- c("Name", "Chr", "Position")
        probe.like.map$Chr <- gsub("chr", "", probe.like.map$Chr)
        write.table(probe.like.map, file.path(all.paths["Inputs"], "MapPenn.txt"), row.names = FALSE, 
            quote = FALSE, sep = "\t")
        
        ## Associate artificial probes to CNVs
        probe.like.map.gr <- GenomicRanges::makeGRangesFromDataFrame(probe.like.map, 
            seqnames.field = "Chr", start.field = "Position", end.field = "Position", 
            keep.extra.columns = TRUE)
        
        ### HELPER - associate probes-like regions with CNVs
        .probeToCNVs <- function(lo, cnv.seq.gr, probe.like.map.gr) {
            ov.pr <- IRanges::subsetByOverlaps(probe.like.map.gr, cnv.seq.gr[lo])
            num.snps <- length(ov.pr)
            start.probe <- ov.pr$Name[1]
            end.probe <- ov.pr$Name[num.snps]
            return(c(num.snps, start.probe, end.probe))
        }
        
        os <- rappdirs:::get_os()
        if (os == "win") param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")
        else param <- BiocParallel::MulticoreParam(workers = n.cor)
        cnvs.probes <- suppressMessages(
            BiocParallel::bplapply(seq_along(cnv.seq.gr), 
            .probeToCNVs, BPPARAM = param, 
            cnv.seq.gr = cnv.seq.gr, probe.like.map.gr = probe.like.map.gr))
      
        cnv.p.df <- t(as.data.frame(cnvs.probes))
        cnv.seq.gr$num.snps <- as.integer(as.character(cnv.p.df[, 1]))
        cnv.seq.gr$start.probe <- as.character(cnv.p.df[, 2])
        cnv.seq.gr$end.probe <- as.character(cnv.p.df[, 3])
        
        cnv.seq <- as.data.frame(cnv.seq.gr)
        cnv.seq <- cnv.seq[, c(1:3, 6:10)]
        colnames(cnv.seq)[1] <- "chr"
        
        cnvs <- cnv.seq
        cnvs <- .convertPenn(cnvs)
    } 
    else stop("Unexpected CNV input format - is it tab delimited?")
    
    return(cnvs)
}

# HELPER - Load multiple CNV inputs @param cnv.file.all CNV file name @param
# cnv.path CNV path name
.loadToMergeCNV <- function(cnv.file.all, cnv.path) {
    ## Load and merge CNV files from multiple populations
    data.table::fread(file.path(cnv.path, cnv.file.all), header = TRUE)
}

# HELPER - Write CNV-genotype in a dosage-like fashion
.writeProbesCNV <- function(lo, all.samples, genofile, CNVsGr, probes.cnv.gr, n) {
    sampleX <- all.samples[lo]
    g <- SNPRelate::snpgdsGetGeno(genofile, sample.id = sampleX, verbose = FALSE)
    g <- drop(g)
    
    ind <- CNVsGr$V5 == sampleX
    CNVSamByState <- split(CNVsGr[ind], CNVsGr[ind]$V4)
    
    for (slo in seq_along(CNVSamByState)) {
        curr <- CNVSamByState[[slo]]
        curr4 <- curr$V4
        state <- unique(curr4)
        
        SamCNVSNP <- suppressWarnings(IRanges::subsetByOverlaps(probes.cnv.gr, curr))
        SamCNVSNP <- SamCNVSNP$snp.id
        g[SamCNVSNP] <- state
    }
    
    g <- as.integer(g)
    gdsfmt::write.gdsn(n, g, start = c(1, lo), count = c(length(g), 1))
}

# HELPER - Write CNV-genotype in a SNP-like fashion @param lo loop number @param
# genofile loaded gds file @param n is the object from \code{gdsfmt::add.gdsn}
# @param ranges.gr is the CNVs of each sample
.replaceSNPtoCNV <- function(lo, all.samples, genofile, CNVsGr, probes.cnv.gr, n) 
{
    sampleX <- all.samples[[lo]]
    g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id = sampleX, verbose = FALSE))
    g <- 1
    
    GenomeInfoDb::seqlevels(CNVsGr) <- GenomeInfoDb::seqlevels(probes.cnv.gr)
    CNVsGr <- CNVsGr[CNVsGr$sample == sampleX]
    states <- paste0("state", c(1,2,5,6), ",cn=", c(0,1,3,4))
    code <- rep(c(0,2), each=2)
    for(i in 1:4)
    {
        s <- states[i]
        co <- code[i]        
        cnvs <- CNVsGr[CNVsGr$type == s]
        cnv.gr <- IRanges::subsetByOverlaps(probes.cnv.gr, cnvs)
        g[cnv.gr$tag.snp] <- co
    }    

    ### Replace the genotype in the gds file
    gdsfmt::write.gdsn(n, g, start = c(1, lo), count = c(length(g), 1))
}

# HELPER - Write CNV-genotype in SNP-like fashion in biallelic loci @param lo
# loop number, based on the number of SNPs per turn @param genofile loaded gds
# file @param n is the object from \code{gdsfmt::add.gdsn} @param chunk
# Intervals of SNP to import and convert @param coding.translate For
# 'CNVgenotypeSNPlike'. If NULL or unrecognized string use only biallelic CNVs.
# If 'all' code multiallelic CNVs as 0 for loss; 1 for 2n and 2 for gain.
.recodeCNVgenotype <- function(lo, genofile, all.samples, n, chunk, coding.translate) 
{
    count.end <- (chunk[lo + 1]) - chunk[lo]
    add <- ifelse(length(chunk) - 1 != lo, 0, 1)
    CNVgenoX <- (g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotypeSNPlike"), 
            start = c(chunk[lo], 1), count = c(count.end + add, length(all.samples))))
    
    ### Put four copies as max
    CNVgenoX[CNVgenoX > 4] <- 4
    CNVgenoX <- as.data.frame(CNVgenoX)
    
    states <- c("z", "o", "t", "th", "f")
    states <- paste0(states, "n")
    CNVgenoX[,states] <- rep(0:4, each=nrow(CNVgenoX))
    
    ### Count genotypes
    genos <- apply(CNVgenoX, 1, table)
    genos <- genos - 1
    genos <- t(genos)
    
    #### Recode genotypes
    sum12 <- rowSums(genos[, 1:2])   
    sum45 <- rowSums(genos[, 4:5])
    indexLoss <- (sum12 > 0) & (sum45 == 0)
    indexGain <- (sum12 == 0) & (sum45 > 0)
    indexExclude <- (sum12 > 0) & (sum45 > 0)
    
    ## Exclude index in the matrix
    CNVgenoX <- CNVgenoX[, -((ncol(CNVgenoX) - 4):ncol(CNVgenoX))]
    
    ### To character
    CNVgenoX <- t(apply(CNVgenoX, 1, paste0, "n"))
    
    ### Replace Gain
    CNVgenoX[indexGain,] <- CNVgenoX[indexGain,] - 2
    
    ### Replace Loss
    CNVgenoX[indexLoss, ] <- 2 - CNVgenoX[indexLoss, ]
    
    ## Non bi-allelic CNVs = no genotype
    if (is.null(coding.translate)) coding.translate <- "biallelic"
    else if (coding.translate == "all") 
    {
        gmap <- c(0, 0, 1, 2, 2)
        # just for orientation: names(gmap) <- c("0n", "1n", "2n", "3n", "4n")
        for(i in indexExclude) CNVgenoX[i,] <- gmap[CNVgenoX[i,] + 1] 
    } 
    else CNVgenoX[indexExclude, ] <- -1
    
    ### Replace the genotype in the gds file
    gdsfmt::write.gdsn(n, CNVgenoX, 
                        start = c(chunk[lo], 1),
                        count = c(count.end + add, 
                        length(all.samples)))
}

# HELPER Wait x seconds # from
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/Sys.sleep.html
testit <- function(x) {
    p1 <- proc.time()
    Sys.sleep(x)
    proc.time() - p1
}

# HELPER - Produce the probes.cnv.gr
.prodProbes <- function(phen.info, lo.phe = 1, freq.cn = 0.01) {
    
    phenotypesSam <- phen.info$phenotypesSam
    samplesPhen <- phen.info$samplesPhen
    FamID <- phen.info$FamID
    SexIds <- phen.info$SexIds
    all.paths <- phen.info$all.paths
    phenotypesSamX <- phenotypesSam[, c(1, (lo.phe + 1))]
    
    ###################### Import CNVs to data-frame
   
    cnv.file <- file.path(all.paths["Inputs"], "CNVOut.txt")
    cnvs <- data.table::fread(cnv.file, sep = "\t", header = FALSE)  ### CNV table 
    CNVs <- .checkConvertCNVs(cnvs, all.paths)
    
    ####################### Check if the chromosomes are numeric
    chr.names <- gsub("chr", "", CNVs$chr)
    stopifnot(all(!is.na(chr.names)))     
    
    CNVs$start <- as.integer(as.character(CNVs$start))
    CNVs$end <- as.integer(as.character(CNVs$end))
    CNVsGr <- GenomicRanges::makeGRangesFromDataFrame(CNVs, keep.extra.columns = TRUE)
    CNVsGr <- CNVsGr[CNVsGr$V5 %in% samplesPhen]  ### Subset CNVs in phenotyped samples
    
    ###################### Import SNP map to data-frame
    map.file <- file.path(all.paths["Inputs"], "MapPenn.txt")
    probes <- data.table::fread(map.file, header = TRUE, sep = "\t")
    probes <- as.data.frame(probes)
    probes$Position <- as.integer(as.character(probes$Position))
    probes <- probes[stats::complete.cases(probes), ]
    probesGr <- GenomicRanges::makeGRangesFromDataFrame(probes, seqnames.field = "Chr", 
        start.field = "Position", end.field = "Position", keep.extra.columns = TRUE)
    all.samples <- unique(as.character(CNVsGr$V5))
    
    ###################### Select probes within CNVs
    GenomeInfoDb::seqlevels(CNVsGr) <- GenomeInfoDb::seqlevels(probesGr)
    probesCNV <- suppressWarnings(as.character(IRanges::subsetByOverlaps(probesGr, 
        CNVsGr)$Name))
    probesCNV <- unique(unlist(probesCNV))
    probes.cnv.gr <- probesGr[probesGr$Name %in% probesCNV]
    counts <- GenomicRanges::countOverlaps(probes.cnv.gr, CNVsGr)
    probes.cnv.gr$freq <- unname(counts)
    
    ##### Subset by frequency
    NumSam <- freq.cn * length(all.samples)
    probes.cnv.gr <- probes.cnv.gr[probes.cnv.gr$freq >= NumSam]
    
    ##### Order the probes
    probes.cnv.gr <- GenomeInfoDb::sortSeqlevels(probes.cnv.gr)
    probes.cnv.gr <- GenomicRanges::sort(probes.cnv.gr)
    probes.cnv.gr$snp.id <- seq_along(probes.cnv.gr)
    
    return(probes.cnv.gr)
}


# HELPER - Internal function of parallel implementation to produce the .gvar file
# (absolute copy number) requested for the CNV-GWAS in PLINK @examples Gens <-
# .prodGvar(lo, genofile, sam.gen, snps, fam.id)
.prodGvar <- function(lo, genofile, sam.gen, snps, fam.id, snp.matrix) {
    g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id = sam.gen[[lo]], 
        verbose = FALSE))
    start <-  c(1, lo)
    count <- c(length(g), 1)
    CNVg <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), 
            start = start, count = count)
    CNVg <- CNVg - 1
    len <- length(CNVg)
    A <- rep("A", len)
    B <- rep("B", len)
    Dose1 <- rep(1, len)
    ind <- CNVg == -1
    Dose1[ind] <- CNVg[ind] <- 0
    B[CNVg == 1] <- "A"
    
    if (!snp.matrix) {
        df <- data.frame(fam.id[[lo]], sam.gen[[lo]], snps, A, CNVg, B, Dose1)
        colnames(df) <- c("FID", "IID", "NAME", "ALLELE1", "DOSAGE1", "ALLELE2", 
            "DOSAGE2")
    } else stop("Option to consider SNP matrix is not implemented yet") 
   
    return(df)
}

# HELPER - Internal function of parallel implementation to produce the .gvar file
# (LRR) requested for the CNV-GWAS in PLINK @examples Gens <- .prodGvarLRR(lo,
# genofile, sam.gen, snps, fam.id)
.prodGvarLRR <- function(lo, genofile, sam.gen, snps, fam.id, snp.matrix) 
{
    ## Transform LRR calclulations into genotypes
    g <- as.numeric(SNPRelate::snpgdsGetGeno(genofile, sample.id = sam.gen[[lo]], 
        verbose = FALSE))
    A <- rep("A", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), 
        start = c(1, lo), count = c(length(g), 1))))
    B <- rep("B", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"), 
        start = c(1, lo), count = c(length(g), 1))))
    LRRg <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "LRR"), start = c(1, lo), 
        count = c(length(g), 1))
    LRRg[is.na(LRRg)] <- 0
    
    CNVg <- LRRg
    Dose1 <- rep(0, length(CNVg))
    
    
    if (!snp.matrix) {
        df <- as.data.frame(cbind(as.character(fam.id[[lo]]), sam.gen[[lo]], snps, 
            A, CNVg, B, Dose1))
        colnames(df) <- c("FID", "IID", "NAME", "ALLELE1", "DOSAGE1", "ALLELE2", 
            "DOSAGE2")
    } else {
        stop("Option to consider SNP matrix is not implemented yet")  
        ## CHECK IF WE HAVE THE SAME SAMPLES
        ## INCLUDE NA FOR SAMPLES WITH ONLY ONE GENOTYPE TYPE? (i.e CNV or SNP)
    }
    return(df)
}

# HELPER - Produce PLINK probe map in disk
.prodPLINKmap <- function(all.paths) 
{
    cnv.gds <- file.path(all.paths["Inputs"], "CNV.gds")
    genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork = TRUE)
    Pmap <- data.frame(
            gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
            gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")), 
            0, gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")))
    
    colnames(Pmap) <- c("Chr", "NAME", "GD", "Position")
    map.file <- file.path(all.paths["PLINK"], "mydata.map")
    write.table(Pmap, map.file, sep = "\t", row.names = FALSE, quote = FALSE)
    SNPRelate::snpgdsClose(genofile)
}

# HELPER - Produce the PLINK .gvar file in disk
.prodPLINKgvar <- function(all.paths, n.cor, snp.matrix, run.lrr = FALSE) 
{
    cnv.gds <- file.path(all.paths["Inputs"], "CNV.gds")
    genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork = TRUE, readonly = FALSE)
    sam.gen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
    snps <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
    fam.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "FamID"))
    
    ## Produce gvar file in parallel processing Identify the SO
    os <- rappdirs:::get_os()
    if(os %in% c("unix", "mac")) param <- BiocParallel::MulticoreParam(workers = n.cor)
    else param <- BiocParallel::SnowParam(workers = n.cor, type = "SOCK")    
        
    if(run.lrr) .fun <- .prodGvarLRR
    else .fun <- .prodGvar

    Gens <- suppressMessages(BiocParallel::bplapply(1:length(sam.gen), .fun, 
                BPPARAM = param, genofile = genofile, sam.gen = sam.gen, 
                snps = snps, fam.id = fam.id, snp.matrix = snp.matrix))
    gentype.all <- data.table::rbindlist(Gens)
    gentype.all1 <- as.data.frame(gentype.all)
    
    gvar.file <- file.path(all.paths["PLINK"], "mydata.gvar")
    write.table(gentype.all1, gvar.file, sep = "\t", row.names = FALSE, quote = FALSE)
    SNPRelate::snpgdsClose(genofile)
}


# HELPER - Produce the PLINK fam file in disk
.prodPLINKfam <- function(all.paths) 
{
    cnv.gds <- file.path(all.paths["Inputs"], "CNV.gds")
    genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork = TRUE, readonly = FALSE)
    SamGen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
    SNPs <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))
    FAMID <- as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "FamID")))
    Sex <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "Sex"))
    Phen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype"))
    
    famdf <- data.frame(FAMID, SamGen, SamGen, "NC", Sex, as.numeric(Phen[, 2]))
    colnames(famdf) <- c("FAID", "IID", "PID", "MID", "Sex", "Phenotype")
    
    fam.file <- file.path(all.paths["PLINK"], "mydata.fam")
    write.table(famdf, fam.file, sep = "\t", row.names = FALSE, quote = FALSE)
    SNPRelate::snpgdsClose(genofile)
}

# HELPER - Run PLINK
.runPLINK <- function(all.paths) 
{
    os <- rappdirs:::get_os()
    if (os %in% c("win", "unix")) {
        plinkPath <- file.path(all.paths["PLINK"], "plink")
        plinkPath <- gsub("\\\\", "/", plinkPath)
        system(paste(plinkPath, "--gfile", file.path(all.paths["PLINK"], "mydata"), paste("--out", 
            file.path(all.paths["PLINK"], "plink")), "--noweb"), wait = TRUE, intern = TRUE)
    }
    else if (os == "mac") {
        plink.path <- file.path(all.paths["PLINK"], "plink")
        mydata <- file.path(all.paths["PLINK"], "mydata")
        mydata <- paste0("'", mydata, "'")
        plink <- file.path(all.paths["PLINK"], "plink")
        plink <- paste0("'", plink, "'")
        args <- c("--gfile", mydata, "--out", plink, "--noweb", "--allow-no-sex")
        system2(plink.path, args = args)
    }
    
    cnv.gds <- file.path(all.paths["Inputs"], "CNV.gds")
    genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork = TRUE, readonly = FALSE)
    sumphen <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype"))
    sumphen <- paste0("plink.gvar.summary.", names(sumphen)[[2]])
   
    from <- file.path(all.paths["PLINK"], "plink.gvar.summary")
    to <- file.path(all.paths["PLINK"], sumphen)
    file.rename(from, to)
    SNPRelate::snpgdsClose(genofile)
}

# HELPER - Produce CNV segments
.prodCNVseg <- function(all.paths, probes.cnv.gr, min.sim, both.up.down) 
{
    cnv.gds <- file.path(all.paths["Inputs"], "CNV.gds")
    genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork = TRUE, readonly = FALSE)
   
    snp.chr <-  gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")) 
    chrx <- data.frame(V1=snp.chr, V2=seq_along(snp.chr))
    g1 <- g2 <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "CNVgenotype"))
    
    chrx$V1 <- factor(chrx$V1, levels = unique(chrx$V1))
    chrx.split <- split(chrx, chrx$V1)
    
    # Find similarity between subsequent probes upstream 
    # g1 <- apply(g1,2,rev) 
    # g2 <- apply(g2,2,rev)
    SimiLar.all <- NULL
    for (lo in seq_along(chrx.split)) {
        snpindex <- as.numeric(as.character(chrx.split[[lo]]$V2))
        
        if (length(snpindex) == 1) {
            SimiLar.all[[lo]] <- "UNIQUE"
            next
        }
        
        g1x <- g1[snpindex, ]
        g1x <- apply(g1x, 2, rev)
        
        if (length(snpindex) == 2) {
            g1x <- rbind(g1x, rep(2, ncol(g1x)))  # If only two SNPs input a 2n row
        }
        
        g1xM <- g1x[-(nrow(g1x)), ]
        
        
        g1CN <- apply(g1xM, 1, function(x) sum(table(x)[names(table(x)) != 2]))
        g1x <- as.data.frame(t(sapply(seq_len(nrow(g1x)), 
            function(i) replace(g1x[i, ], g1x[i, ] == 2, paste0(i, "R")))))
        
        g2x <- g2[snpindex, ]
        g2x <- apply(g2x, 2, rev)
        
        # If only two SNPs input a 2n row
        if (length(snpindex) == 2) g2x <- rbind(g2x, rep(2, ncol(g2x)))

        g2x <- as.data.frame(t(sapply(seq_len(nrow(g2x)), function(i) replace(g2x[i, ], 
            g2x[i, ] == 2, paste0(i, "R")))))
        g1x <- g1x[-nrow(g1x), ]
        g2x <- g2x[-1, ]
        ftma <- g1x == g2x
        SimiLar <- rowSums(ftma, na.rm=TRUE)
        
        simX <- SimiLar / g1CN
        if (length(snpindex) == 2) simX <- simX[-2]
        SimiLar.all[[lo]] <- simX  ### CNV genotype similarity between subsequent SNPs
    }
    SimiLar.all.Up <- SimiLar.all
    
    ### Find similarity between subsequent probes Downstream
    SimiLar.all <- NULL
    for (lo in seq_along(chrx.split)) 
    {
        snpindex <- as.numeric(as.character(chrx.split[[lo]]$V2))
        
        if (length(snpindex) == 1) {
            SimiLar.all[[lo]] <- "UNIQUE"
            next
        }
        
        g1x <- g1[snpindex, ]
        # If only two SNPs input a 2n row
        if (length(snpindex) == 2) g1x <- rbind(g1x, rep(2, ncol(g1x)))  
        
        g1xM <- g1x[-(nrow(g1x)), ]
        g1CN <- apply(g1xM, 1, function(x) sum(table(x)[names(table(x)) != 2]))
        g1x <- as.data.frame(t(sapply(seq_len(nrow(g1x)), function(i) replace(g1x[i, ], 
            g1x[i, ] == 2, paste0(i, "R")))))
        
        g2x <- g2[snpindex, ]
        # If only two SNPs input a 2n row
        if (length(snpindex) == 2) g2x <- rbind(g2x, rep(2, ncol(g2x)))  
        
        g2x <- as.data.frame(t(sapply(seq_len(nrow(g2x)), function(i) replace(g2x[i, ], 
            g2x[i, ] == 2, paste0(i, "R")))))
        g1x <- g1x[-(nrow(g1x)), ]
        g2x <- g2x[-1, ]
        ftma <- g1x == g2x
        SimiLar <- rowSums(ftma, na.rm=TRUE)
        
        simX <- SimiLar / g1CN
        if (length(snpindex) == 2) simX <- simX[-2]

        # CNV genotype similarity between subsequent SNPs
        SimiLar.all[[lo]] <- simX  
    }
    SimiLar.all.Down <- SimiLar.all
    
    if (both.up.down) {
        for (lo in seq_along(SimiLar.all)) {
            SimiLar.all[[lo]] <- pmin(SimiLar.all.Up[[lo]], SimiLar.all.Down[[lo]])
        }
    }
    
    ### Produce segments
    all.segs <- NULL
    for (ch in seq_along(SimiLar.all)) {
        SimiLarX <- SimiLar.all[[ch]]
        if (all(SimiLarX == "UNIQUE")) {
            all.segs[[ch]] <- "start"
            next
        }
        nextP <- NULL
        all.segsch <- NULL
        for (se in seq_along(SimiLarX)) 
        {
            nextP[[se]] <- SimiLarX[[se]] >= min.sim
            
            if (se == 1) all.segsch[[se]] <- "start"
            else all.segsch[[se]] <- ifelse(nextP[[se - 1]], "cont", "start")
                
            if (se == length(SimiLarX)) 
                all.segsch[[se + 1]] <- ifelse(nextP[[se]], "cont", "start")
        }
        all.segs[[ch]] <- all.segsch
    }
    
    probes.cnv.gr$starts.seg <- unlist(all.segs)
    
    pstarts <- probes.cnv.gr[probes.cnv.gr$starts.seg == "start"]
    pstarts.split = split(pstarts, GenomeInfoDb::seqnames(pstarts))
    
    all.segs.gr <- GenomicRanges::GRangesList()
    for (chx in seq_along(pstarts.split)) {
        pstartsx <- pstarts.split[[chx]]
        if (!length(pstartsx)) {
            all.segs.gr[[chx]] <- GenomicRanges::reduce(pstartsx)  #Empty if no segments
            next
        }
        all.segs.grX <- GenomicRanges::GRangesList()
        for (st in seq_along(pstartsx)) {
            prx <- pstartsx[st]
            if (st < length(pstartsx)) {
                LimPrx <- GenomicRanges::start(pstartsx[st + 1]) - 1
            }
            if (st == length(pstartsx)) {
                sel.probes <- probes.cnv.gr[GenomeInfoDb::seqnames(probes.cnv.gr) == 
                  GenomeInfoDb::seqnames(prx)]
                LimPrx <- max(GenomicRanges::start(sel.probes))
            }
            seqLar <- GenomicRanges::GRanges(as.character(GenomeInfoDb::seqnames(prx)), 
                IRanges::IRanges(GenomicRanges::start(prx), LimPrx))
            GenomeInfoDb::seqlevels(seqLar) <- GenomeInfoDb::seqlevels(probes.cnv.gr)
            seqPrbs <- suppressWarnings(IRanges::subsetByOverlaps(probes.cnv.gr, 
                seqLar))
            seqx <- GenomicRanges::GRanges(as.character(GenomeInfoDb::seqnames(prx)), 
                IRanges::IRanges(GenomicRanges::start(prx), GenomicRanges::start(seqPrbs[length(seqPrbs)])))
            all.segs.grX[[st]] <- seqx
        }
        suppressWarnings(all.segs.gr[[chx]] <- unlist(all.segs.grX))
    }
    
    SNPRelate::snpgdsClose(genofile)
    
    return(all.segs.gr)
    
}


# HELPER - Associate probes with CNV segments and draw p-values
.assoPrCNV <- function(all.paths, all.segs.gr, phenotypesSamX, method.m.test, 
    probes.cnv.gr, assign.probe = "min.pvalue", correct.inflation) {
    
    segs.pvalue.gr <- unlist(all.segs.gr)
    segs.pvalue.gr$SegName <- seq_along(segs.pvalue.gr)
    
    cnv.gds <- file.path(all.paths["Inputs"], "CNV.gds")
    genofile <- SNPRelate::snpgdsOpen(cnv.gds, allow.fork = TRUE, readonly = FALSE)
    
    gvar.name <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "phenotype"))
    gvar.name <- names(gvar.name)[[2]]
    sum.file <- paste0("plink.gvar.summary.", gvar.name)
    results <- read.table(file.path(all.paths["PLINK"], sum.file), header = TRUE)
    mydata.map <- read.table(file.path(all.paths["PLINK"], "mydata.map"), header = TRUE)
    resultsp <- merge(results, mydata.map, by.x = "NAME", by.y = "NAME", all.x = TRUE, 
        all.y = FALSE)
    
    # Choose the method
    # TODO: change if SNPMatrix implemented
    # allow to get other p-values
    resultsp <- resultsp[grep("P\\(CNP\\)", resultsp$FIELD), ]
    resultsp$VALUE <- as.numeric(as.character(resultsp$VALUE))
    resultsp$VALUE <- round(resultsp$VALUE, digits = 5)
    
    #### Correct for genomic inflation
    all.samples <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
    if (correct.inflation) {
        #### Calculating chi-square distribution based on original P-values
        chisq.n <- qchisq(resultsp$VALUE, 1, lower.tail = FALSE)
        ### Estimating genomic inflation factor
        lambda <- estlambda(resultsp$VALUE, filter = FALSE)$estimate
        ### correcting the qui-square distribution with lambda
        chisq_corr = chisq.n/lambda
        #### re-calculating corrected P values
        resultsp$VALUE = pchisq(chisq_corr, 1, lower.tail = FALSE)
    }
    
    resultsp <- resultsp[order(resultsp$VALUE), ]
    resultspGr <- GenomicRanges::makeGRangesFromDataFrame(resultsp, seqnames.field = "Chr", 
        start.field = "Position", end.field = "Position", keep.extra.columns = TRUE)
    
    ######## Assign the lowest p-value to the CNV segment
    
    values.all <- vector("list", length(segs.pvalue.gr))
    names.all <- vector("list", length(segs.pvalue.gr))
    freqs.all <- vector("list", length(segs.pvalue.gr))
    
    if (assign.probe == "min.pvalue") {
        for (se in seq_along(segs.pvalue.gr)) {
            Prbx <- suppressWarnings(IRanges::subsetByOverlaps(resultspGr, segs.pvalue.gr[se]))
            if (length(Prbx)) {
                values.all[[se]] <- min(Prbx$VALUE)
                names.all[[se]] <- as.character(Prbx[Prbx$VALUE %in% values.all[[se]]]$NAME[[1]])
                probe.sel <- Prbx[order(Prbx$VALUE)][1]
                probe.sel <- suppressWarnings(IRanges::subsetByOverlaps(probes.cnv.gr, probe.sel))
                freqs.all[[se]] <- as.character(probe.sel$freq[1])
            } else {
                probe.sel <- suppressWarnings(IRanges::subsetByOverlaps(probes.cnv.gr, segs.pvalue.gr[se])[1])
                values.all[[se]] <- 1
                names.all[[se]] <- as.character(probe.sel$Name)
                freqs.all[[se]] <- as.character(probe.sel$freq[1])
            }
        }
    }
    
    if (assign.probe == "high.freq") {
        for (se in seq_along(segs.pvalue.gr)) {
            GenomeInfoDb::seqlevels(segs.pvalue.gr) <- GenomeInfoDb::seqlevels(probes.cnv.gr)
            probe.sel <- suppressWarnings(IRanges::subsetByOverlaps(probes.cnv.gr, 
                segs.pvalue.gr[se]))
            probe.sel <- probe.sel[rev(order(probe.sel$freq))][1]
            ## Take the first probe if the position is duplicated
            GenomeInfoDb::seqlevels(probe.sel) <- GenomeInfoDb::seqlevels(resultspGr)
            resultX <- suppressWarnings(IRanges::subsetByOverlaps(resultspGr, probe.sel)[1])
            values.all[[se]] <- resultX$VALUE
            names.all[[se]] <- as.character(probe.sel$Name)
            freqs.all[[se]] <- as.character(probe.sel$freq[1])
        }
    }
    
    if (assign.probe == "median") {
        for (se in seq_along(segs.pvalue.gr)) {
            GenomeInfoDb::seqlevels(segs.pvalue.gr) <- GenomeInfoDb::seqlevels(probes.cnv.gr)
            probe.sel <- suppressWarnings(IRanges::subsetByOverlaps(probes.cnv.gr, 
                segs.pvalue.gr[se]))
            absdiff <- abs(probe.sel$freq - median(probe.sel$freq))
            probe.sel <- probe.sel[which.min(absdiff)]
            ## Take the first probe if the position is duplicated
            GenomeInfoDb::seqlevels(probe.sel) <- GenomeInfoDb::seqlevels(resultspGr)
            resultX <- suppressWarnings(IRanges::subsetByOverlaps(resultspGr, probe.sel)[1])
            values.all[[se]] <- resultX$VALUE
            names.all[[se]] <- as.character(probe.sel$Name)
            freqs.all[[se]] <- as.character(probe.sel$freq[1])
        }
    }
    
    if (assign.probe == "mean") {
        for (se in seq_along(segs.pvalue.gr)) {
          GenomeInfoDb::seqlevels(segs.pvalue.gr) <- GenomeInfoDb::seqlevels(probes.cnv.gr)  
            probe.sel <- suppressWarnings(IRanges::subsetByOverlaps(probes.cnv.gr, 
                segs.pvalue.gr[se]))
            absdiff <- abs(probe.sel$freq - mean(probe.sel$freq))
            probe.sel <- probe.sel[which.min(absdiff)]
            ## Take the first probe if the position is duplicated
            GenomeInfoDb::seqlevels(probe.sel) <- GenomeInfoDb::seqlevels(resultspGr)  
            resultX <- suppressWarnings(IRanges::subsetByOverlaps(resultspGr, probe.sel)[1])
            values.all[[se]] <- resultX$VALUE
            names.all[[se]] <- as.character(probe.sel$Name)
            freqs.all[[se]] <- as.character(probe.sel$freq[1])
        }
    }
    
    segs.pvalue.gr$MinPvalue <- unlist(values.all)
    segs.pvalue.gr$NameProbe <- unlist(names.all)
    segs.pvalue.gr$Frequency <- unlist(freqs.all)
    
    segs.pvalue.gr$MinPvalueAdjusted <- stats::p.adjust(segs.pvalue.gr$MinPvalue, 
        method = method.m.test, n = length(segs.pvalue.gr))
    segs.pvalue.gr$Phenotype <- names(phenotypesSamX)[[2]]
    SNPRelate::snpgdsClose(genofile)
    
    return(segs.pvalue.gr)
    
}

# HELPER - Function to apply during LRR/BAF import @param lo loop number, based
# on the number of SNPs per turn @param list.filesLo Data-frame with two columns
# where the (i) is the path + file name with signals and (ii) is the
# correspondent name of the sample in the gds file @param genofile loaded gds
# file @param all.samples All samples in the gds file @param nLRR Connection to
# write LRR values @param nBAF Connection to write BAF values @param
# snps.included All SNP probe names to be included @param verbose Print the
# samples while importing
.freadImport <- function(lo, list.filesLo, genofile, all.samples, nLRR = NULL, nBAF = NULL, 
    snps.included, verbose) {
    
    if (verbose) message(paste0("sample ", lo, " of ", length(all.samples)))
    file.x <- c(as.character(list.filesLo[lo, 1]), as.character(list.filesLo[lo,2]))
    
    ### Import the sample file with fread
    sig.x <- data.table::fread(file.x[1], skip = 1L, header = FALSE)
    sig.x <- as.data.frame(sig.x)
    colnames(sig.x) <- c("name", "chr", "position", "lrr", "baf")
    
    ### Order the signal file as in the gds
    ind <- as.vector(sig.x$name) %in% snps.included
    sig.x <- sig.x[ind, ]
    
    ### Include missing SNPs
    missing.snps <- snps.included[!(snps.included %in% sig.x$name)]
    if (length(missing.snps) > 0) {
        missing.snps <- data.frame(missing.snps, chr = "NoN", position = "NoN", lrr = NA, 
            baf = NA)
        colnames(missing.snps)[1] <- "name"
        sig.x <- rbind(sig.x, missing.snps)
    }
    
    sig.x <- sig.x[order(match(sig.x[, 1], snps.included)), ]
    
    ### Extract the LRR
    LRR.x <- sig.x$lrr
    
    ### Extract the BAF
    BAF.x <- sig.x$baf
    
    ### Find the number of the sample
    sam.num <- which(all.samples == file.x[2])
    
    if (!is.null(nLRR)) {
        ### Write LRR value for sample x
        gdsfmt::write.gdsn(nLRR, LRR.x, start = c(1, sam.num), count = c(length(snps.included), 
            1))
    } else if (!is.null(nBAF)) {
        ### Write BAF value for sample x
        gdsfmt::write.gdsn(nBAF, BAF.x, start = c(1, sam.num), count = c(length(snps.included), 
            1))
    }
    
}

# HELPER - Estimate lambda - From GeneAbel package that was temporarily removed from CRAN
estlambda <- function(data, plot = FALSE, proportion = 1, method = "regression", 
          filter = TRUE, df = 1, ...) 
{
  data <- data[which(!is.na(data))]
  if (proportion > 1 || proportion <= 0) 
    stop("proportion argument should be greater then zero and less than or equal to one")
  ntp <- round(proportion * length(data))
  if (ntp < 1) 
    stop("no valid measurements")
  if (ntp == 1) {
    warning(paste("One measurement, lambda = 1 returned"))
    return(list(estimate = 1, se = 999.99))
  }
  if (ntp < 10) 
    warning(paste("number of points is too small:", ntp))
  if (min(data) < 0) 
    stop("data argument has values <0")
  if (max(data) <= 1) {
    data <- qchisq(data, 1, lower.tail = FALSE)
  }
  if (filter) {
    data[which(abs(data) < 1e-08)] <- NA
  }
  data <- sort(data)
  ppoi <- stats::ppoints(data)
  ppoi <- sort(qchisq(ppoi, df = df, lower.tail = FALSE))
  data <- data[seq_len(ntp)]
  ppoi <- ppoi[seq_len(ntp)]
  out <- list()
  if (method == "regression") {
    s <- summary(stats::lm(data ~ 0 + ppoi))$coeff
    out$estimate <- s[1, 1]
    out$se <- s[1, 2]
  }
  else if (method == "median") {
    out$estimate <- median(data, na.rm = TRUE) / qchisq(0.5, df)
    out$se <- NA
  }
  else {
    stop("'method' should be either 'regression' or 'median'!")
  }
  if (plot) {
    lim <- c(0, max(data, ppoi, na.rm = TRUE))
    oldmargins <- graphics::par()$mar
    graphics::par(mar = oldmargins + 0.2)
    plot(ppoi, data, xlab = expression("Expected " ~ chi^2), 
         ylab = expression("Observed " ~ chi^2), ...)
    graphics::abline(a = 0, b = 1)
    graphics::abline(a = 0, b = out$estimate, col = "red")
    graphics::par(mar = oldmargins)
  }
  out
}
