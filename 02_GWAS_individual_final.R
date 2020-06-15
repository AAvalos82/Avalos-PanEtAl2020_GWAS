#The purpose of this script is to conduct a genome wide association analysis 
#of the 2017 worker data set focusing on individual phenotype

  #input:
    #worker_mvcall.phasedimputed.snp.gds
    #kin_KING.RDS
    #pcair.RDS
    #pcrelate.RDS
    #grm.RDS
    #snpset.RDS
    #all210_idx.RDS
    #worker_mvcall.phasedimputed.snp.gds
    #worker_mvcall.phasedimputed.snp_subset.gds
    #individual_out.RDS
    #Assembly_Report.txt

  #output:
    #LD_prunning_results.txt
    #ldpreport.txt
    #snpset.RDS
    #worker_mvcall.phasedimputed.snp_subset.gds
    #individual_out.RDS
    
#------------------------------------------------------------------------------- 
# Load analysis options and libraries
#------------------------------------------------------------------------------- 

#libraries
pkgs <- c("UpSetR", "tidyverse", "magrittr", "RColorBrewer", "data.table", 
          "igraph", "Rtsne", "plotly", "MASS", "car", "viridis", "mclust", 
          "adegenet", "SeqArray", "SNPRelate", "igraph", "Rtsne", "plotly", 
          "MASS", "car", "viridis", "mclust", "BiocParallel", "GenomicRanges", 
          "GENESIS", "GWASTools")
invisible(lapply(pkgs, library, character.only = T))
rm(pkgs); gc()

#options
options(stringsAsFactors = F, scipen = 9999)
set.seed(12345)

#------------------------------------------------------------------------------- 
# Conduct formal GWAS analyses
#------------------------------------------------------------------------------- 

#read in the data bases
snp.file     <- snpgdsOpen("./worker_mvcall.phasedimputed.snp.gds")
snp.king     <- readRDS("./kinship/kin_KING.RDS")
snp.pcair    <- readRDS("./kinship/pcair.RDS")
snp.pcrelate <- readRDS("./kinship/pcrelate.RDS")
snp.grm      <- readRDS("./kinship/grm.RDS")
nsnp.set     <- readRDS("./snpset.RDS")
meta         <- readRDS("./all210_idx.RDS")

#correct the numeric variables (for plotting purposes)
meta$colony %<>% as.numeric %>% {ifelse(. < 10, paste0("0", .), .)}
meta$ind    %<>% as.numeric %>% {ifelse(. < 10, paste0("0", .), .)}

#subset the metadata to samples of interest
meta.subset <- meta %>% subset(subset = .$include_alt)

#create required indices
meta.color  <- meta.subset$colony %>% 
  unique %>%
  {setNames(object = brewer.pal(n = length(.), "Paired"), 
            nm     = .)}
meta.colony <- meta.subset %>% 
  {setNames(object = .$colony, 
            nm     = .$seqID)}
meta.sample <- meta.subset %>% 
  {setNames(object = paste(.$colony, .$ipheno, .$ind, sep = "."), 
            nm     = .$seqID)}

#create an analysis directory
dir.create("./gwas_individual", showWarnings = F)

#filter out homogeneous heterozygotes
#extract the genotype for the samples of interest
snp.poly <- snp.file     %>% 
  index.gdsn("genotype") %>% 
  read.gdsn              %>%
  subset(., subset = meta$include_alt)
#create a boolean vector identifying the polymorphic sites
snp.poly.set <- sapply(X   = 0:2, 
                       FUN = function(i) {
                         tmp = colSums(snp.poly == i)
                         out = tmp == dim(snp.poly)[[1]]
                       }) %>% 
  rowSums                 %>% 
  `==` (0)
#extract the snp IDs
snp.poly.set <- snp.file %>% 
  index.gdsn("snp.id")   %>% 
  read.gdsn              %>%
  `[` (snp.poly.set)
#isolate the concordant SNPs within the set
snp.poly.set %<>% .[. %in% nsnp.set]
#also apply a two-tail frequency cutoff 
#necessary because software equates non-reference with minor
snp.af <- snpgdsSNPRateFreq(snp.file) %>% 
  with(AlleleFreq <= 0.95 & AlleleFreq >= 0.05)
snp.af <- snp.file     %>% 
  index.gdsn("snp.id") %>% 
  read.gdsn            %>%
  `[` (snp.af)
#isolate the concordant SNPs within the set
snp.set <- snp.poly.set[snp.poly.set %in% snp.af]
#clean-up intermediates
rm(snp.poly, snp.poly.set, snp.af); gc()

#use the SNP set to generate the analysis file
snpgdsClose(snp.file)
snpgdsCreateGenoSet(
  src.fn    = "./worker_mvcall.phasedimputed.snp.gds",
  dest.fn   = "./gwas_individual/worker_mvcall.phasedimputed.snp_subset.gds",
  snp.id    = snp.set, 
  sample.id = names(meta.sample)
  )
snp.file <- GdsGenotypeReader(
  filename     = "./gwas_individual/worker_mvcall.phasedimputed.snp_subset.gds",
  autosomeCode = 1:16L
  )

#build metadata object
#prepare the data set annotation
scan.df <- data.frame(scanID = meta.subset$seqID,
                      snp.pcair$vectors[, 1:(length(snp.pcair$values) - 1)],
                      snp.pca$eigenvect[, 1:(ncol(snp.pca$eigenvect) - 1)],
                      sample_pheno = ifelse(
                        meta.subset$ipheno == "2", 1L, 0L
                      ),
                      colony = meta.subset$colony,
                      colony_pheno = meta.subset$MDS1_Aggression)
names(scan.df) <- c("scanID", paste0("pcair", 1:(length(snp.pcair$values) - 1)), 
                    paste0("pc", 1:(ncol(snp.pca$eigenvect) - 1)),
                    "sample_pheno", "colony", "colony_pheno")
scanAnnot <- ScanAnnotationDataFrame(scan.df)

#construct the annotated data base
geno <- GenotypeData(data      = snp.file, 
                     scanAnnot = scanAnnot)

#determine the number of PCs to include as covariates
#using elbow method with AIC
#set up the parallelization and function
aic.fun   <- function(i, pheno, kin, cov_pc) {
  #load libraries
  require(GENESIS)
  #establish the conditional
  if(i == 0){
    tmp = tryCatch(
      fitNullModel(x          = pheno, 
                   outcome    = "sample_pheno",
                   cov.mat    = list(GRM = kin),
                   family     = "binomial", 
                   verbose    = F), 
      error = function(e){})
    #extract the Akaike information criteria
    out = tmp$AIC
  } else {
    tmp = tryCatch(
      fitNullModel(x          = pheno, 
                   outcome    = "sample_pheno",
                   covars     = c(paste0(cov_pc, 1:i)),
                   cov.mat    = list(GRM = kin),
                   family     = "binomial", 
                   verbose    = F), 
      error = function(e){})
    #extract the Akaike information criteria
    out = tmp$AIC
    return(out)
  }
}
aic.fun   <- compiler::cmpfun(aic.fun)

#the grm is not organized in the same direction as the annotation data frame
#this is to fix this issue
mygrm <- as.matrix(snp.grm)[rownames(scan.df), rownames(scan.df)]

#conduct the ascertainment of optimal PC numbers
aic.pcair <- bplapply(X       = c(0, 1:ncol(snp.pcair$vectors)),
                      BPPARAM = SerialParam(),
                      FUN     = aic.fun,
                      pheno   = scan.df,
                      kin     = mygrm,
                      cov_pc  = "pcair") %>% 
  unlist(use.names = F)
#NOTE: This analysis showed that inclusion of these variables did not 
#improve model fit and as such individual analysis only includes the GRM
#as this was sufficient to account for kinship and structure variance.

#fit the null model
nullmod <- fitNullModel(x       = scan.df,
                        outcome = "sample_pheno", 
                        cov.mat = list(GRM = mygrm), 
                        family  = "binomial")

#conduct the correlation test for each SNP in the data set
#establish the iterator data object
geno.iter <- GenotypeBlockIterator(genoData   = geno, 
                                   snpBlock   = 50000,
                                   snpInclude = snp.set)
#run the analysis
assoc     <- assocTestSingle(gdsobj      = geno.iter, 
                             null.model  = nullmod, 
                             test        = "Score", verbose = T)

#save the output
saveRDS(object   = assoc, 
        file     = "./gwas_individual/individual_out.RDS", 
        compress = T)
