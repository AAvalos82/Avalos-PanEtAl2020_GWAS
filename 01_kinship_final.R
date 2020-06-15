#The purpose of this script is to estimate the kinship and genetic relationship matrices.

  #input:
    #worker_mvcall.phasedimputed.snp.gds
    #snpset.RDS
    #all210_idx.RDS

  #output:
    #kin_report.txt
    #kin_snpset.RDS
    #pca.RDS
    #kin_KING.RDS
    #pcair.RDS
    #pcrelate.RDS
    #grm.RDS

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
# Conduct kinship analysis and correction
#------------------------------------------------------------------------------- 

#read in the data bases
snp.file <- snpgdsOpen("./worker_mvcall.phasedimputed.snp.gds")
nsnp.set <- readRDS("./snpset.RDS")
meta     <- readRDS("./all210_idx.RDS")

#correct the numeric variables (for plotting purposes)
meta$colony %<>% as.numeric %>% {ifelse(. < 10, paste0("0", .), .)}
meta$ind    %<>% as.numeric %>% {ifelse(. < 10, paste0("0", .), .)}

#subset the metadata to samples of interest
meta.subset <- meta %>% subset(subset = .$include_alt)

#create required indices
meta.color <- meta.subset$colony %>% 
  unique                         %>%
  {setNames(object = brewer.pal(n = length(.), "Paired"), 
            nm     = .)}
meta.colony <- meta.subset %>% 
  {setNames(object = .$colony, 
            nm     = .$seqID)}
meta.sample <- meta.subset %>% 
  {setNames(object = paste(.$colony, .$ipheno, .$ind, sep = "."), 
            nm     = .$seqID)}

#create a housing directory
dir.create("./kinship", showWarnings = F)
dir.create("./kinship/figures", showWarnings = F)

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
#clean-up intermediates
rm(snp.poly); gc()

#isolate the concordant SNPs within the set
snp.poly.set %<>% .[. %in% nsnp.set]

#conduct LD prunning with MAF filter set at 20% for structure assessment
capture.output(
  snp.set <- snpgdsLDpruning(gdsobj         = snp.file, 
                             ld.threshold   = 0.3,
                             snp.id         = snp.poly.set,
                             sample.id      = meta.subset$seqID,
                             maf            = 0.20,
                             remove.monosnp = T,
                             autosome.only  = F, 
                             method         = "dprime",
                             slide.max.bp   = 100, 
                             num.thread     = (detectCores() / 2) - 1),
  file = "./kinship/kin_report.txt"
  )
#convert to vector
snp.set %<>% unlist(use.names = F)
#save the object
saveRDS(object   = snp.set, 
        file     = "./kinship/kin_snpset.RDS", 
        compress = T)

#construct the pca object
snp.pca <- snpgdsPCA(gdsobj        = snp.file, 
                     snp.id        = snp.set,
                     sample.id     = meta.subset$seqID, 
                     autosome.only = F, 
                     eigen.cnt     = 0,
                     eigen.method  = "DSPEV",
                     num.thread    = (detectCores() / 2) - 1)

#save the object
saveRDS(object   = snp.pca, 
        file     = "./kinship/pca.RDS", 
        compress = T)

#construct the kinship matrices
#using the KING method
snp.king <- snpgdsIBDKING(gdsobj        = snp.file, 
                          sample.id     = meta.subset$seqID, 
                          snp.id        = snp.set, 
                          autosome.only = F,
                          type          = "KING-robust",
                          num.thread    = (detectCores() / 2) - 1)
#save the object
saveRDS(object   = snp.king, 
        file     = "./kinship/kin_KING.RDS", 
        compress = T)

#------------------------------------------------------------------------------- 
# Conduct the GRM analysis
#------------------------------------------------------------------------------- 

#close the data base and re-open using the GENESIS accessor function
snpgdsClose(snp.file)
snp.file <- GdsGenotypeReader("./worker_mvcall.phasedimputed.snp.gds")

#extract the genotype data
snp.geno <- GenotypeData(snp.file)

#prepare the kinship matrices format to work with PC-AIR
KINGmat <- snp.king$kinship
rownames(KINGmat) <- colnames(KINGmat) <- snp.king$sample.id

#conduct the PC-AiR analysis to isolate unrelated individuals
snp.pcair <- pcair(gdsobj         = snp.geno, 
                   kinobj         = KINGmat, 
                   divobj         = KINGmat,
                   snp.include    = snp.set, 
                   sample.include = meta.subset$seqID, 
                   autosome.only  = F, 
                   eigen.cnt      = 0, 
                   eigen.method   = "DSPEV",
                   num.cores      = (detectCores() / 2) - 1)

#save the object
saveRDS(object   = snp.pcair, 
        file     = "./kinship/pcair.RDS", 
        compress = T)

#conduct the relatedness estimation adjustment for the PCA
geno.iter    <- GenotypeBlockIterator(genoData   = snp.geno, 
                                      snpBlock   = 50000,
                                      snpInclude = snp.set)
snp.pcrelate <- pcrelate(gdsobj         = geno.iter, 
                         pcs            = snp.pcair$vectors[, 1:7],
                         sample.include = meta.subset$seqID,
                         training.set   = snp.pcair$unrels)

#save the object
saveRDS(object   = snp.pcrelate, 
        file     = "./kinship/pcrelate.RDS", 
        compress = T)

#construct the genetic relatedness matrix from the output of the adjustment
snp.grm <- pcrelateToMatrix(pcrelobj       = snp.pcrelate,
                            sample.include = meta.subset$seqID)

#save the object
saveRDS(object   = snp.grm, 
        file     = "./kinship/grm.RDS", 
        compress = T)
