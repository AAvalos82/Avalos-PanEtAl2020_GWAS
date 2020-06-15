#The purpose of this script is to conduct a genome wide association analysis using per-colony phenotype and minor allele.

  #input:
    #worker_mvcall.phasedimputed.snp_subset.gds
    #all210_idx.RDS
    #kin_snpset.RDS

  #output:
    #SNP_minor_allele_index.RDS
    #colony_covariate_matrix.RDS
    #colony_whole_glm.RDS

#------------------------------------------------------------------------------- 
# Load analysis options and libraries
#------------------------------------------------------------------------------- 

#libraries
pkgs <- c("tidyverse", "magrittr", "SeqArray", "RColorBrewer", "data.table", 
          "SNPRelate", "igraph", "Rtsne", "plotly", "MASS", "car", "viridis", 
          "mclust", "BiocParallel", "GenomicRanges", "GENESIS", "GWASTools",
          "compiler")
invisible(lapply(pkgs, library, character.only = T))
rm(pkgs); gc()

#options
options(stringsAsFactors = F, scipen = 9999)
set.seed(12345)

#------------------------------------------------------------------------------- 
# Conduct formal GWAS analysis
#------------------------------------------------------------------------------- 

#read in the data bases
snp.file <- snpgdsOpen(
  "./gwas_individual/worker_mvcall.phasedimputed.snp_subset.gds"
)
meta     <- readRDS("./all210_idx.RDS")
snp.set  <- readRDS("./kinship/kin_snpset.RDS")

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
dir.create("./gwas_colony", showWarnings = F)

#create the phenotype indices
col.agg <- setNames(object = meta.subset$MDS1_Aggression %>% unique,
                    nm     = meta.subset$colony %>% unique)

#construct the working GRanges
snp.gr          <- GRanges(seqnames = snp.file            %>% 
                             index.gdsn("snp.chromosome") %>%
                             read.gdsn,
                           ranges   = IRanges(start = snp.file             %>% 
                                                index.gdsn("snp.position") %>%
                                                read.gdsn, 
                                              width = 1))
snp.gr$ID       <- snp.file %>% 
  index.gdsn("snp.id")      %>% 
  read.gdsn
snp.gr$genotype <- snp.file %>% 
  index.gdsn("genotype")    %>% 
  read.gdsn                 %>% 
  t
colnames(snp.gr$genotype) <- snp.file %>% 
  index.gdsn("sample.id")             %>% 
  read.gdsn                           %>% 
  meta.colony[.]

#calculate allele frequency per colony
#determine whether the minor allele frequency is same as reference
snp.maf <- snpgdsSNPRateFreq(snp.file) %>% with(AlleleFreq <= MinorFreq) 
#create a housing object
snp.gr.maf    <- snp.gr[, -(1:2)]
snp.gr.maf$GT <- snp.file %>% index.gdsn("snp.allele") %>% read.gdsn()
snp.gr.maf$MA <- snp.maf
#save this index for reference
saveRDS(snp.gr.maf, "./gwas_colony/SNP_minor_allele_index.RDS", compress = T)
#remove the extraneous index and clean-up
rm(snp.gr.maf); gc()
#calculate the frequency for the minor allele in the data set
snp.gr$cfreq <- meta.colony %>% 
  unique                    %>% 
  sort                      %>% 
  sapply(FUN = function(i, gr, maf) {
    tmp = gr$genotype[, colnames(gr$genotype) == i]
    out = ifelse(maf, rowSums(tmp) / (ncol(tmp) * 2), 
                 1 - (rowSums(tmp) / (ncol(tmp) * 2)))
    return(out)
  },
         gr  = snp.gr, 
         maf = snp.maf)

#construct the colony covariate matrix from the minor allele frequency
#create the index of target SNPs from the kinship analysis
tsnp.idx <- snp.gr$ID %in% snp.set
#transpose and center the matrix
col.covar <- snp.gr[tsnp.idx]$cfreq %>% 
  t                                 %>% 
  scale(center = T)                 %>% 
  #matrix multiply by its transposition and divide by the number of markers
  tcrossprod %>% `/` (length(snp.gr[tsnp.idx]))
#save the object
saveRDS(object   = col.covar, 
        file     = "./gwas_colony/colony_covariate_matrix.RDS",
        compress = T)
#decompose via PCA
col.covar.pca <- prcomp(x = col.covar, center = F, scale. = F)

#correlate colony phenotype with variant allele frequency
#construct the function and compile
glm_fit <- function(i, freq, pheno, grm) {
  #fit the test model
  null     = glm(pheno ~ 1 + grm, family = gaussian)
  tmp      = glm(pheno ~ 1 + freq[i, ] + grm, family = gaussian)
  #extract metrics
  out      = data.frame(
    beta   = summary(tmp)$coefficients[2, 1],
    pvl    = anova(null, tmp, test = "LRT")$`Pr(>Chi)`[2]
  )
  #produce output
  return(out)
}
glm_fit <- cmpfun(glm_fit)
#set up the workers
cl <- MulticoreParam(detectCores() - 2)
#if on Windows
#cl <- SnowParam(detectCores() - 2)
#conduct the analysis
col.glm <- bplapply(X       = 1:length(snp.gr), 
                    BPPARAM = cl, FUN = glm_fit,
                    freq    = snp.gr$cfreq, 
                    pheno   = col.agg[colnames(snp.gr$cfreq)],
                    grm     = col.covar.pca$x[, 1]) %>% 
                    {do.call(rbind, .)}
#add to the existing object
snp.gr$cglm <- col.glm

#save the object
saveRDS(object   = snp.gr, 
        file     = "./gwas_colony/colony_whole_glm.RDS", 
        compress = T)

