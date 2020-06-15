# Avalos-PanEtAl2020_GWAS
Repository containing raw code for kinship estimate and also individual-, and colony-level GWAS conducted in the Avalos &amp; Pan et al., 2020 study

# File information
01_kinship_final.R         
- The purpose of this script is to estimate the kinship and genetic relationship matrices.

02_GWAS_individual_final.R 
- The purpose of this script is to conduct a genome wide association analysis of individual phenotype.

03_GWAS_colony_final.R     
- The purpose of this script is to conduct a genome wide association analysis using per-colony phenotype and minor allele.

all210_idx.RDS             
- This is a file containing per-sample meta data for individuals used in this study.

snpset.RDS                 
- This is a file containing IDs of markers concordant between two independent variant-calling pipelines (see source publication for details).

# Note
In the code file: "worker_mvcall.phasedimputed.snp.gds" is referenced, this is a dataset file deposted in Dryad see:

Avalos, Arián (2020), Genomic regions influencing aggressive behavior in honey bees are defined by colony allele frequencies, Dryad, Dataset, https://doi.org/10.5061/dryad.q573n5tg8
