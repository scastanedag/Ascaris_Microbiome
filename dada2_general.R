## DADA2 Analysis Pipeline
## This code represents an analysis pipeline for processing and analyzing DNA sequencing data using the DADA2 package in R. 
## The pipeline includes steps for quality filtering, denoising, merging paired-end reads, chimera removal, and taxonomic assignment.
## Load required library
old <- Sys.time() # get start time

set.seed(531) 
library(dada2)
packageVersion('dada2')

# Defined Packages
## CRAN
cran_packages <- c("bookdown", "knitr", "tidyverse", "plyr", "grid", "gridExtra", "kableExtra", "xtable", "ggpubr")
## ioconductor
bioc_packages <- c("phyloseq", "dada2", "DECIPHER", "phangorn", "ggpubr", "BiocManager","DESeq2", "microbiome", "philr")
## GitHub
git_source <- c("twbattaglia/btools", "gmteunisse/Fantaxtic", "MadsAlbertsen/ampvis2", "opisthokonta/tsnemicrobiota")
# name 
git_packages <- c("btools", "fantaxtic", "ampvis2", "tsnemicrobiota") # package name

#Second, install the packages defined above using the function corresponding to each repository.
# Install CRAN packages
.inst <- cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(cran_packages[!.inst])
}
# Install packages from BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
.inst <- bioc_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(bioc_packages[!.inst])
}
# nstall packages from  GitHub
.inst <- git_source %in% installed.packages()
if(any(!.inst)) {
  devtools::install_github(git_source[!.inst])
}

# Load packages
sapply(c(cran_packages, bioc_packages, git_packages), require, character.only = TRUE)

# Load packages
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(kableExtr)
library(xtable)
library(ggpubr)
library(phyloseq)
library(dada2)
library(DECIPHER)
library(phangorn)
library(ggpubr)
library(BiocManager)
library(DESeq2)
library(microbiome)
library(philr)
library(btools)
library(fantaxtic)
library(ampvis2)
library(tsnemicrobiota)

## Specify the path for input FASTQ files (replace with the appropriate path after downloading the files from ENA/EGA)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dirname(rstudioapi::getActiveDocumentContext()$path)

# replace by directory with fastq files
path <- paste(dirname(rstudioapi::getActiveDocumentContext()$path),'/fastq', sep = '')
fns <- list.files(path)
fastqs <- fns[grepl('.fastq.gz$', fns)]
fastqs <- sort(fastqs)

## Separate forward and reverse read files and extract sample names
fnFs <- fastqs[grepl('_1.fastq.gz', fastqs)]
fnRs <- fastqs[grepl('_2.fastq.gz', fastqs)]
sample.names <- sapply(strsplit(fnFs, '_1.fastq.gz'), `[`, 1)

## Fully specify file paths for forward and reverse reads
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

## Set up filtered file paths and perform filtering and trimming
filt_path <- file.path(path, 'filtered')
filtFs <- file.path(filt_path, paste0(sample.names, '_F_filt.fastq.gz'))
filtRs <- file.path(filt_path, paste0(sample.names, '_R_filt.fastq.gz'))

## Plot overall quality
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

## MiSeq (LCPM cohort) and HiSeq (FGFP cohort) tested for my primer constructs
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(260, 230), trimLeft = c(5, 5), 
                     maxN = 0, maxEE = c(2, 2), truncQ = 11, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

head(out)
## Create a data frame with filtering results
DFout <- data.frame(out)
## Select samples with at least 50 reads
sample.names0 <- sapply(strsplit(rownames(subset(DFout, reads.out > 50)), '_1.fastq.gz'), `[`, 1)
filtFs0 <- file.path(filt_path, paste0(sample.names0, '_F_filt.fastq.gz'))
filtRs0 <- file.path(filt_path, paste0(sample.names0, '_R_filt.fastq.gz'))

## Plot overall quality
plotQualityProfile(filtFs0[1:2])
plotQualityProfile(filtRs0[1:2])

## Learn forward and reverse error rates 
errF <- learnErrors(filtFs0, nbases = 1e+05, multithread = TRUE) # minimum nbases = 1e+08,
errR <- learnErrors(filtRs0, nbases = 1e+05, multithread = TRUE) # minimum nbases = 1e+08,

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
## Dereplicate reads
derepRs <- derepFastq(filtRs0, verbose = TRUE)
derepFs <- derepFastq(filtFs0, verbose = TRUE)

## Assign names to dereplicated objects
names(derepFs) <- sample.names0
names(derepRs) <- sample.names0

## Perform denoising and merge paired-end reads
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

## Merge the denoised reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

## Create a sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
seqtabHS <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

## Save the sequence table
saveRDS(seqtabHS, "dd2.seqtab.rds")

###### 

seqtab_all <- seqtabHS## Remove chimeras
seqtab_all_no_chimeras <- removeBimeraDenovo(seqtab_all, method = "consensus", multithread = TRUE)
saveRDS(seqtab_all_no_chimeras, "dd2.seqtab_all_no_chimeras.rds")

## Assign taxonomy

## Assign taxonomy using the RDP training set 18 (replace with your reference file)

# Load RDP training set
rdp_file <- paste(dirname(rstudioapi::getActiveDocumentContext()$path), '/database/rdp_train_set_18.fa.gz', sep = '')
ASV_seq_rdp_set18 <- assignTaxonomy(seqtab_all_no_chimeras, rdp_file, multithread = TRUE)
saveRDS(ASV_seq_rdp_set18, "ASV_seq_rdp_set18.rds")

dim(seqtab)
## [1]  20 293
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
sum(seqtabHS)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtabHS))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
filereport_read_run_PRJNA525566_tsv <- data.frame(read.delim("filereport_read_run_PRJNA525566_tsv.txt"))

rownames(filereport_read_run_PRJNA525566_tsv)
filereport_read_run_PRJNA525566_tsv$run_accession

rownames(filereport_read_run_PRJNA525566_tsv) <- 
filereport_read_run_PRJNA525566_tsv$run_accession

# Phyloseq obj

library(phyloseq)

Phy_obj_RDP <- merge_phyloseq(sample_data(filereport_read_run_PRJNA525566_tsv), 
               otu_table(seqtab_all_no_chimeras, taxa_are_rows = FALSE),
               tax_table(ASV_seq_rdp_set18))
Phy_obj_RDP
save(Phy_obj_RDP, file="Phy_obj_RDP")
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) 
