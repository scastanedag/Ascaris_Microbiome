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

# Change this path to the one that corresponds to where the reads are located

miseq_path <- "/Samples"
list.files(miseq_path)

Fs_path <- sort(list.files(miseq_path, pattern="_1.fastq.gz", full.names =
                             TRUE))
Rs_path <- sort(list.files(miseq_path, pattern="_2.fastq.gz", full.names =
                             TRUE))

# Create directories for processed forward and reverse reads
Fs_path_filtered <- file.path(miseq_path, "filtered_Fs")
Rs_path_filtered <- file.path(miseq_path, "filtered_Rs")

# Cargar el paquete stringr
library(stringr)

sample_names <- str_replace(string = basename(Fs_path),
                                    pattern = "_1\\.fastq.gz",
                                    replacement = "")

#quality profile

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)

# Create a directory to put the "clean" reads in
#filt_path <- file.path(miseq_path, "filtered") 
#if(!file_test("-d", filt_path)) dir.create(filt_path)
#filtFs <- file.path(filt_path, paste0(sampleNames, "_1_filt.fastq.gz"))
#filtRs <- file.path(filt_path, paste0(sampleNames, "_2_filt.fastq.gz"))


# Finally we proceed with the quality control.
out <- filterAndTrim(fwd=Fs_path,
                     filt=Fs_path_filtered,
                     rev=Rs_path,
                     filt.rev=Rs_path_filtered,
                     truncLen=c(200,200), # forward and reverse read
                     maxEE=c(2,2),
                     truncQ=2,
                     maxN=0,
                     rm.phix=TRUE,
                     compress=TRUE,
                     verbose=TRUE,
                     multithread=TRUE)

# Get list of filtered sequences
Fs_filt <- list.files(Fs_path_filtered, full.names = TRUE, pattern = "fastq.gz")
Rs_filt <- list.files(Rs_path_filtered, full.names = TRUE, pattern = "fastq.gz")

# Create names
names(Fs_filt) <- sample_names
names(Rs_filt) <- sample_names

# If you are using Windows, set multithread=FALSE
derepFs <- derepFastq(Fs_path_filtered, verbose=TRUE)
derepRs <- derepFastq(Rs_path_filtered, verbose=TRUE)

# We add the names of the samples to the de-replicated object
names(derepFs) <- sample_names
names(derepRs) <- sample_names


errF <- learnErrors(Fs_path_filtered, multithread=TRUE)
errR <- learnErrors(Rs_path_filtered, multithread=TRUE)

# We plot the errors for each pair
plotErrors(errF)
plotErrors(errR)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)

# Merged
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# We generate a sequence table
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))

# We remove the chemical sequences
seqtabNoC <- removeBimeraDenovo(seqtabAll)

num_chim_removed <- 1 - (sum(seqtabNoC)/sum(seqtabAll))

num_chim_removed

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtabNoC))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
write.csv(track, "preprocesamiento.csv")
#kableExtra::kable(track)


# Taxonomic Assignment
fastaRef <- "/silva_nr99_v138.1_train_set.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)

# In the case of wanting to add the taxonomic range of species, 
#we simply use an extra database, which contains this information. 
#taxTabExtra <- addSpecies(taxTab, "silva_species_assignment_v138.1.fa.gz", verbose=TRUE)

unname(head(taxTab)) -> tabla
colnames(tabla) <- c("Kingdom", "Phylum", "Order", "Class", "Family", "Genus")

write.csv(taxTab, "taxTab2.csv")
#write.csv(taxTabExtra, "taxTabExtra.csv")
write.csv(seqtabNoC, "seqtabNoC2.csv")
