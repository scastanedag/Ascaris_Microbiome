library(dada2)
library(phyloseq)
library(ggpubr)
library(tibble)
library(dplyr)
library(microViz)
library(vegan)
library(pairwiseAdonis)
library(rstatix)
library(data.table)
library(tidyr)
library(microbiome)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(gridExtra)
library(scales)
library(vcd)
library(stats)
library(car)
library(Hmisc)
library(magrittr)
library(microbiomeMarker)
library(microbiomeutilities)
library(jeevanuDB)
library(ANCOMBC)
library(tidyverse)
library(nlme)
library(tidyverse)
library(compositions)
#source("programs/ancom.R")
library(ANCOMBC)
#install.packages("miaViz")
library(mia)
library(dplyr)
library(tidyr)
library(microbiomeMarker)
#library("treeio")
#install.packages("treeio")
#detach("ggtree", unload = TRUE)
#detach("treeio", unload = TRUE)
#detach("microbiomeMarker")
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
library(cowplot)
library(ggpubr)
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
#install.packages("tsnemicrobiota")
library("phyloseq")
library("ggplot2")
library("readxl")       # necessary to import the data from Excel file
library("dplyr")  # filter and reformat data frames
library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
#install_github("hallucigenia-sparsa/seqtime")
library(seqtime)

#install.packages("remotes")
#remotes::install_github("Russel88/MicEco")
library("MicEco")
#install.packages("vegan")
library("vegan")
#vignette("phyloseq-basics")
library(biome)
library(phyloseq)
library(metagenomeSeq)
library(iNEXT)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(data.table)
library(decontam); packageVersion("decontam")


# data nececsary to phyloseq object
otu_mat<- read_excel ("/seqtabNoC.xlsx")
tax_mat<- read_excel("/taxTab.xlsx")
samples_df <- read_excel("/data_sequencing_samples.xlsx")

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu")

tax_mat <- tax_mat %>%
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>%
  tibble::column_to_rownames("sample")

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)


# phyloseq object
ascaris <- phyloseq(OTU, TAX, samples)
ascaris

any(taxa_sums(ascaris) == 0)
sum(taxa_sums(ascaris) == 0)


PS <- ascaris
summarize_phyloseq(PS)

## HOW many ASVs for off-target eukaryotes and archaea
table(tax_table(PS)[, "Kingdom"], exclude = NULL) 

## HOW many reads for off-target eukaryotes and archaea
by((otu_table(PS)), tax_table(PS)[, "Kingdom"], sum) 

##Filtering
## Eliminate "empty" samples
PS <- prune_samples(sample_sums(PS)>0, PS)
summarize_phyloseq(PS)
##Sample filtering: Filtering samples with low counts
PS <- prune_samples(sample_sums(PS)>=2000, PS)
summarize_phyloseq(PS)

##Taxa filtering
PS<- subset_taxa(PS, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "NA"))
PS <- subset_taxa(PS, !Genus %in% c("Stenotrophomonas", "Stenotrophomonas 1", "Stenotrophomonas 2"))


##General check
plot_richness(PS, x= "Organ", color = "Organ" , measures = c("Observed","Simpson", "Shannon")) +
  geom_jitter(alpha= 0.005)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

## Remove low prevalent taxa
##Create a prevalence dataframe
Prevdf<- apply(X = otu_table(PS),
               MARGIN = 1,
               FUN = function(x){sum(x > 0)})

##Add taxonomy and total read counts to this data.frame
Prevdf<- data.frame(Prevalence = Prevdf,
                    TotalAbundance = taxa_sums(PS),
                    tax_table(PS))

plyr::ddply(Prevdf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

##Filter the Archaea phylum
phyla2Filter<- c("Euryarchaeota")
PS<- subset_taxa(PS, !Phylum %in% phyla2Filter)

Prevdf<- apply(X = otu_table(PS),
               MARGIN = 1,
               FUN = function(x){sum(x > 0)})
Prevdf<- data.frame(Prevalence = Prevdf,
                    TotalAbundance = taxa_sums(PS),
                    tax_table(PS))

ggplot(Prevdf, aes(TotalAbundance, Prevalence / nsamples(PS),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Log 10 Total Reads") + ylab("Prevalence [Prop. of Samples]") +
  theme_bw()+
  facet_wrap(~Phylum) + theme(legend.position="none")

#Remove low prevalent ASVs
##Remove ASVs that do not show appear more than 5 times in more than 10% of the samples
wh0 <- genefilter_sample(PS, filterfun_sample(function(x) x > 5), A=0.01*nsamples(PS))
PS<- prune_taxa(wh0, PS)


##Normalization of proportions
##Normalization transformation to an even sample size
PS.Norm<- transform_sample_counts(PS, function(x) 1E6 * x/sum(x))


##Normalization transformation to an even sample size
PS.Norm<- transform_sample_counts(PS, function(x) 1E6 * x/sum(x))

##Transform to even sampling depth
## Rarefy without replacement
vegan::rarecurve(t(otu_table(PS)), step=1000, cex=0.6)

## Rarefy without replacement to the min sequencing depth
PS.Rare<- rarefy_even_depth(PS, rngseed=2020, sample.size=min(sample_sums(PS)), replace=F)
PS.Rare #used for alpha diversity analysis

##Check how many samples ended after filtering
table(PS@sam_data$Organ, PS@sam_data$Groups_Inf_NI)

## Merge ASVs that have the same taxonomy at a certain taxonomic rank
PS.Fam<-  tax_glom(PS, "Family", NArm = F)

PS.Gen<-  tax_glom(PS, "Genus", NArm = T)

PS.Phy<-  tax_glom(PS, "Phylum", NArm = F)

#General Plot
plot_bar(PS.Gen, fill="Genus") +
  facet_wrap(~Groups_Inf_NI, scales= "free_x", nrow=1) +
  theme(legend.position = "none")

##Alpha diversity (rarefied)
plot_richness(PS.Rare, x= "Organ", color = "Groups_Inf_NI" , measures = c("Observed","Simpson", "Shannon")) +
  geom_jitter(alpha= 0.005)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

##Estimate alpha diversity for individual samples (no-rarefied)
alphadiv<- estimate_richness(PS)

##Estimate alpha diversity (rarefied)
alphadiv.rare<- estimate_richness(PS.Rare)
as.data.frame(PS@sam_data)->tmp
alphadiv.rare<-cbind(alphadiv.rare, tmp)

##Add sample data into a single data frame
as.data.frame(PS@sam_data)->tmp

alphadiv<-cbind(alphadiv, tmp)

table(alphadiv$Organ, alphadiv$Group_2, alphadiv$Groups_Inf_NI) 

saveRDS(PS, "/Users/scast/Desktop/SERGIO/Cursos_Extra/Barcelona_Training/Barcelona_Course/Ascaris_Analysis/PS.PA.Rds") ##Total
saveRDS(PS.Norm, "/Users/scast/Desktop/SERGIO/Cursos_Extra/Barcelona_Training/Barcelona_Course/Ascaris_Analysis/PS.Norm.Rds") ##Anormalized
saveRDS(PS.Rare, "/Users/scast/Desktop/SERGIO/Cursos_Extra/Barcelona_Training/Barcelona_Course/Ascaris_Analysis/PS.PA.rare.Rds") #for alpha diversity plots


##Alpha diverisity tables
saveRDS(alphadiv, "/Ascaris_Analysis/alphadiv.rds")
saveRDS(alphadiv.rare, "/alphadiv.PA.rare.rds")


#Check packages 
library(phyloseq)
library(microbiome)
library(tidyverse)
require(ggpubr)
require(RColorBrewer)
require(rstatix)
library(cowplot)
library(gridExtra)
library(grid)
library(ggsci)
library(microbiome)
library(lme4)

##Load data
PS.PA<- readRDS("Ascaris_Analysis/PS.PA.Rds") 
PS.PA.Norm<- readRDS("Ascaris_Analysis/PS.PA.Norm.Rds")
PS.PA.rare<- readRDS("Ascaris_Analysis/PS.PA.rare.Rds") 

##Alpha diversity tables with sample information
alphadiv.PA.rare<- readRDS("Ascaris_Analysis/alphadiv.PA.rare.rds")

##Color palette for compartment and system ##

pal.compartment <- c("Ascaris"="#1B9E77","Feces"= "#D95F02","Intestine"= "#7570B3",
                     "Liver"= "#E7298A","Lung"= "#66A61E")

pal.infection<- c("Infected"= "#E31A1C","Non-Infected" = "#1B9E77")

##Functions
##Find dominant taxa per samples
find.top.asv <- function(x, taxa, num){
  require(phyloseq)
  require(magrittr)

  top.taxa <- tax_glom(x,taxa)
  otu <- as(otu_table(top.taxa), "matrix")
  # transpose if necessary
  if(taxa_are_rows(top.taxa)){otu <- t(otu)}
  otu <- otu_table(otu, taxa_are_rows = F)
  tax <- tax_table(top.taxa)
  # Coerce to data.frame
  n <- as.data.frame(tax)
  n%>%
    rownames_to_column()%>%
    dplyr::rename(ASV = rowname)-> n

  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) 
  j2 <- lapply(j1,'[[',"x") 

  m <- data.frame(unlist(j2))
  m%>%
    rownames_to_column()%>%
    dplyr::filter(unlist.j2.!=0)%>%
    dplyr::mutate(rowname = gsub(".ASV", "_ASV", rowname))%>%
    separate(rowname, sep = "_", c("Replicate","ASV"))%>%
    dplyr::group_by(Replicate)%>%
    slice_max(order_by = unlist.j2., n = num)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")->m

  rm(top.taxa, otu, tax, j1, j2, n)
  return(m)
}

##Get data frame for bar plot at genus level
count.high.genus <- function(x, num){
  require(phyloseq)
  require(magrittr)
  #x is a phyloseq object glomed to Genus
  #num is the threshold of Relative abundance desired
  otu <- as(otu_table(x), "matrix")
  # transpose if necessary
  if(taxa_are_rows(x)){otu <- t(otu)}
  otu <- otu_table(otu, taxa_are_rows = F)
  tax <- tax_table(x)
  # Coerce to data.frame
  n <- as.data.frame(tax)
  n%>%
    rownames_to_column()%>%
    dplyr::rename(ASV = rowname)-> n

  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"x") # select for Names

  m <- data.frame(unlist(j2))
  m%>%
    rownames_to_column()%>%
    dplyr::filter(unlist.j2.!=0)%>%
    dplyr::mutate(rowname = gsub(".ASV", "_ASV", rowname))%>%
    separate(rowname, sep = "_", c("Replicate","ASV"))%>%
    dplyr::group_by(Replicate)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")%>%
    mutate(Main_taxa= Abundance>= num)%>%
    dplyr::mutate(Genus= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(Genus)))%>%
    arrange(Replicate, desc(Genus))->m

  m$Genus[is.na(m$Genus)]<- "Unassigned" ##Change NA's into Unassigned
  m$Species<- NULL

  rm(otu, tax, j1, j2, n)
  return(m)
}

##Transform abundance into relative abundance
Rel.abund_fun <- function(df){
  df2 <- sapply(df, function(x) (x/1E6)*100)
  colnames(df2) <- colnames(df)
  rownames(df2) <- rownames(df)
  df2<- as.data.frame(df2)
  return(df2)
}



#Alpha diversity by organ INF-NI vs Ascaris

#Statistical evaluation of diferences between groups, regarging Shannon Index
#Replace according with analysis

alphadiv.PA.rare%>%
  dplyr::filter(Organ%in%c("Ascaris"))%>%
  mutate(Group_2 = fct_relevel(Group_2, "Larva_InVitro",
                             "LarLiv_4DPI", "LarLun_8DPI"))%>%
  wilcox_test(Shannon ~ Group_2)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Group_2") %>% print(n = 30)#

#Plot
alphadiv.PA.rare%>%
  mutate(Compartment = fct_relevel(Organ, 
                                   "Ascaris", "Liver", "Lung",
                                   "Intestine", "Feces"))%>%
  dplyr::filter(Compartment%in%c("Ascaris", "Intestine"))%>%
  ggplot(aes(x= Compartment, y= Shannon, color= Groups_Inf_NI, fill= Groups_Inf_NI))+
  geom_boxplot(aes(),outlier.shape=NA)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black", "black"))+
  #scale_fill_manual(values = c("#D55E00","#009E73","#E69F00"), labels = c("Infected", "Non infected", "Ascaris"))+
  xlab("")+
  ylab("ASV Diversity (Shannon Index)")+
  labs(fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  scale_y_continuous(limits=c(0, 4.5))


#Plot to compare organ by all days vs Ascaris Stage
#Statistical Test
alphadiv.PA.rare%>%
  dplyr::filter(Group_2%in% c("Naïve Lungs", "Lungs_4DPI",
                              "Lungs_8DPI", "Lungs_14DPI",
                              "Larva_InVitro",
                              "LarLiv_4DPI", "LarLun_8DPI"))%>%
  mutate(Organ = fct_relevel(Organ, "Lung",
                                   "Ascaris"))%>%
  wilcox_test(Shannon ~ Group_2)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Group_2") %>% print(n = 30)#No difference

##Plot
alphadiv.PA.rare%>%
  mutate(Group_2 = fct_relevel(Group_2,
                               "Naïve Liver", "Liver_4DPI",
                               "Liver_8DPI", "Liver_14DPI",
                               "Larva_InVitro",
                               "LarLiv_4DPI", "LarLun_8DPI"))%>%
  mutate(Group = fct_relevel(Group,
                             "2.1. Feces Naïve",
                             "2.2. Feces Naïve",
                             "2.3. Feces Naïve",
                             "2.4. Feces_4DPI",
                             "2.5. Feces_4DPI",
                             "2.6. Feces_4DPI",
                             "2.7. Feces_8DPI",
                             "2.8. Feces_8DPI",
                             "2.9. Feces_8DPI",
                             "3.10. Feces_14DPI",
                             "3.11. Feces_14DPI",
                             "3.12. Feces_14DPI",
                             "4.1. Naïve Intestine",
                             "4.2. Naïve Intestine",
                             "4.3. Intestine 4DPI",
                             "4.4. Intestine 4DPI",
                             "4.5. Intestine 8DPI",
                             "4.6. Intestine 8DPI",
                             "4.7. Intestine 14DPI",
                             "4.8. Intestine 14DPI",
                             "5.1. Naïve Liver",
                             "5.2. Naïve Liver",
                             "5.3. Liver 4DPI",
                             "5.4. Liver 4DPI",
                             "5.5. Liver 8DPI",
                             "5.6. Liver 8DPI",
                             "5.7. Liver 14DPI",
                             "5.8. Liver 14DPI",
                             "6.1. Naïve Lungs",
                             "6.2. Naïve Lungs",
                             "6.3. Lungs 4DPI",
                             "6.4. Lungs 4DPI",
                             "6.5. Lungs 8DPI",
                             "6.6. Lungs 8DPI",
                             "6.7. Lungs 14DPI",
                             "6.8. Lungs 14DPI",
                             "7.1. Larva_InVitro",
                             "8.1. LarLun_8DPI",
                             "7.2. Larva_InVitro",
                             "7.3. LarLiv_4DPI",
                             "7.4. LarLiv_4DPI",
                             "7.5. LarLiv_4DPI",
                             "7.6. LarLiv_4DPI",
                             "7.7. LarLun_8DPI",
                             "7.8. LarLun_8DPI",
                             "7.9. LarLun_8DPI"))%>%
  dplyr::filter(Organ%in%c("Liver", "Ascaris"))%>%
  ggplot(aes(x= Group_2, y= Shannon, color= Group_2, fill= Group_2))+
  geom_boxplot(aes(),outlier.shape=NA)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black", "black",
                                "black", "black", "black",
                                "black", "black", "black"))+
  #"Infected"= "#E31A1C","Non-Infected" = "#1B9E77"
  #scale_fill_manual(values = c("#E31A1C", "#1B9E77","#E69F00"), labels = c("Infected", "Non-Infected", "Ascaris"))+
  xlab("")+
  ylab("ASV Diversity (Shannon Index)")+
  labs(tag= "", fill= "Group_2")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  scale_y_continuous(limits=c(0, 4.5))-> Plot_Organ_Stage

Plot_Organ_Stage




