
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


#Beta Diversity Analysys

#plot using ggplot
##extract sample id from the phyloseq object
##extract sample id

otu_tab <- t(abundances(PS))
rare_curve <- vegan::rarecurve(otu_tab, col = colors_rare, step = 50, ylab = "Amplicon sequence variants",
                               label = FALSE, sample = min(rowSums(otu_tab),cex = 0.6))


PS_rare <- rarefy_even_depth(PS, sample.size =  10000) 

dist_bray = phyloseq::distance(PS_rare, method="bray")

#to get help about a function of a package
#?phyloseq::distance

##Perform Prinicipal Coordinates Analysis (PCoA) on the rarefied phyloseq object
ordination_bray = ordinate(PS_rare, method="PCoA", distance=dist_bray)


sampleid_PS_rare <- row.names(sample_data(PS_rare))

##add sample id to the data frame and extract four axis
df_PS_rare_bray <- rbind(data.frame(ordination_bray$vectors[,1:4],
                                          sample_name = sampleid_PS_rare))

##Get the metadata out as separate object
PS.meta <- meta((PS_rare))

##Check variables
PS.meta$Organ

# Add the rownames as a new colum for easy integration later.
PS.meta$sample_name <- rownames(PS.meta)


#merge sample data and ordination
df_PS_rare_bray_complete <- merge(df_PS_rare_bray,
                                           PS.meta, by = 'sample_name')

# Add percentages to the data frame
df_PS_rare_bray_complete$Axis.1_percent <- df_PS_rare_bray_complete$Axis.1 * 100
df_PS_rare_bray_complete$Axis.2_percent <- df_PS_rare_bray_complete$Axis.2 * 100

##Plot PCoA plot using ggplot2
ggplot(data = df_PS_rare_bray_complete, aes(Axis.1, Axis.2,color = factor(Organ))) +
  theme_classic() +
  geom_point() +
  stat_ellipse(aes(group = factor(Organ))) +
  theme(legend.position = "right",  # Remove the automatic legend positioning
        legend.justification = "center",
        legend.title = element_text(hjust = 0.5),
        legend.position.inside = c(0.1, 0.85),  # Adjust these values to move the legend
        legend.margin = margin(t = -0, r = 0, b = 0, l = -20))
ggsave("figure_output/bmi_braycurtis.tiff", bg = "white", dpi = 600)


per_brayasv_bmi <- adonis2(dist_bray ~ Group_2, data = PS.meta, permutations = 999)
per_brayasv_bmi


#Abundance analysis

desired_order_phylum <- c("Ascaris",
                   "Lung",
                   "Liver",
                   "Intestine",
                   "Feces")

desired_order <- c("LarLun_8DPI",
                   "LarLiv_4DPI",
                   "Larva_InVitro",
                   "Lungs_14DPI",
                   "Lungs_8DPI",
                   "Lungs_4DPI",
                   "Naïve Lungs",
                   "Liver_14DPI",
                   "Liver_8DPI",
                   "Liver_4DPI",
                   "Naïve Liver",
                   "Intestine 14DPI",
                   "Intestine 8DPI",
                   "Intestine 4DPI",
                   "Naïve Intestine",
                   "Feces_14DPI",
                   "Feces_8DPI",
                   "Feces_4DPI",
                   "Feces Naïve")


# Plot taxa abundance
PS_rare %>% tax_fix(unknowns = c("NA")) %>%
    ps_filter(Organ != "") %>%
    comp_barplot(tax_level = "Phylum", n_taxa = 7) +
    #facet_wrap(~Organ, scales= "free", ncol=1) #+
    coord_flip()

# To get names no "Other"
taxa_to_keep <- rownames(tax_table(PS_rare))[!grepl("Other", rownames(tax_table(PS_rare)))]

# Filter phyloseq
PS_rare_filtered <- subset_taxa(PS_rare, taxa_names(PS_rare) %in% taxa_to_keep)

# 15 most abundant taxa
top_taxa <- names(sort(taxa_sums(PS_rare_filtered), decreasing = TRUE)[1:20])

# Filter low taxa
PS_rare_filtered <- prune_taxa(top_taxa, PS_rare_filtered)

#General (include other)
PS_rare %>% 
  tax_fix(unknowns = c("NA")) %>%  # Avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Organ") %>% 
  comp_barplot(tax_level = "Phylum", n_taxa = 10, bar_width = 0.8) + 
  #facet_wrap(~ Group_2, scales = "free_y", ncol = 1) +  # Custom order for Group_2
  coord_flip() + labs(x = NULL, y = NULL) +
  scale_x_discrete(limits = desired_order_phylum)

#"Other" filtered
PS_rare_filtered %>% 
  tax_fix(unknowns = c("NA")) %>%  # Avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Group_2") %>% 
  comp_barplot(tax_level = "Class", n_taxa = 20, bar_width = 0.8) + 
  #facet_wrap(~ Group_2, scales = "free_y", ncol = 1) +  # Custom order for Group_2
  coord_flip() + labs(x = NULL, y = NULL) +
  scale_x_discrete(limits = desired_order)


#differential abundance analysis
##Show significantly different taxa at ASV level between positive and negative blasto sample

#By organ
lefse_model <- run_lefse(
  PS_rare,
  taxa_rank = "all",
  wilcoxon_cutoff = 0.01,
  group = "Organ",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 4.5
)


#By INF-NI Group
lefse_2 <- run_lefse(
  PS_rare,
  taxa_rank = "all",
  wilcoxon_cutoff = 0.01,
  group = "Groups_Inf_NI",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 4
)



#plot lefese at the asv level

#barplot
plot_ef_bar(lefse_2)
plot_ef_do(lefse_model)

#heatmap
plot_heatmap(lefse_2, transform = "log10p", group = "Groups_Inf_NI",
             cluster_marker = T, cluster_sample = F)


#Abundance (Boxplot)
p_abd <- plot_abundance(lefse_2, group = "Groups_Inf_NI")
p_abd

#Cladogram
plot_cladogram(lefse_model, color = c(Ascaris = "orange", Feces = "brown",
                                       Intestine = "darkgreen" , Liver = "blue" , Lung = "purple" )) +
  theme(plot.margin = margin(0, 0, 0, 0))

#Cladogram
plot_cladogram(lefse_2, color = c("Non-Infected" = "darkgreen", Ascaris = "orange")) +
  theme(plot.margin = margin(0, 0, 0, 0))


##Core microbiome analysis
##
###Infected - Non Infected - Ascaris
infst<- c("Infected", "Non-Infected", "Ascaris")

list_core <- c() # an empty object to store information

for (n in infst){ # for each variable n in InfectionStatus
  
  tmp<- subset_samples(PS_rare, Groups_Inf_NI==n)
  tmp <- microbiome::transform(tmp, "compositional")
  core_m <- core_members(tmp, 
                         detection = 0.0001, # 0.001 in atleast 90% samples 
                         prevalence = 0.3)
  
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each InfectionStatus.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}


#Venn diagram
plot(eulerr::venn(list_core),
     fills = c("seagreen3", "dodgerblue","salmon2"),
     alpha= 0.8,
     labels = c("Infected", "Non infected", "Ascaris"))-> D

cowplot::plot_grid(
  D, labels = "", 
  label_fontfamily = "sans",
  label_fontface = "plain",
  label_size = 20)-> venn_diagram

venn_diagram


ggsave(file = "core_taxa.pdf", plot = venn_diagram, width = 8, height = 8, dpi = 600)
ggsave(file = "core_taxa.png", plot = venn_diagram, width = 8, height = 8, dpi = 600)
ggsave(file = "core_taxa.svg", plot = venn_diagram, width = 8, height = 8, dpi = 600)

# get the taxonomy data
tax.mat <- tax_table(PS_rare)
tax.df <- as.data.frame(tax.mat)

# add the OTUs to last column
tax.df$ASV <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.inf <- dplyr::filter(tax.df, rownames(tax.df) %in% list_core[["Infected"]])
core.taxa.Ninf <- dplyr::filter(tax.df, rownames(tax.df) %in% list_core[["Non-Infected"]])
core.taxa.ascaris <- dplyr::filter(tax.df, rownames(tax.df) %in% list_core[["Ascaris"]])

core.asv<- union(list_core[["Infected"]], list_core[["Non-Infected"]])
core.asv <- dplyr::filter(tax.df, rownames(tax.df) %in% list_core[["Infected"]])

core.asv<- union(core.asv, list_core[["Ascaris"]])

##Store "core" ASVs
write.csv(core.taxa.inf, "Core_ASVs_Inf.csv")
write.csv(core.taxa.Ninf, "Core_ASVs_NonInf.csv")
write.csv(core.taxa.ascaris, "Core_ASVs_Ascaris.csv")
write.csv(core.asv, "Core_ASVs.csv")

# common Inf-Ascaris
common_ascaris_inf <- inner_join(core.taxa.ascaris, core.taxa.inf, by = "ASV")

#Differente common_ascaris_inf vs NI
common_inf_ascaris <- anti_join(common_ascaris_inf, core.taxa.Ninf, by = "ASV")

print(common_inf_worm)
n_common_inf_worm <- nrow(common_inf_worm)

print(n_common_inf_worm)
#8
