library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)
library(vegan)
library(ggplot2)

#### Load in RData ####
load("Alpha:beta diversity/metadata-based_filtered_phyloseq.RData")

rarecurve(t(as.data.frame(otu_table(mpt_filtered_phyloseq))), cex=0.1)
mpt_filtered_rare <- rarefy_even_depth(mpt_filtered_phyloseq, rngseed = 1, sample.size = 12000)

#980 OTUs removed after random subsampling

#### Alpha diversity ######

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(mpt_filtered_rare)), phy_tree(mpt_filtered_rare),
                 include.root=F) 

# add PD to metadata table
sample_data(mpt_filtered_rare)$PD <- phylo_dist$PD

#Add alpha diversity metrics into the phyloseq
samp_dat_wdiv <- data.frame(sample_data(mpt_filtered_rare), estimate_richness(mpt_filtered_rare))

#Kruskal-wallis on alpha diversity metrics including PD

Shannon_plot <- ggplot(samp_dat_wdiv, aes(x=`Group`, y=Shannon)) +
  geom_boxplot() +
  geom_point()

Observed_plot <- ggplot(samp_dat_wdiv, aes(x=`Group`, y=Observed)) +
  geom_boxplot() +
  geom_point()

ggsave(filename = "plot_Shannon.png"
       , Shannon_plot
       , height=4, width=6)

ggsave(filename = "plot_Observed.png"
       , Observed_plot
       , height=4, width=6)

kruskal_obs <- kruskal.test(Observed ~ `Group`, data = samp_dat_wdiv)

kruskal_Shannon <- kruskal.test(Shannon ~ `Group`, data = samp_dat_wdiv)

kruskal_Chao1 <- kruskal.test(Chao1 ~ `Group`, data = samp_dat_wdiv)

kruskal_ACE <- kruskal.test(ACE ~ `Group`, data = samp_dat_wdiv)

kruskal_Simpson <- kruskal.test(Simpson ~ `Group`, data = samp_dat_wdiv)

kruskal_Fisher <- kruskal.test(Fisher ~ `Group`, data = samp_dat_wdiv)

kruskal_PD <- kruskal.test(PD ~ `Group`, data = samp_dat_wdiv)

#None are significant

###Beta diversity####

#Bray-curtis#

bc_dm <- distance(mpt_filtered_rare, method="bray")
pcoa_bc <- ordinate(mpt_filtered_rare, method="PCoA", distance=bc_dm)

gg_pcoa_bray <- plot_ordination(mpt_filtered_rare, pcoa_bc, color = "Group") +
  labs(col = "Group")
gg_pcoa

ggsave("plot_bray_pcoa.png" 
       , gg_pcoa_bray
       , height=4, width=5)

adonis2(bc_dm ~ `Group`, data=samp_dat_wdiv)

gg_pcoa_bray_ellipse <- plot_ordination(mpt_filtered_rare, pcoa_bc, color = "Group") +
  stat_ellipse(type = "norm")

ggsave("plot_bray_pcoa_ellipse.png" 
       , gg_pcoa_bray_ellipse
       , height=4, width=5)

##Bray-curtis distance is significant

#Weighted unifrac

# Use phyloseq to calculate weighted Unifrac distance matrix
dm_unifrac <- UniFrac(mpt_filtered_rare, weighted=TRUE)

ord.unifrac <- ordinate(mpt_filtered_rare, method="PCoA", distance="unifrac")
plot_unifrac <- plot_ordination(mpt_filtered_rare, ord.unifrac, color="Group")

adonis2(dm_unifrac ~ `Group`, data=samp_dat_wdiv)

plot_unifrac_ellipse <- plot_ordination(mpt_filtered_rare, ord.unifrac, color = "Group") +
  stat_ellipse(type = "norm")

ggsave("plot_unifrac.png" 
       , plot_unifrac
       , height=4, width=5)

ggsave("plot_unifrac_ellipse.png" 
       , plot_unifrac_ellipse
       , height=4, width=5)

# Jaccard is also significant

dm_jaccard <- vegdist(t(otu_table(mpt_filtered_rare)), method="jaccard")
adonis2(dm_jaccard ~ `Group`, data=samp_dat_wdiv)

