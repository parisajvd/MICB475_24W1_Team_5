library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)
library(vegan)
library(ggplot2)

#### Load in RData ####
load("../Metadata-based filtering/metadata-based_filtered_phyloseq.RData")

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

ggsave(filename = "Alpha_beta_diversity/Plots/plot_Shannon.png"
       , Shannon_plot
       , height=4, width=6)

ggsave(filename = "Alpha_beta_diversity/Plots/plot_Observed.png"
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

##Started with overall comparisons of the 3 groups

#Bray-curtis#

bc_dm <- distance(mpt_filtered_rare, method="bray")
pcoa_bc <- ordinate(mpt_filtered_rare, method="PCoA", distance=bc_dm)
adonis2(bc_dm ~ `Group`, data=samp_dat_wdiv)

gg_pcoa_bray <- plot_ordination(mpt_filtered_rare, pcoa_bc, color = "Group") +
  labs(col = "Group", title = "Bray-Curtis Dissimilarity, p = 0.019")
gg_pcoa_bray

ggsave("Alpha_beta_diversity/Plots/plot_bray_pcoa.png" 
       , gg_pcoa_bray
       , height=4, width=5)

gg_pcoa_bray_ellipse <- plot_ordination(mpt_filtered_rare, pcoa_bc, color = "Group", title = "Bray-Curtis Dissimilarity, p = 0.019") +
  stat_ellipse(type = "norm")

ggsave("Alpha_beta_diversity/Plots/plot_bray_pcoa_ellipse.png" 
       , gg_pcoa_bray_ellipse
       , height=4, width=5)

#Weighted unifrac

# Use phyloseq to calculate weighted Unifrac distance matrix
dm_unifrac <- UniFrac(mpt_filtered_rare, weighted=TRUE)
adonis2(dm_unifrac ~ `Group`, data=samp_dat_wdiv)

ord.unifrac <- ordinate(mpt_filtered_rare, method="PCoA", distance="wunifrac")
plot_unifrac <- plot_ordination(mpt_filtered_rare, ord.unifrac, color="Group", title = "Weighted UniFrac, p = 0.02")

plot_unifrac_ellipse <- plot_ordination(mpt_filtered_rare, ord.unifrac, color = "Group", title = "Weighted UniFrac, p = 0.02") +
  stat_ellipse(type = "norm")

ggsave("Alpha_beta_diversity/Plots/plot_unifrac.png" 
       , plot_unifrac
       , height=4, width=5)

ggsave("Alpha_beta_diversity/Plots/plot_unifrac_ellipse.png" 
       , plot_unifrac_ellipse
       , height=4, width=5)

# Jaccard is also significant

dm_jaccard <- vegdist(t(otu_table(mpt_filtered_rare)), method="jaccard")
adonis2(dm_jaccard ~ `Group`, data=samp_dat_wdiv)

ord.jaccard <- ordinate(mpt_filtered_rare, method="PCoA", distance="jaccard")
gg_pcoa_jaccard <- plot_ordination(mpt_filtered_rare, ord.jaccard, color = "Group", title = "Jaccard Distance, p = 0.013")

gg_pcoa_jaccard_ellipse <- plot_ordination(mpt_filtered_rare, ord.jaccard, color = "Group", title = "Jaccard Distance, p = 0.013") +
  stat_ellipse(type = "norm")

ggsave("Alpha_beta_diversity/Plots/plot_jaccard.png" 
       , gg_pcoa_jaccard
       , height=4, width=5)

ggsave("Alpha_beta_diversity/Plots/plot_jaccard_ellipse.png" 
       , gg_pcoa_jaccard_ellipse
       , height=4, width=5)


#Unweighted Unifrac
dm_unweighted_unifrac <- UniFrac(mpt_filtered_rare, weighted=FALSE)
adonis2(dm_unweighted_unifrac ~ `Group`, data=samp_dat_wdiv)
  #Significant

ord.unweight <- ordinate(mpt_filtered_rare, method="PCoA", distance="uunifrac")
gg_pcoa_unweight <- plot_ordination(mpt_filtered_rare, ord.unweight, color = "Group", title = "Unweighted UniFrac, p = 0.015")

gg_pcoa_unweight_ellipse <- plot_ordination(mpt_filtered_rare, ord.unweight, color = "Group", title = "Unweighted UniFrac, p = 0.015") +
  stat_ellipse(type = "norm")

ggsave("Alpha_beta_diversity/Plots/plot_unweighted.png" 
       , gg_pcoa_unweight
       , height=4, width=5)

ggsave("Alpha_beta_diversity/Plots/plot_unweighted_ellipse.png" 
       , gg_pcoa_unweight_ellipse
       , height=4, width=5)


#Filter to make pairwise comparisons between Control and other groups

mpt_con_ec <- subset_samples(mpt_filtered_rare, Group != 'TS')
mpt_con_ts <- subset_samples(mpt_filtered_rare, Group != 'EC')
mpt_ec_ts <- subset_samples(mpt_filtered_rare, Group != 'Con')

samp_dat_wdiv_con_ec <- filter(samp_dat_wdiv, Group != 'TS')
samp_dat_wdiv_con_ts <- filter(samp_dat_wdiv, Group != 'EC')
samp_dat_wdiv_ec_ts <- filter(samp_dat_wdiv, Group != 'Con')

#Looking at beta diversity among different pairs using Bray-Curtis

con_vs_ec_dm <- vegdist(t(otu_table(mpt_con_ec)), method="bray")
adonis2(con_vs_ec_dm ~ `Group`, data=samp_dat_wdiv_con_ec)
#Not significant

con_vs_ts_dm <- vegdist(t(otu_table(mpt_con_ts)), method="bray")
adonis2(con_vs_ts_dm ~ `Group`, data=samp_dat_wdiv_con_ts)
#Significant

ec_vs_ts_dm <- vegdist(t(otu_table(mpt_ec_ts)), method="bray")
adonis2(ec_vs_ts_dm ~ `Group`, data=samp_dat_wdiv_ec_ts)
#Significant

##Using the Tobacco column to consider effects of tobacco smoking on diversity and visualize PcoA plots again

#Bray-curtis
adonis2(bc_dm ~ `Tobacco`, data=samp_dat_wdiv)
  #Significant, p = 0.018

gg_pcoa_bray_tobacco_ellipse <- plot_ordination(mpt_filtered_rare, pcoa_bc, color = "Tobacco", title = "Bray-Curtis Dissimilarity, p = 0.018") +
  stat_ellipse(type = "norm")
gg_pcoa_bray_tobacco_ellipse

ggsave("Alpha_beta_diversity/Plots/plot_bray_tobacco_ellipse.png" 
       , gg_pcoa_bray_tobacco_ellipse
       , height=4, width=5)

#Jaccard
adonis2(dm_jaccard ~ `Tobacco`, data=samp_dat_wdiv)
  #Significant, p = 0.024

gg_pcoa_jaccard_tobacco_ellipse <- plot_ordination(mpt_filtered_rare, ord.jaccard, color = "Tobacco", title = "Jaccard Distance, p = 0.024") +
  stat_ellipse(type = "norm")
gg_pcoa_jaccard_tobacco_ellipse

ggsave("Alpha_beta_diversity/Plots/plot_jaccard_tobacco_ellipse.png" 
       , gg_pcoa_jaccard_tobacco_ellipse
       , height=4, width=5)

#Weighted Unifrac
adonis2(dm_unifrac ~ `Tobacco`, data=samp_dat_wdiv)
  #Significant, p = 0.028

gg_pcoa_unifrac_tobacco_ellipse <- plot_ordination(mpt_filtered_rare, ord.unifrac, color = "Tobacco", title = "Weighted UniFrac, p = 0.028") +
  stat_ellipse(type = "norm")
gg_pcoa_unifrac_tobacco_ellipse

ggsave("Alpha_beta_diversity/Plots/plot_weighted_tobacco_ellipse.png" 
       , gg_pcoa_unifrac_tobacco_ellipse
       , height=4, width=5)


#Unweighted Unifrac
dm_unweighted_unifrac <- UniFrac(mpt_filtered_rare, weighted=FALSE)
adonis2(dm_unweighted_unifrac ~ `Tobacco`, data=samp_dat_wdiv)
  #Not significant

gg_pcoa_unweight_ellipse <- plot_ordination(mpt_filtered_rare, ord.unweight, color = "Tobacco", title = "Unweighted UniFrac, p = non-significant") +
  stat_ellipse(type = "norm")
gg_pcoa_unweight_ellipse

ggsave("Alpha_beta_diversity/Plots/plot_unweighted_tobacco_ellipse.png" 
       , gg_pcoa_unweight_ellipse
       , height=4, width=5)


