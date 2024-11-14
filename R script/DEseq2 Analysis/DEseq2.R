#load libraries 
library(tidyverse)
library(phyloseq)
library(DESeq2)

#load metadata filtered phyloseq object
load("~/Desktop/MICB 475 R/metadata-based_filtered_phyloseq.RData")

#convert phyloseq to deseq
mpt_deseq <- phyloseq_to_deseq2(mpt_filtered_phyloseq, ~`Tobacco`)

#add '1' count to all reads
mpt_plus1 <- transform_sample_counts(mpt_filtered_phyloseq, function(x) x+1)
mpt_deseq <- phyloseq_to_deseq2(mpt_plus1, ~`Tobacco`)
DESEQ_mpt <- DESeq(mpt_deseq)
res <- results(DESEQ_mpt, tidy=TRUE)
View(res)

#visualize Volcano plot
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

#color by whether it is significant and has a large fold change
vol_plot <- res %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot

#save Volcano plot
ggsave(filename="vol_plot.png",vol_plot)

# To get table of results
sigASVs <- res %>% 
  filter(padj<0.05 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)

# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
mpt_DESeq <- prune_taxa(sigASVs_vec,mpt_filtered_phyloseq)
sigASVs <- tax_table(mpt_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

#Visualize barplot
bar_plot <- ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
bar_plot

#save barplot
ggsave(filename="bar_plot.png",bar_plot)
