#Indicator Genera Analysis using the metadata-based filtering phyloseq object to determine ASVs and genera that are unique to tobacco smokers or non-smokers 
#scoring ASVs on prevalence and abundance

#load libraries 
library(tidyverse)
library(phyloseq)
library(indicspecies)
library(ggplot2)

#load metadata-based filtered phyloseq object
load("../Metadata-based filtering/metadata-based_filtered_phyloseq.RData")

#### Indicator Taxa Analysis ####
# group/glom data to Genus level, do not remove NAs
mpt_genus <- tax_glom(mpt_filtered_phyloseq, "Genus", NArm = FALSE)
#convert counts in phyloseq table to relative abundance
mpt_genus_RA <- transform_sample_counts(mpt_genus, fun=function(x) x/sum(x))
#use multipatt to cluster samples by group using, transpose OTU table, cluster based on tobacco usage
isa_mpt <- multipatt(t(otu_table(mpt_genus_RA)), cluster = sample_data(mpt_genus_RA)$`Tobacco`)
#view results and filter based on significance of p-value less than or equal to 0.05
summary(isa_mpt)
#tabulate by combining taxonomic information to indicator species list by extracting tax_table as data frame, separate ASVs as own column
taxtable <- tax_table(mpt_filtered_phyloseq) %>% as.data.frame() %>% rownames_to_column(var="ASV")
#create comprehensive summary table 
isa_final <- isa_mpt$sign %>%
  rownames_to_column(var="ASV") %>% #convert ASVs to its own column
  left_join(taxtable) %>% #left join taxtable based on ASVS
  filter(p.value<0.05) #filter for p < 0.05

#save summary table as csv file
write.csv(isa_final, file = "indicator_genera_analysis.csv", row.names = TRUE)

#### Bar Plot Visualization ####
#reorder the Genus factor based on the 'stat' variable in ascending order
isa_final$Genus <- factor(isa_final$Genus, levels = isa_final$Genus[order(isa_final$stat, decreasing = FALSE)])

#create a new column to label the groups: "Non-smokers" or "Tobacco smokers" and label first 3 rows as "Non-smokers" and last row as "Tobacco Smokers"
isa_final$group_label <- c(rep("Non-smokers", 3), "Tobacco smokers", rep(NA, nrow(isa_final) - 4))

#visualize indicator genera using a bar plot
indicator_genera_plot <- ggplot(isa_final, aes(x = Genus, y = stat, fill = Genus)) + #plot Genus on x-axis, the stat value (association value) on y-axis and color the bars based on Genus
  geom_bar(stat = "identity") +  #use the stat values as the height of the bars
  geom_text(aes(label = paste("p =", format(p.value, digits =4))), vjust = 1.5) + #add p-value above the bars, rounded to 3 decimal places, adjust positioning of label
  geom_text(aes(label = group_label), vjust = -0.5, color = 'black', size = 4, fontface = "italic") + #add individual labels above each bar with respective smoking status, adjust position, size and font
  labs(title = "Indicator Genera Analysis", x = "Genus", y = "Indicator Statistic (Stat)") + #label the graph
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #adjust the positioning of the x-axis labels so everything can be read

#save the bar graph
ggsave(indicator_genera_plot, file = "indicator_genera_analysis_plot.png", width = 8, height = 7)
