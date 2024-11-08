# This script is designed to rank multiple variables in terms of their contribution to microbial diversity using a PERMANOVA and looking at the
# R- squared value and pvalue. The R-squared statistic will indicate how much any one variable contributes to driving diversity.


#import libraries 
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)

### Create phyloseq object ###

#Load in the metadata,  OTU table, taxonomy file, and phylogenetic tree

metafp <- "uk_metadata.tsv"
meta <- read_delim(metafp, delim="\t")
class(meta)

saliva_meta <- filter(meta, SampleType == "Saliva")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read_tree(phylotreefp)

#Adjust files to be read into a phyloseq object. Make the phyloseq object.

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(saliva_meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- saliva_meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
mpt <- phyloseq(OTU, SAMP, TAX, phylotree)



#Filter for Chloroplast/Mitochondrial sequences
mpt_filt <- subset_taxa(mpt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
#Filter to remove ASV's with counts fewer than 5
mpt_filt_nolow <- filter_taxa(mpt_filt, function(x) sum(x)>5, prune = TRUE)
#Filter to remove samples with fewer than 100 reads
mpt_final <- prune_samples(sample_sums(mpt_filt)>100, mpt_filt_nolow)

#Rarefy samples to a depth of 20034.
rarecurve(t(as.data.frame(otu_table(mpt_final))), cex=0.1)
mpt_rare <- rarefy_even_depth(mpt_final, rngseed = 1, sample.size = 25000)

#2 samples removed due to rarefaction, 610 OTUs no lonfer in any sample after random subsampling.

#Save rarefied phyloseq object
save(mpt_rare, file="mpt_rare_saliva.RData")

###Preparing the modelling table###

#Take out the ASV matrix from the phyloseq object and convert it to a dataframe.
ASV_table = otu_table(mpt_rare)
ASV_table = as.data.frame(otu_table(mpt_rare))#convert to data.frame 

#Take out the metadata from the phyloseq object and convert it to a dataframe.
meta = sample_data(mpt_rare)
meta_df <- data.frame(meta)

#Filter the metadata to remove columns that are not relevant

meta_filt <- meta_df [, c("Age", "Race", "Gender", "Diet", "Ecig", "Tobacco", "Nicotine", "Marijuana", 
                          "Alcohol","CO_ppm", "CO_percent", "Weight", "BMI")] 

MF = colnames(meta_filt)

adonis.res = list()   #Build an empty list that will be filled up by the loop

for (i in 1:length(MF)){ #COMPLETE the for loop argument. You need to to loop through as many variables that are present in "nutrients". Use a range, IE (1:?)
  print(i)#Printing i to keep track of loop progress
  meta_no_missing_data <- meta_filt[complete.cases(meta_filt[, i]), ]#Remove the rows in metadata that contain missing data for the i'th variable
  
  samples = rownames(meta_no_missing_data) #Create a vector of all the samples in the metadata after removing NA's
  ASV_mat = ASV_table[,samples] #Filter the ASV_table to remove the same individuals that were filtered out for NA values 
  #This is important because we need the number of individuals represented in the ASV table and metadata to be the same.
  
  ASV_mat = t(ASV_mat +1) # Transpose the ASV matrix so that the sample IDs become rownames. This is required for  adonis2()
  dis = vegdist(ASV_mat, method = "bray",diag =T, upper = T) #Create the distance matrix based on a bray-curtis dissimilarity model
  adonis.res[[i]] = vegan::adonis2(as.formula(paste("dis~",MF[[i]],sep = "")), data = meta_no_missing_data) #This line runs the PERMANOVA test and puts it into the empty list we made above.
}

#Create an empty table to import the R-squared and pvalue.
result = matrix(NA, nrow = length(MF), ncol =2)

#Loop though each variable and generate a table that can be plotted.
for (i in 1:length(MF)){ #This for loop argument will be the same as before
  
  result[i,1] = adonis.res[[i]][1,3] #Grab the R-squared statistic
  result[i,2] = adonis.res[[i]][1,5]#Grab the pvalue
  
}

rownames(result) = c(MF)   #Convert the rowmanes to variables 
colnames(result) = c("R2", "Pvalue")#Change the column names TO "R2" AND "Pvalue"
result = data.frame(result, stringsAsFactors = F) #Convert it to a data.frame (easiest to work with when plotting)
result$Padjust = p.adjust(result$Pvalue, method = "fdr") #Generate an adjusted pvalue to correct for the probability of false positives
result$Factor =  c("Age", "Race", "Gender", "Diet", "Ecig", "Tobacco", "Nicotine", "Marijuana", 
                   "Alcohol","CO_ppm", "CO_percent", "Weight", "BMI")
#Create another column with variable names
View(result)

###############################PLOTTING

#Try and generate the plot yourself using ggplot.
#Here is a skeleton of the plotting code. I wrote "FILLOUT" where information needs to be added.

#Also, filter the results table to only include significant variables with a pvalue<0.05

result_filtered_Padjust = subset(result, Padjust < 0.05)#Write solution here

#####Saving######
save(result, file = "Modelling_result_saliva.Rdata")
save(result_filtered_Padjust, file = "Modelling_result_filtered_saliva.Rdata")

