##metadata-based filtering after deciding that our project will focus on comparing the fecal microbiome or tobacco smokers vs non-smokers
##filtering parameters based on the phyloseq object generated during the modelling process which showed only 1 vegetarian, 1 vegan and the rest of the participants consuming meat; 2 females and the rest male
#BMI filtered based on guidelines for a healthy and overweight adult BMI (18.5-29.9)

#load libraries
library(phyloseq)
library(dplyr)

#load taxonomy-based and feature-based filtered fecal sample phyloseq object from Modeling (Aim 1)
load("R script/Modelling/fecal_samples_phyloseq.Rdata")

#filter phyloseq object to contain only samples from participants with a BMI between 18.5 and 29.9 (Healthy and Overweight samples), those whose diet is meat, and males
mpt_filtered_phyloseq <- mpt_final %>%
  subset_samples(BMI >= 18.5 & BMI <= 29.9 & Diet == "Meat" & Gender == "Male") 

#remove Height column from metadata
##extract the metadata from the phyloseq object
no_height_metadata <- as(sample_data(mpt_filtered_phyloseq), "data.frame")
##Remove the "Height" column from metadata
no_height_metadata <- no_height_metadata %>%
  select(-Height)
##Reassign the modified metadata back to the phyloseq object
sample_data(mpt_filtered_phyloseq) <- no_height_metadata

#save metadata-based filter phyloseq as .RData file
save(mpt_filtered_phyloseq, file = "R script/Metadata-based filtering/metadata-based_filtered_phyloseq.RData")

