# Results of Diversity analyses
- Florrie: DEseq2: compares which ASVs are upregulated in tobacco users vs non-users.
  -   Volcano plot: Blue is significant (Florrie has plotted both p-values 0.01 and 0.05). Fold change greater than 0, upregulated, less than 0, downregulated. The corners are most significantly different -> to look at the relative abundance.
  -   Bar plots: takes significant ASVs and plots it. the higher the bar the more significant it is.
### Note: We should be consistent with the paths used to load our files. So that they can load the files from the GitHub 

- Monica: Alpha-beta diversity: Alpha diversity: no significance in any of them (bar plot, Shannon) an this is consistent with the main paper. Beta diversity: we found significance, bray Curtis and the pcoA plot ellipses showed a bit of significance. Unifrac also shows a bit of significance (pcoa plot ellipses visually helps)
### Should we keep the metadata groups as: control, EC users, and tobacco users? 
Maybe this would be good for the start of the paper. But then further in the paper as we found tobacco as significant, we can only keep control and tobacco users.

- Raychal: Indicator genera analysis: genera that could be indicative of tobacco vs non-smokers based on the ASVs. 4 that were classified as non-smoking vs tobacco (shown in a table). The first 3 in non-tobacco,  Last one in the table was for tobacco smokers.
  - Bar graph: height is the stat -> closer to 1 means it is more associated with one group than the other. Tobacco-associated one is the most significant. Protevella genus is the one that can indicate the tobacco smoker!
  - Maybe stick to p-value of 0.05 for all of the analysis.

-  Parisa: Taxonomy / taxabar plot analysis: When resolved to genus level, the barplot gets a bit complicated and messy. But looks like we have more of the purple color genra in tobacco smokers vs non. -> associated with Protevella again! As taxabar plots are to accompany DEseq, this might still be useful, maybe need to narrow down the genus, or just make the graph visually better with ggplot
  - Also resolved to phylum level -> doesn't look too different between tobacco vs non, might not be the most useful for our paper because all other analyses are resolving to genus level.
