# Presenting our diversity and functional analyses
- 2024-11-22 meeting presentation: https://docs.google.com/presentation/d/1M0-2-cZqUso1gHC1aZEmNixnxrNG9UIhkUna8_c7evA/edit#slide=id.p
# Meeting notes:
#Going over results:
- Modelling showed that Tobacco was the only significant variable, only in Fecal samples
- Alpha diversity: no significance
- In beta diversity, we saw a significance, shown in our presentation bray-curtis. Tobacco yes, no, which one is most significant. It might be abundance so maybe we see it in bray curtis -> p-values 
- For Stats, we need to break it down to 2 groups
- Taxabar plot at the genus level and phylum level-> use both to show the enrichments. Look into some of the genera to see if it is associated with genera
- Deseq2 kinda says the same thing as taxabarplot. Bacteriods go down, Provetella goes up
- Volcano plot of deseq2 has a lot of significance
- Indicator taxa: if changes are in abundance and not in species, it makes sense that we don't get a lot of indicators -> worth putting in supplemental
- Picrust: this shows everything, not just fecal. Abundance data as well as the metadata filtering for fecal -> maybe use qiime2 to filter
- Z more associated (red). Yes vs no
- Our data is very abundance-driven, maybe functional analysis is not the most meaningful -> use Network plot from Hans group, how these abundances affect the interactions with each other. Links them together -> can use phyloseq object for fecal (has everything we need). Circle plot and nodes
- If we don't figure the network analysis out, maybe go with the picrust
- For the network/coordination diagrams: maybe do it on phylum/genus level
- In paper: discuss fecal quite a bit, abundance changes is the main thing.
- Check the literature to see why fecal is changing more than oral.

  # Presentation notes:
  - Titles should be descriptive for each of the results (not just like: beta diversity)
  - For the analyses that are not the main part of the course, try to put enough info for them on the slides.
  - Keep the slides simple. Include enough just to tell the story. 
