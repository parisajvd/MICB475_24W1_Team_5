# Team Meeting 3 Oct 18, 2024
proposal: https://docs.google.com/document/d/1aESWHjTN2cO4Zw7szx_nXRlyc-A9cc_rX5Cf7a7IMEU/edit?usp=sharing

## Questions for Ritu
- our dataset used primers 515F and 806R for V4 region, confirm: we can use the silva-138-99-515-806-nb-classifier.qza already on the server?
 - can use this! :)
- Looking at table-no-choloroplast-no-mitochondria.qzv file vs looking at alpha-rarefication.qzv file, what is the best rarefaction depth to choose?
 - better to go with 20,000 because it shows more features (important because we are doing modeling)
 - pick the depth before looking at the no-chloroplast-no-mitochondria file
 -  can go with this for now and then when making our model, we can choose another depth
 -  unless you have a very good data set, some samples will be discarded regardless
  - should ensure that these samples being discarded are not just from one specific category, but otherwise it is fine 
- help with the approach for the modeling (aim 1)
 - may be good to go over references in sex inference paper, easier to understand <-- follow univariant regression analysis for aim 1 
 - do taxonomic classification, train classifier then go to methods, modeling part done in R 
- where to include metadata-based filtering, done in qiime2 or R?
- do we need to specify which type of statistical testing we will use for our diversity metrics, and if so, would it be ANOVA (alpha) and PERMANOVA (beta)?
 - good to state it, but not bound to it
 - just a metric you choose, can say we want to choose these specific metrics
 - if cannot find references for sample type in case of smoking, can say that they found it for something else but we want to see if it applies to smoking 
 

## General Questions for the team
- how to refer to these metadata categories? As "personal factors" "characteristics" etc?
 - need to describe what 'lifestyle and demographics' means in the introduction, don't need to mention all categories
 - stick with "lifestyle and demographics"
 - sample type not really within their control, it is based on where technician samples 
- Finalize a title for the proposal
 - adjust title because right now it seems like we are saying exploring ie. race and smoking on the microbiome
 - exploring the effects of different lifestyles and demographics on microbial diversity on a cohort of e-cigarette users or tobacco smokers
  - don't include any specific variables
 - looking at effects of (BMI/race/sample type), and how they individually affect gut microbiome of each smoker, e-cigarette user and control groups then comparing them
  - compare within each category and across each category   
