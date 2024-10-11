# Meeting 2 Friday Oct 11, 2024

## Pre Meeting Notes
- will be using UK Vaping Dataset 
- doing a modeling project, everyone: look at the papers Ritu sent for Friday
- make copy of dataset into data folder 
- systems:
  - qiime2: denoising, declustering (up to table and rep-seq file), then we will move to R
- Proposal logistics:
  - split up work: writing vs coding?
    - Make a more detailed plan for the proposal this week !!
    - parisa: coding up to table and rep-seq part, then we need Ritu's help
    - start prepping an introduction by tuesday/wednesday preferably
    - ensure we are referencing alternate grading rubric for proposal

# Meeting Notes
- google docs file for collaborating on proposal: https://docs.google.com/document/d/1aESWHjTN2cO4Zw7szx_nXRlyc-A9cc_rX5Cf7a7IMEU/edit?usp=sharing
- alternative proposal is maybe not the most fitting for ours, the original proposal might be better
  - updates: USING ORIGINAL PROPOSAL RUBRIC INSTEAD OF THE ALTERNATIVE RUBRIC !! 
- use a google document for better collaboration, jot down points/an outline for each section
- finish writing draft potentially by Thursday and then Ritu can look at it and give us feedback before submitting, can add her to the Google document  
  - if it doesn't work, let her know in advance
- Proposal:
  - use the same formatting of the aims from the modeling papers sent
  - can pick like 3 variables that we find interesting, say we are doing modeling project to see whether these are significant and then look     into the literature to talk about what has been found regarding their microbiome
  - processing overview (Parisa):
    - make a copy into data folder
    - manifest file step
    - demultiplexing, look at table.qzv file on qiimeview
    - denoising and declustering step (no preference over DADA2 vs deblur) = ASV/OTU file
    - look at qzv of those, trimming
    - final files: rep-seq, feature table or proposal
    - phylo tree and taxonomic classification after proposal
    - then modeling (after taxonomic classification): there is an R code in the papers, will run into every category, calculate alpha and beta diversity, give us p-value and then we need to decide
    - when we mention dataset, cite the age group
    - AIMS: 2 main ones
      - 1. modeling
        2. taxonomic analysis - subsections: diversity metrics (JUST STICK TO ONE EITHER TAXONOMIC OR FUNCTIONAL, BUT CAN CHANGE LATER)
    - research question: which factors affect the microbiome
    - approach: in form of table
- Next Steps:
  - modeling gives us p-value to find significant factors
  - depending on how many are significant, focus on those particular ones and focus taxonomic classification with just that particular variable
    - redo taxonomic analysis step with just chosen significant variables  
  - can do multiple variables if there is time or can focus on one
