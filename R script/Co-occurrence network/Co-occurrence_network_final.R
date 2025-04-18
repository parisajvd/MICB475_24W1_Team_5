### Lets try building some co-occurrence networks
library(SpiecEasi)
library(igraph)
library(phyloseq)
library(ggplot2)
library(viridisLite)
library(RColorBrewer)
library(dplyr)

load("../Metadata-based filtering/metadata-based_filtered_phyloseq.RData")

# Step 0: Filter to include only samples with 'Tobacco' equal to 'No/Yes'
tobacco_no <- subset_samples(mpt_filtered_phyloseq, Tobacco == "No")
tobacco_yes <- subset_samples(mpt_filtered_phyloseq, Tobacco == "Yes")


# Step 1: Filter by genus with a minimum abundance threshold
filter_by_genus <- function(physeq_obj, min_count = 1, min_samples = 1) {
  # Aggregate to genus level
  physeq_genus <- tax_glom(physeq_obj, taxrank = "Genus")
  
  # Remove low-abundance genera
  physeq_genus <- prune_taxa(taxa_sums(physeq_genus) > min_count, physeq_genus)
  
  # Keep only genera that appear in at least 'min_samples' samples
  physeq_genus <- filter_taxa(physeq_genus, function(x) sum(x > 0) >= min_samples, prune = TRUE)
  
  return(physeq_genus)
}

# Step 2: Filter the dataset at the genus level using the function above
tobacco_no_filtered <- filter_by_genus(tobacco_no)



genera_to_label <- c(
  "g_Bacteroides", "g__Blautia", "g_Bacteroides.1", "g__Dialister", "g__Faecalibacterium",
  "g_Lachnoclostridium", "g_Barnesiella", "g_Butyricicoccus", "g_[Eubacterium]_coprostanoligenes_group",
  "g_Sutterella", "g_Agathobacter", "g_Bacteroides.2", "g_Lactococcus", "g_Bacteroides.3",
  "g_Incertae_Sedis", "g_Alistipes", "g_Barnesiella.1", "g_[Ruminococcus]_gnavus_group",
  "g_Christensenellaceae_R-7_group", "g_NK4A214_group", "g_Colidextribacter", "g_[Eubacterium]_eligens_group",
  "g_Streptococcus", "g_Clostridium_sensu_stricto_1", "g_UBA1819", "g_UCG-002", "g_Marvinbryantia",
  "g_Incertae_Sedis.1", "g_Turicibacter", "g_Bifidobacterium", "g_[Eubacterium]_hallii_group",
  "g_[Eubacterium]_siraeum_group", "g_Veillonella", "g_Collinsella", "g_Moryella", "g_Ruminococcus",
  "g_Escherichia-Shigella", "g_Clostridia_vadinBB60_group", "g_Marvinbryantia.1", "g_Bilophila",
  "g_[Eubacterium]_coprostanoligenes_group.1", "g_[Eubacterium]_siraeum_group.1", "g_UCG-002.1",
  "g_Clostridia_UCG-014", "g_Sutterella.1", "g_[Ruminococcus]_gauvreauii_group", "g_Desulfovibrio",
  "g_Subdoligranulum", "g_Terrisporobacter", "g_Blautia.1", "g_Coprococcus", "g_Dialister.1",
  "g_Blautia.2", "g_Mogibacterium", "g_Bacteroides.4", "g_Prevotella", "g_[Ruminococcus]_torques_group",
  "g_Holdemanella", "g_Phascolarctobacterium", "g_Clostridia_UCG-014.1", 
  "g_[Eubacterium]_coprostanoligenes_group.2", "g_Prevotella.1", "g_Prevotella.2"
)


generate_network_igraph <- function(physeq_obj, seed = 42, genera_to_label) {
  # Step 1: Apply data transformations and generate SpiecEasi network
  print("Applying data transformations...")
  otu_data <- as.matrix(otu_table(physeq_obj))
  se.mb <- spiec.easi(otu_data, method = 'mb', pulsar.params = list(rep.num = 40))
  print("Selecting model with pulsar using stars...")
  print("Fitting final estimate with mb... done")
  print("SpiecEasi analysis completed")
  
  # Step 2: Convert to igraph object
  adj_matrix <- getRefit(se.mb)
  ig_network <- adj2igraph(adj_matrix, vertex.attr = list(name = taxa_names(physeq_obj)))
  print("igraph Network Created")
  
  # Step 3: Add taxonomic information to vertices
  taxa_info <- as.data.frame(tax_table(physeq_obj))
  V(ig_network)$Genus <- taxa_info$Genus[match(V(ig_network)$name, rownames(taxa_info))]
  V(ig_network)$Phylum <- taxa_info$Phylum[match(V(ig_network)$name, rownames(taxa_info))]
  
  # Step 4: Check abundance data
  otu_abundance <- taxa_sums(physeq_obj)
  valid_names <- intersect(V(ig_network)$name, names(otu_abundance))
  otu_abundance <- otu_abundance[valid_names]
  
  # Step 5: Assign sizes to vertices
  V(ig_network)$size <- sapply(V(ig_network)$name, function(x) {
    value <- otu_abundance[x]
    if (is.list(value)) value <- unlist(value)[1]
    if (!is.null(value) && !is.na(value) && is.numeric(value)) {
      return(log10(value + 1) * 4)
    } else {
      return(1)
    }
  })
  
  # Ensure sizes are numeric
  if (any(sapply(V(ig_network)$size, is.list))) {
    V(ig_network)$size <- sapply(V(ig_network)$size, function(x) if (is.list(x)) unlist(x)[1] else x)
  }
  V(ig_network)$size <- as.numeric(V(ig_network)$size)
  
  # Step 6: Generate colors for the phyla
  unique_phyla <- unique(V(ig_network)$Phylum)
  phylum_colors <- rainbow(length(unique_phyla))
  names(phylum_colors) <- unique_phyla
  V(ig_network)$color <- phylum_colors[V(ig_network)$Phylum]
  
  # Step 7: Assign label and border colors for genera of interest and their neighbors
  V(ig_network)$label_color <- "black"
  V(ig_network)$border_color <- "black"
  
  # Step 8: Label nodes of interest in red and their neighbors in black
  V(ig_network)$label <- NA
  
  # Label the specified genera in red
  V(ig_network)$label[V(ig_network)$Genus %in% genera_to_label] <- V(ig_network)$Genus[V(ig_network)$Genus %in% genera_to_label]
  V(ig_network)$label_color[V(ig_network)$Genus %in% genera_to_label] <- "red"
  
  # Find nodes directly connected to the specified genera
  connected_to_labeled <- unique(unlist(sapply(genera_to_label, function(g) {
    if (g %in% V(ig_network)$Genus) {
      neighbor_vertices <- neighbors(ig_network, V(ig_network)$Genus == g)
      return(V(ig_network)[neighbor_vertices]$name)
    }
  })))
  
  # Label the neighbors in black
  V(ig_network)$label[V(ig_network)$name %in% connected_to_labeled] <- V(ig_network)$Genus[V(ig_network)$name %in% connected_to_labeled]
  V(ig_network)$label_color[V(ig_network)$name %in% connected_to_labeled] <- "black"
  
  # Step 9: Ensure degrees are numeric
  node_degrees <- degree(ig_network)
  if (any(sapply(node_degrees, is.list))) {
    node_degrees <- sapply(node_degrees, function(x) if (is.list(x)) unlist(x)[1] else x)
  }
  node_degrees <- as.numeric(node_degrees)
  
  # Remove nodes with zero degree (unconnected)
  ig_network <- delete_vertices(ig_network, V(ig_network)[node_degrees == 0])
  
  # Step 10: Plot the network
  set.seed(seed)
  layout <- layout_with_fr(ig_network)
  
  plot(
    ig_network,
    layout = layout,
    vertex.label = V(ig_network)$label,
    vertex.label.color = V(ig_network)$label_color,
    vertex.label.cex = 0.6,
    vertex.size = V(ig_network)$size,
    vertex.color = V(ig_network)$color,
    vertex.frame.color = V(ig_network)$border_color,
    edge.width = 0.5,
    edge.color = "gray50",
    main = "Microbial Network by Phylum"
  )
  
  # Add legend
  legend(
    "topleft",
    legend = names(phylum_colors),
    col = phylum_colors,
    pch = 19,
    cex = 0.8,
    title = "Phylum"
  )
  
  print("Plotting completed")
}
# Run the function
generate_network_igraph(tobacco_no_filtered, genera_to_label = genera_to_label)
