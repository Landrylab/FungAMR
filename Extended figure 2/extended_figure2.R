### Created on 31-07-2024
### @author: aliciapageau

library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(edgebundleR)
library(ggplot2)


#### Load functions ####
abbreviate_species <- function(name) {
  # Abbreviate Species Name
  #
  # This function abbreviates a species name by keeping the first letter of the genus 
  # and the full species name, separated by an underscore (ˍ). It handles both general 
  # cases and special cases where the species name includes the string "f. sp.".
  #
  # Args:
  #   name: A string representing the full species name. The name can include 
  #         special forms denoted by "f. sp.".
  #
  # Returns:
  #   A string with the abbreviated species name.
  #
  # Examples:
  #   # General case
  #   abbreviate_species("Homo sapiens")
  #   # Returns "Hˍsapiens"
  #
  #   # Special case
  #   abbreviate_species("Fusarium oxysporum f. sp. lycopersici")
  #   # Returns "Fˍoxysporum"
  
  if (grepl("f\\. sp\\.", name)) {
    # Handle the special case
    sub("(^\\w)\\w+\\s(\\w+).*", "\\1ˍ\\2", name)
  } else {
    # Handle the general case
    sub("(^\\w)\\w+\\s(\\w+)", "\\1ˍ\\2", name)
  }
}

class_colors <- c("drug target" = "#8C0250", 
                  "multidrug transporter" = "#238CCD", 
                  "transcription factor" = "#C7B3FF", 
                  "other" = "#616161", 
                  "human&animal" = "#DC91BC", 
                  "human&plant" = "#B0E1FF", 
                  "plant" = "#61A184",
                  "saprotroph" = "white",
                  "SDH inhibitors" = "#DA8632",
                  "Qo inhibitors" = "#E3C398",
                  "Pyrimidine analogs" = "#959595",
                  "Polyenes" = "#71AED4",
                  "Echinocandins"= "#0F767A",
                  "Dicarboximides" = "#749603",
                  "Tubulin polymerization" = "#C1EB3A",
                  "Clinical azoles" = "#420DA5",
                  "Agricultural azoles" = "#420DA5",
                  "Squalene epoxidase" = "#106BA7")

#### Import files ####
setwd("/home/aliciapageau/Documents/antifungal_project/mardy2.0/clean_code")

# Read the CSV file
database <- read.csv("FungAMR.csv", stringsAsFactors = FALSE, row.names = 1, na.strings = 'NA')

drugs_class <- read.csv("drugs_class.csv", stringsAsFactors = FALSE) %>% 
  filter(usage != 'Experimental') %>% 
  select(drug,class) %>% 
  mutate(ID = paste0(class, ".", drug)) 

genes_class <- read.csv("genes_class.csv", stringsAsFactors = FALSE) %>%
  filter(class != 'Unknown') %>%
  mutate(ID = paste0(class, ".", ortho_homolog))


#### Prep data  ####
# Adjust gene and drug annotation
database <- database %>% 
  mutate(ortho_homolog = str_replace_all(ortho_homolog, "Cytochrome b", "CytB")) %>%
  mutate(ortho_homolog = str_replace_all(ortho_homolog, "Squalene epoxidase", "Sqle"))%>%
  mutate(ortho_homolog = str_replace_all(ortho_homolog, "Beta-tubulin", "β-tubulin"))%>%
  mutate(drug = str_replace_all(drug, "Amphotericin B", "AmphotericinB"))

# Filter for resistance mutations (positive degree of evidence)
# Filter for single gene mutations
single_gene <- database %>%
  filter(degree.of.evidence > 0) %>%
  filter(!grepl("\\|", mutation))

## Count for edgebundle plot
# Count by gene, specie and drug
# Filter for ortho_homolog present more than 5 times in the database
# Filter for species present more than 5 times in the database
gene_species_drug_counts <- single_gene %>%
  group_by(ortho_homolog) %>%
  filter(n() > 5) %>%
  group_by(species) %>%
  filter(n() > 5) %>%
  group_by(ortho_homolog, species, drug) %>%
  summarise(count = n(),
            evidence = min(best.degree.of.evidence))

# Filter out unclear drug and gene
# Filter out hybrids species
# Abbreviate species names
gene_species_drug_counts <- gene_species_drug_counts %>%
  filter(drug != 'No drug' & drug != 'Unknown' &
           !ortho_homolog %in% c("Hypothetical protein",
                                 "Predicted regulatory region of PbCyp51",
                                 "Putative fggy-family carbohydrate kinase",
                                 "Upstream of Cyp51", 'Pmr5 ABC transporter')) %>%
  filter(!species %in% c('Tspp.','Cryptococcus neoformans/Cryptococcus deneoformans',
                         'Saccharomyces cerevisiae/Saccharomyces paradoxus')) %>%
  mutate(species = sapply(species, abbreviate_species)) %>%
  select(-count,-evidence) %>%
  distinct() %>% ungroup()

#### Extended Fig 2 - Edgebundle plot genes, species, drugs ####
# Combine gene and species info into a new column
df <- gene_species_drug_counts %>% mutate(gene_species = paste0(ortho_homolog, "_", species))

# Create gene-species classification df and attribute IDs
genes_species <- df %>% select('ortho_homolog','species') %>% distinct()
genes_specie_class <- genes_class %>% left_join(genes_species) %>%
  mutate(gene_species = paste0(ortho_homolog, "_", species),
         ID = paste0(class, ".", gene_species)) %>%
  filter(!is.na(species))
genes_specie_class <- genes_specie_class[, c('gene_species','class','ID')]

# Convert genes and species into IDs with their respective classifications
edge <- df %>% 
  # Join with genes-species classification to replace genes-species names with IDs
  inner_join(genes_specie_class, by = 'gene_species') %>% 
  mutate(gene_species = ID) %>% 
  select(-ID, -class) %>% 
  
  # Join with species classification to replace species names with IDs
  inner_join(drugs_class, by = 'drug') %>% 
  mutate(drug = ID) %>% 
  select(-ID, -class)

# Select only relevant columns for the edge list
edge <- edge[, c('drug','gene_species')]

# Filter for genes-species that are associated with more than one drug
# This helps in focusing on cross-resistance patterns
edge <- edge %>% group_by(gene_species) %>% filter(n() > 1) 

## Make a grouping dataframe for the edgebundle graph (must contain all unique IDs present in the edges)
# Filter the groups to include only those IDs that are in the edges
groups <- full_join(genes_specie_class, drugs_class)
groups <- groups %>% filter((ID %in% unique(edge$gene_species)) | (ID %in% unique(edge$drug)))
groups <- groups[, c('ID', 'class')]

## Create the edgebundle plot
# Convert the edges and groups into a graph object
graph <- graph_from_data_frame(edge, directed = TRUE, vertices = groups)

# Set vertex color based on class using a predefined color palette
V(graph)$color <- class_colors[V(graph)$class]

# Set vertex size based on the degree (number of connections), scaled by 0.5
V(graph)$size <- degree(graph) * 0.5

# Generate the edgebundle plot with specified aesthetics
edgebundle(graph, tension = 0.6, fontsize = 20, width = 1200, padding = 260)


#### Legends ####

# Define the class colors dictionary
class_colors_gene <- c(
  "Drug target" = "#8C0250",
  "Efflux pumps" = "#238CCD",
  "Transcription factor" = "#C7B3FF",
  "Other" = "#616161"
)

class_colors_drug <- c(
  "Echinocandins"= "#0F767A",
  "Pyrimidine analogs" = "#959595",
  "Squalene epoxidase" = "#106BA7",
  "Polyenes" = "#71AED4",
  "Clinical azoles" = "#420DA5",
  "Agricultural azoles" = "#420DA5",
  "SDH inhibitors" = "#DA8632",
  "Qo inhibitors" = "#E3C398",
  "Dicarboximides" = "#749603",
  "Tubulin polymerization" = "#C1EB3A"
)


# Function to create custom legends with dots and title on the left in bold
create_dot_legend <- function(labels, colors, title) {
  plot.new()
  
  # Manually place the title on the left in bold
  text(x = 0, y = 1, labels = title, pos = 4, font = 2, cex = 1.2)  
  
  # Adjust legend to the right of the title
  legend(x = 0, y = 1, legend = labels, pch = 19, col = colors, bty = "n", pt.cex = 2, cex = 1.2)
}

# Genes class legend
svg(filename = "../figures/final/legend_genes.svg", width = 5, height = 5)
create_dot_legend(names(class_colors_gene), class_colors_gene, "Gene class")
dev.off()

# Drugs class legend
svg(filename = "../figures/final/legend_drugs.svg", width = 5, height = 5)
create_dot_legend(names(class_colors_drug), class_colors_drug, "Drug class")
dev.off()

