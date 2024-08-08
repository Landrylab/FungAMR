### Created on 31-07-2024
### @author: aliciapageau

library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(edgebundleR)
library(ggplot2)


#### Load functions ####
custom_combinations <- function(drugs) {
  # Custom Combinations of Drugs
  #
  # This function generates all possible pairs of combinations from a given vector 
  # of drugs. The pairs are returned as a data frame with two columns: "from" and "to".
  # If the length of the input vector is less than 2, the function returns NULL.
  #
  # Args:
  #   drugs: A vector of strings, where each string represents a drug.
  #
  # Returns:
  #   A data frame with two columns, "from" and "to", representing all possible 
  #   pairs of drugs. If the input vector has less than two elements, the function 
  #   returns NULL.
  #
  # Examples:
  #   custom_combinations(c("Aspirin", "Ibuprofen", "Paracetamol"))
  #   # Returns a data frame with pairs:
  #   #     from         to
  #   # 1 Aspirin  Ibuprofen
  #   # 2 Aspirin Paracetamol
  #   # 3 Ibuprofen Paracetamol
  #
  #   custom_combinations(c("Aspirin"))
  #   # Returns NULL
  
  l <- length(drugs)
  if (l > 1) {
    comb <- combn(drugs, 2)
    comb <- as.data.frame(t(comb))
    colnames(comb) <- c("from", "to")
    return(comb)
  } else {
    return(NULL)
  }
}

class_colors <- c("drug target" = "#8C0250", 
                  "multidrug transporter" = "#238CCD", 
                  "transcription factor" = "#C7B3FF", 
                  "other" = "#616161", 
                  ####
                  "human&animal" = "#DC91BC", 
                  "human&plant" = "#B0E1FF", 
                  "plant" = "#61A184",
                  "saprotroph" = "white",
                  ####
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

# Filter for single mutation
single_mut <- single_gene %>%
  filter(!grepl(",", mutation))

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
gene_species_drug_counts <- gene_species_drug_counts %>%
  filter(drug != 'No drug' & drug != 'Unknown' &
           !ortho_homolog %in% c("Hypothetical protein",
                                 "Predicted regulatory region of PbCyp51",
                                 "Putative fggy-family carbohydrate kinase",
                                 "Upstream of Cyp51", 'Pmr5 ABC transporter')) %>%
  filter(!species %in% c('Tspp.','Cryptococcus neoformans/Cryptococcus deneoformans',
                         'Saccharomyces cerevisiae/Saccharomyces paradoxus')) %>%
  select(-count,-evidence) %>%
  distinct() %>% ungroup()

# Count by mutation, gene, and drug
# Filter for ortho_homolog present more than 5 times in the database
mut_counts <- single_mut %>%
  group_by(ortho_homolog, ortho_mut, drug) %>%
  summarise(count = n(),
            evidence = min(best.degree.of.evidence),
            mutation = first(mutation)) 

#### Figure 5A - Edgebundle plot drug to drug ####
## Create from-to table establishing unique mutation cross-resistance between drugs
# Filter out rows with NA orthologous mutations and group by ortho_mut
from_to <- mut_counts %>%
  filter(!is.na(ortho_mut)) %>%
  group_by(ortho_mut) %>%
  # For each unique mutation, create a list of drug combinations where cross-resistance is observed
  summarize(drug_combinations = list(custom_combinations(drug))) %>%
  unnest(drug_combinations) %>%
  ungroup() #%>%

# Remove the ortho_mut column as it is no longer needed
#select(-ortho_mut)

# Filter for drugs present in the drugs_class dataframe
# Ensures that only classified drugs are included in the analysis
from_to <- from_to %>% 
  filter((from %in% unique(drugs_class$drug)) & (to %in% unique(drugs_class$drug)))

## Convert 'from' and 'to' columns into drug IDs using the drugs_class dataframe
# Match drug names with their corresponding IDs and replace the names with IDs
from <- from_to %>% select(from) %>% rename(drug = from) %>%
  left_join(drugs_class) %>% select(ID) %>% rename(from = ID)
to <- from_to %>% select(to) %>% rename(drug = to) %>%
  left_join(drugs_class) %>% select(ID) %>% rename(to = ID)

# Combine 'from' and 'to' columns into a single dataframe for edge representation
edges <- bind_cols(from, to)

## Create a grouping dataframe for the edgebundle graph
# Ensure that the groups dataframe contains only the IDs present in the edges
groups <- drugs_class %>% filter((ID %in% unique(edges$from)) | (ID %in% unique(edges$to)))
groups <- groups[, c("ID", "class", "drug")]

## Create the edgebundle plot
# Convert the edges and groups into a graph object
graph <- graph_from_data_frame(edges, directed = TRUE, vertices = groups)

# Set vertex color based on the class using a predefined color palette
V(graph)$color <- class_colors[V(graph)$class]

# Set vertex size based on the degree (number of connections), scaled by 0.5
V(graph)$size <- degree(graph) * 0.5

# Generate the edgebundle plot with specified aesthetics
edgebundle(graph, tension = 0.6, fontsize = 20, width = 1000, padding = 150)


#### Figure 5B - Edgebundle plot genes and drugs ####
## Create edges from genes to species in the database
# Select distinct combinations of genes and drugs, excluding the species column
df <- gene_species_drug_counts %>% 
  select(-species) %>% distinct()

# Convert genes and species into IDs with their respective classifications
edges <- df %>% 
  # Join with gene classification to replace gene names with IDs
  inner_join(genes_class, by = 'ortho_homolog') %>% 
  mutate(ortho_homolog = ID) %>% 
  select(-ID, -class) %>% 
  
  # Join with species classification to replace species names with IDs
  inner_join(drugs_class, by = 'drug') %>% 
  mutate(drug = ID) %>% 
  select(-ID, -class)

# Select only relevant columns for the edge list
edges <- edges[, c("drug", "ortho_homolog")]

# Filter for genes that are associated with more than one species
# This helps in focusing on convergence patterns
edges <- edges %>% 
  group_by(ortho_homolog) %>% 
  filter(n() > 1) 

## Make a grouping dataframe for the edgebundle graph (must contain all unique IDs present in the edges)
# Filter the groups to include only those IDs that are in the edges
groups <- full_join(genes_class, drugs_class)
groups <- groups %>% filter((ID %in% unique(edges$ortho_homolog)) | (ID %in% unique(edges$drug)))
groups <- groups[, c("ID", "class")]

## Create the edgebundle plot
# Convert the edges and groups into a graph object
graph <- graph_from_data_frame(edges, directed = TRUE, vertices = groups)

# Set vertex color based on class using a predefined color palette
V(graph)$color <- class_colors[V(graph)$class]

# Set vertex size based on the degree (number of connections), scaled by 0.5
V(graph)$size <- degree(graph) * 0.5

# Generate the edgebundle plot with specified aesthetics
edgebundle(graph, tension = 0.6, fontsize = 20, width = 1000, padding = 150)


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
