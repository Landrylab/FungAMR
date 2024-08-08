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

convert_pos <- read.csv("database_convert_pos_convergence.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1) %>%
  select(first.author.name, species, gene.or.protein.name, mutation, degree.of.evidence,
         ortho_homolog, wt_AA, position, alt_AA, alignment_pos, ortho_mut, ortho_res, 
         convert_pos, convert_wt, conversion_species) %>% distinct()

host <- read.csv("species_class.csv", stringsAsFactors = FALSE) %>% rename(host = class)
species_class <- read.csv("species_class.csv", stringsAsFactors = FALSE) %>%
  mutate(species = sapply(species, abbreviate_species),
         ID = paste0(class, ".", species)) %>% distinct

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

convert_pos <- convert_pos %>% 
  mutate(ortho_homolog = str_replace_all(ortho_homolog, "Cytochrome b", "CytB")) %>%
  mutate(ortho_homolog = str_replace_all(ortho_homolog, "Squalene epoxidase", "Sqle"))%>%
  mutate(ortho_homolog = str_replace_all(ortho_homolog, "Beta-tubulin", "β-tubulin"))

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

## Convergence bar plot
convergence <- single_gene %>% left_join(convert_pos) %>% 
  left_join(select(genes_class, -ID)) %>% 
  left_join(host)

#### Figure 4A - Edgebundle plot genes and species ####
## Create edges from genes to species in the database
# Select distinct combinations of genes and species, excluding the drug column
df <- gene_species_drug_counts %>% 
  select(-drug) %>% distinct()

# Convert genes and species into IDs with their respective classifications
edges <- df %>% 
  # Join with gene classification to replace gene names with IDs
  inner_join(genes_class, by = 'ortho_homolog') %>% 
  mutate(ortho_homolog = ID) %>% 
  select(-ID, -class) %>% 
  
  # Join with species classification to replace species names with IDs
  inner_join(species_class, by = 'species') %>% 
  mutate(species = ID) %>% 
  select(-ID, -class)

# Select only relevant columns for the edge list
edges <- edges[, c("species", "ortho_homolog")]

# Filter for genes that are associated with more than one species
# This helps in focusing on convergence patterns
edges <- edges %>% 
  group_by(ortho_homolog) %>% 
  filter(n() > 1) 

## Make a grouping dataframe for the edgebundle graph (must contain all unique IDs present in the edges)
# Filter the groups to include only those IDs that are in the edges
groups <- full_join(genes_class, species_class)
groups <- groups %>% filter((ID %in% unique(edges$ortho_homolog)) | (ID %in% unique(edges$species)))
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

#### Figure 4B - Bar plot convergence genes ####
# Define genes to exclude from the plot
excluded_gene <- c('Hmg2', 'Hap1', 'Myosin-1')

# Filter and process convergence data
convergence_gene <- convergence %>%
  group_by(ortho_homolog) %>%                  # Group by each ortho_homolog gene
  filter(n() > 5) %>%                          # Keep only genes with more than 5 entries
  mutate(
    n_species = n_distinct(species),           # Count distinct species for each ortho_homolog
    count = 1                                  # Set count to 1 for each row
  ) %>%
  select(ortho_homolog, species, n_species, count, host, class) %>%  # Select relevant columns
  distinct() %>%  ungroup() %>%                # Remove duplicate rows
  filter(
    n_species >= 2,                            # Keep genes associated with at least 2 species
    !ortho_homolog %in% excluded_gene          # Exclude specified genes that are not in edgebundle
  ) %>%
  arrange(n_species) %>%                       # Arrange by number of species
  ungroup() %>%
  mutate(ortho_homolog = factor(ortho_homolog, levels = unique(ortho_homolog))) # Convert ortho_homolog to factor

# Plot
custom_labels <- c('drug target' = "Drug \ntarget ", 'multidrug transporter' = "Efflux \npump",
                   'other' = "Other", 'transcription factor' = "Transcription \nfactor")
ggplot(convergence_gene, aes(x = count, y = ortho_homolog, fill = host)) + 
  geom_col() +
  xlab("Number of species") +
  ylab("Gene") +
  facet_grid(rows = vars(class), scales = "free_y", space = "free",
             labeller = labeller(class = as_labeller(custom_labels))) +
  theme(strip.text.y = element_text(angle = 0, size = 10),
        strip.placement = "outside",
        strip.background = element_rect(fill="gray92", color = "white", size = 0.5),
        axis.text.x = element_text(),
        legend.position = 'none') +
  scale_fill_manual(values=class_colors)
#ggsave("../figures/final/nb_species_genes_v4.png", device = "png", dpi = 600)
#ggsave("../figures/final/nb_species_genes_v4.svg", device = "svg", dpi = 600)

#### Figure 4C - Bar plot convergence mutations ####
# Filter and process convergence data
convergence_mut <- convergence %>%
  filter(!is.na(ortho_res)) %>%          # Remove rows with NA orthologous residue
  group_by(ortho_res) %>%                # Group by orthologous residue
  mutate(
    n_species = n_distinct(species),     # Count unique species
    count = 1,                           # Set count to 1 for each row
    mutation = paste(
      first(na.omit(convert_wt)),        # Get first non-NA convert_wt
      first(na.omit(convert_pos)),       # Get first non-NA convert_pos
      sep = ""
    ),
    numbering_sp = first(na.omit(conversion_species)), # Get first non-NA conversion_species
    convert_pos = first(na.omit(convert_pos))          # Get first non-NA convert_pos
  ) %>%
  select(ortho_homolog, mutation, numbering_sp, species, n_species, count, host, convert_pos) %>%
  distinct() %>% ungroup() %>%           # Remove duplicate rows and ungroup the data
  filter(n_species > 2) %>%              # Keep mutations present in more than 2 species
  arrange(desc(convert_pos), .by_group = TRUE) %>%  # Arrange by convert_pos within each group
  ungroup() %>%                          # Ungroup the data
  mutate(mutation = factor(mutation, levels = unique(mutation)))  # Convert mutation to a factor

# Plot
ggplot(convergence_mut, aes(x = count, y = mutation, fill = host)) + 
  geom_col() +
  xlab("Number of species") +
  ylab("Residu mutated") +
  facet_grid(rows = vars(ortho_homolog), scales = "free_y", space = "free") +
  theme(strip.text.y = element_text(angle = 0, size = 10),
        strip.placement = "outside",
        strip.background = element_rect(fill="gray92", color = "white", size = 0.5),
        axis.text.x = element_text(),
        legend.position = 'none') +
  scale_fill_manual(values=class_colors)
#ggsave("../figures/final/nb_species_mutations.png", device = "png", dpi = 600)
#ggsave("../figures/final/nb_species_mutations.svg", device = "svg", dpi = 600)

#### Legends ####

# Define the class colors dictionary
class_colors_gene <- c(
  "Drug target" = "#8C0250",
  "Efflux pumps" = "#238CCD",
  "Transcription factor" = "#C7B3FF",
  "Other" = "#616161"
)

class_colors_species <- c(
  "Human&animal" = "#DC91BC",
  "Human&plant" = "#B0E1FF",
  "Plant" = "#61A184",
  "Saprotroph" = "white"
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

# Species class legend
svg(filename = "../figures/final/legend_species.svg", width = 5, height = 5)
create_dot_legend(names(class_colors_species), class_colors_species, "Specie class")
dev.off()

