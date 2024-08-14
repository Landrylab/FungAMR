rm(list=ls())
dir()


library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggtree)
library(ape)
library(RRphylo)
library(castor)

sessionInfo()
#R version 4.2.2 (2022-10-31 ucrt)
#ggtree_3.7.2
#tidyr_1.3.0
#dplyr_1.1.4
#ape_5.6-2
#castor_1.8.0
#RRphylo_2.8.0






#################
###--- DATA---###
#################


### Reading genomic tree from Li et al 2021
tree1 <- read.tree("1672taxa_290genes_bb_1.treefile")
length(tree1$tip.label)

# Metadata with species names from the database
df <- read.csv("01_species_for_tree_with_meta.tsv", sep="\t", header=TRUE)
head(df)
dim(df)

# Selectd outgroup species from the Li et al 2021
dout <- read.csv("phylogeny_li_root.txt", sep="\t", header=FALSE)
outgroup <- dout$V1
outgroup

# Filtering big tree
species <- df$label
length(species)
tree1.sub <- keep.tip(tree1, species)
rtree1 <- root(tree1.sub, outgroup = outgroup, resolve.root = TRUE)

# Renaming some of the species
rtree1$tip.label[rtree1$tip.label == "Pyrenophora_tritici-repentis"] = "Pyrenophora_tritici_repentis"
rtree1$tip.label[rtree1$tip.label == "Ramularia_collo-cygni"] = "Ramularia_collo_cygni"
length(rtree1$tip.label)


# ITS tree
#tree2 <- read.tree("clean_combined.fasta.treefile")
#tree2$tip.label
#tree2 <- drop.tip(tree2, "u-Fusarium_solani")
#rtree2 <- root(tree2, outgroup = "its-Mucor_circinelloides", resolve.root = TRUE)
#rtree2



#######################
###--- FUNCTIONS ---###
#######################




getTreeSpecs <- function(tree) {
  
  # Tree nodes
  node.labs <- tree$node.label
  tips.labs <- tree$tip.label
  all.labs <- c(tips.labs, node.labs[c(2:length(node.labs))])
  
  # Table with edge lengths
  tree.edges <- data.frame(tree$edge)
  colnames(tree.edges)=c("parent", "node")
  
  # Calc distance from root
  root.dist <- get_all_distances_to_root(tree, as_edge_count=FALSE)
  root.dist.ordered = root.dist[tree.edges$node]
  all.labs.ordered <- all.labs[tree.edges$node]
  edge.length <- tree$edge.length
  edge.length.ordered <- edge.length[tree.edges$node]
  
  tree.edges$node_name <- all.labs.ordered
  tree.edges$root.dist <- root.dist.ordered
  tree.edges$edge.length <- edge.length.ordered
  tree.edges$num <- c(1:length(tree.edges$node))
  
  return(tree.edges)
}







getMergedTree <- function(focal, sister.id, out.id,
                          focal.name, sister.name, out.name,
                          tree_its, genomic_tree_merged) {
  
  p2 <- ggtree(tree_its) +
    geom_tiplab()
  
  if (length(sister.id) == 1) {
    sister <- sister.id
  } else if (length(sister.id) == 2) {
    sister <- MRCA(p2, sister.id)
  } 
  if (length(out.id) == 1) {
    out <- out.id
  } else if (length(out.id) == 2) {
    out <- MRCA(p2, out.id)
  }
  
  base <- MRCA(p2, sister, out)
  mrca_fs <-  MRCA(p2, focal, sister)
  
  ix <- get_pairwise_distances(tree_its, mrca_fs, focal)
  iy1 <- get_pairwise_distances(tree_its, base, mrca_fs)
  iy2 <- get_pairwise_distances(tree_its, sister, mrca_fs)
  iC <- get_pairwise_distances(tree_its, base, focal)
  iB <- get_pairwise_distances(tree_its, base, sister)
  
  ix.iB <- ix/iB
  iy1.iB <- iy1/iB
  iy2.iB <- iy2/iB
  
  # Genomic
  #genomic.tree <- merged.tree.2
  p1 <- ggtree(genomic_tree_merged) + geom_tiplab()
  
  if (length(sister.name) == 1) {
    gen.sister <- sister.name
  } else if (length(sister.name) == 2) {
    gen.sister <- MRCA(p1, sister.name)
  } 
  if (length(out.name) == 1) {
    gen.out <- out.name
  } else if (length(out.name) == 2) {
    gen.out <- MRCA(p1, out.name)
  }  
  
  gen.base <- MRCA(p1, gen.sister, gen.out)
  
  gB <- get_pairwise_distances(genomic_tree_merged, gen.base, gen.sister)
  gx <- gB * ix.iB
  
  ### Merging
  
  sister.names <- paste(sister.name, collapse="-")
  dato = data.frame("bind" = c(focal.name),
                    "reference" = c(sister.names),
                    "column" = c(FALSE))
  
  ### Getting node names, IDs, and edge lengths for species
  merged.tree <- tree.merger(backbone=genomic_tree_merged, data=dato)
  dmerg <- getTreeSpecs(merged.tree)
  
  ###
  pm <- ggtree(merged.tree) + geom_tiplab()
  
  focal_sister <- c(focal.name, sister.name)
  node.node <- MRCA(pm, focal_sister)
  if (length(sister.name) == 1) {
    sis.node <- sister.name
    sis.num <- dmerg %>% filter(node_name == sister.name) %>% pull(num)
  } else if (length(sister.name) == 2) {
    sis.node <- MRCA(pm, sister.name)
    sis.num <- dmerg %>% filter(node == sis.node) %>% pull(num)
  } 
  tip.node <- dmerg %>% filter(node_name == focal.name) %>% pull(node)
  out.node <- dmerg %>% filter(node == node.node) %>% pull(parent)
  
  # Node ordered numbers in the tree
  tip.num <- dmerg %>% filter(node == tip.node) %>% pull(num)
  node.num <- dmerg %>% filter(node == node.node) %>% pull(num)
  
  # New edge length
  y1 <- get_pairwise_distances(merged.tree, out.node, node.node)
  y2 <- get_pairwise_distances(merged.tree, node.node, sis.node)
  #y3 <- get_pairwise_distances(merged.tree, out.node, sis.node)
  tot.len <- y1 + y2
  #print(c(y1,y2,y3))
  
  # Transformed tree
  merged.tree.2 <- merged.tree
  merged.tree.2$edge.length[tip.num] <- gx
  merged.tree.2$edge.length[node.num] <- tot.len * iy1.iB
  merged.tree.2$edge.length[sis.num] <- tot.len * iy2.iB
  
  return(c(merged.tree,merged.tree.2))
  
}




################################################################################



### Checking tree edges from the genomic species tree
rtree1.edges <- getTreeSpecs(rtree1)
head(rtree1.edges)

### Table with branch lengths fro the genomic tree
x1 <- as_tibble(rtree1)
x1

# Plotting genomic tree for check
p0 <- ggtree(rtree1) +
  geom_tiplab() +
  xlim(0, 2.5)
p0


################################################################################

######################    Adding branches one-by-one

################################################################################





##############################
#      Adding  branches      #
##############################


#  |-------- outgroup
#  |
#--|
#  |   |--- focal
#  |---|
#      |----- sister
#
#       ---   x (distance mrca_fs - focal)
#       ----- y2 (distance mrca_fs - sister)
#   ---       y1 (distance base - mrca_fs)
#   -------   C (distance base - focal)
#   ----------B (distance base - sister)
#
# base    = mrcs(outgroup, focal, sister)
# mrca_fs = mrca(focal, sister) 



################################################################################

### candida
# ITS tree
tree2 <- read.tree("candida.fasta.treefile")
rtree2 <- root(tree2, outgroup = "its-Mucor_circinelloides", resolve.root = TRUE)
p2 <- ggtree(rtree2) +   geom_tiplab()
p2

### Candida_metapsilosis
focal <- "u-Candida_metapsilosis"
sister <- "its-Candida_orthopsilosis"
out <- "its-Candida_parapsilosis"

focal.name <- "Candida_metapsilosis"
sister.name <- "Candida_orthopsilosis"
out.name <- "Candida_parapsilosis"

new.tree <- getMergedTree(focal, sister, out,
                     focal.name, sister.name, out.name,
                     rtree2, rtree1)
new.tree.out <- new.tree[[2]]

ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)


################################################################################

### fusarium
# ITS tree
tree2 <- read.tree("fusarium.fasta.treefile")
rtree2 <- root(tree2, outgroup = "its-Aspergillus_fumigatus", resolve.root = TRUE)
p9 <- ggtree(rtree2) +   geom_tiplab() + xlim(0, 1)
p9

### Colletotrichum cereale

focal <- "ncbi-Colletotrichum_cereale"
sister <- c("its-Colletotrichum_nymphaeae","its-Colletotrichum_acutatum")
out <- c("its-Colletotrichum_gloeosporioides")

focal.name <- "Colletotrichum_cereale"
sister.name <- c("Colletotrichum_nymphaeae","Colletotrichum_acutatum")
out.name <- c("Colletotrichum_gloeosporioides")

new.tree <- getMergedTree(focal, sister, out,
                     focal.name, sister.name, out.name,
                     rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]


### Colletotrichum_asianum
p9
focal <- "its-Colletotrichum_asianum"
sister <- c("its-Colletotrichum_fructicola","its-Colletotrichum_gloeosporioides")
out <- "its-Colletotrichum_truncatum"

focal.name <- "Colletotrichum_asianum"
sister.name <- c("Colletotrichum_fructicola","Colletotrichum_gloeosporioides")
out.name <- "Colletotrichum_truncatum"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
new.tree.out <- new.tree[[2]]

ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)


### Colletotrichum_siamense
p9
focal <- "its-Colletotrichum_siamense"
sister <- "its-Colletotrichum_gloeosporioides"
out <- "its-Colletotrichum_fructicola"

focal.name <- "Colletotrichum_siamense"
sister.name <- "Colletotrichum_gloeosporioides"
out.name <- "Colletotrichum_fructicola"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
new.tree.out <- new.tree[[2]]

ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)




### Gibberella_fujikuroi
p9

# Fusarium_incarnatum
p9
focal <- "ncbi-Fusarium_incarnatum"
sister <- c("its-Fusarium_asiaticum","its-Fusarium_pseudograminearum")
out <- "ncbi-Gibberella_fujikuroi"

focal.name <- "Fusarium_incarnatum"
sister.name <- c("Fusarium_asiaticum","Fusarium_pseudograminearum")
out.name <- "Fusarium_fujikuroi"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]



### Diaporthe_eres
focal <- "its-Diaporthe_eres"
sister <- c("its-Beauveria_bassiana","its-Colletotrichum_acutatum")
out <- "its-Magnaporthe_grisea"

focal.name <- "Diaporthe_eres"
sister.name <- c("Beauveria_bassiana","Colletotrichum_acutatum")
out.name <- "Magnaporthe_grisea"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]



### Diaporthe_unshiuensis
p9
focal <- "u-Diaporthe_unshiuensis"
sister <- "its-Diaporthe_eres"
out <- "its-Magnaporthe_grisea"

focal.name <- "Diaporthe_unshiuensis"
sister.name <- "Diaporthe_eres"
out.name <- "Magnaporthe_grisea"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]


############################
### Tapesia_yallundae
p9
focal <- "ncbi-Tapesia_yallundae"
sister <- c("u-Uncinula_necator","nc-Blumeria_graminis23")
out <- "u-Botrytis_cinerea"

focal.name <- "Tapesia_yallundae"
sister.name <- c("Erysiphe_necator","Blumeria_graminis")
out.name <- "Botrytis_cinerea"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]



### Tapesia_acuformis
p9
focal <- "ncbi-Tapesia_acuformis"
sister <- c("ncbi-Tapesia_yallundae")
out <- c("u-Uncinula_necator","nc-Blumeria_graminis23")

focal.name <- "Tapesia_acuformis"
sister.name <- c("Tapesia_yallundae")
out.name <- c("Erysiphe_necator","Blumeria_graminis")

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]



### Pyrenopeziza_brassicae
p9
focal <- "ncbi-Pyrenopeziza_brassicacae"
sister <- c("ncbi-Tapesia_yallundae","ncbi-Tapesia_acuformis")
out <- c("u-Uncinula_necator","nc-Blumeria_graminis23")

focal.name <- "Pyrenopeziza_brassicae"
sister.name <- c("Tapesia_yallundae","Tapesia_acuformis")
out.name <- c("Erysiphe_necator","Blumeria_graminis")

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]



################################################################################

### aspergillus
# ITS tree
tree2 <- read.tree("aspergillus.fasta.treefile")
rtree2 <- root(tree2, outgroup = "its-Mucor_circinelloides", resolve.root = TRUE)
rtree2
p8 <- ggtree(rtree2) +   geom_tiplab() + xlim(0, 1)
p8


### Trichophyton_indotineae
focal <- "its-Trichophyton_indotineae"
sister <- c("its-Trichophyton_interdigitale","ncbi-Trichophyton_mentagrophytes")
out <- "its-Arthroderma_vanbreuseghemii"

focal.name <- "Trichophyton_indotineae"
sister.name <- c("Trichophyton_interdigitale","Trichophyton_mentagrophytes")
out.name <- c("Trichophyton_tonsurans")

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]




################################################################################

### cryptococcus
# ITS tree
tree2 <- read.tree("cryptococcus.fasta.treefile")
rtree2 <- root(tree2, outgroup = "its-Mucor_circinelloides", resolve.root = TRUE)
rtree2
p7 <- ggtree(rtree2) +   geom_tiplab() + xlim(0, 2)
p7


### Cryptococcus_deneoformans
focal <- "ncbi-Cryptococcus_deneoformans"
sister <- "its-Cryptococcus_neoformans"
out <- "its-Cryptococcus_gattii"

focal.name <- "Cryptococcus_deneoformans"
sister.name <- "Cryptococcus_neoformans"
out.name <- "Cryptococcus_gattii_VGIV"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)

ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]


### Rhizoctonia_cerealis
p7
focal <- "ncbi-Rhizoctonia_cerealis"
sister <- c("ncbi-Rhizoctonia_solani")
out <- c("ncbi-Cryptococcus_deneoformans", "its-Cryptococcus_gattii")

focal.name <- "Rhizoctonia_cerealis"
sister.name <- c("Rhizoctonia_solani")
out.name <- c("Cryptococcus_deneoformans","Cryptococcus_gattii_VGIV")

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]


### Phakopsora_pachyrhizi
p7
focal <- "u-Phakopsora_pachyrhizi"
sister <- c("ncbi-Puccinia_triticina","u-Puccinia_horiana")
out <- c("ncbi-Cryptococcus_deneoformans","ncbi-Rhizoctonia_cerealis")

focal.name <- "Phakopsora_pachyrhizi"
sister.name <- c("Puccinia_triticina","Puccinia_horiana")
out.name <- c("Cryptococcus_deneoformans","Rhizoctonia_cerealis")

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]




################################################################################

### zymoseptoria
# ITS tree
tree2 <- read.tree("zymoseptoria.fasta.treefile")
rtree2 <- root(tree2, outgroup = "its-Candida_albicans", resolve.root = TRUE)
rtree2
p6 <- ggtree(rtree2) +   geom_tiplab() + xlim(0, 2)
p6

### Septoria_glycines
p6
focal <- "ncbi-Septoria_glycines"
sister <- c("its-Cercospora_sojina")
out <- "its-Cercospora_beticola"

focal.name <- "Septoria_glycines"
sister.name <- c("Cercospora_sojina")
out.name <- "Cercospora_beticola"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]


### Venturia_nashicola
p6
focal <- "ncbi-Venturia_nashicola"
sister <- "its-Venturia_inaequalis"
out <- "its-Passalora_fulva"

focal.name <- "Venturia_nashicola"
sister.name <- "Venturia_inaequalis"
out.name <- "Passalora_fulva"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]



### Phaeosphaeria_nodorum
p6
focal <- "ncbi-Phaeosphaeria_nodorum"
sister <- c("u-Alternaria_alternata","u-Pyrenophora_tritici-repentis")
out <- "nc-Corynespora_cassiicola14"

focal.name <- "Phaeosphaeria_nodorum"
sister.name <- c("Alternaria_alternata","Pyrenophora_tritici_repentis")
out.name <- "Corynespora_cassiicola"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]


### Didymella_bryoniae
p6
focal <- "ncbi-Didymella_bryoniae"
sister <- "ncbi-Phaeosphaeria_nodorum"
out <- "u-Alternaria_alternata"

focal.name <- "Didymella_bryoniae"
sister.name <- "Phaeosphaeria_nodorum"
out.name <- "Alternaria_alternata"

new.tree <- getMergedTree(focal, sister, out,
                          focal.name, sister.name, out.name,
                          rtree2, new.tree.out)
ggtree(new.tree[[1]]) +
  geom_tiplab() +
  xlim(0, 3)

ggtree(new.tree[[2]]) +
  geom_tiplab() +
  xlim(0, 3)
new.tree.out <- new.tree[[2]]



################################################################################






##############################
#   WRITING DOWN THE TREE    #
##############################



write.tree(new.tree.out, "adding_branches.tree")



plot(new.tree.out)
nodelabels(bg="white")
tiplabels(bg="white")


