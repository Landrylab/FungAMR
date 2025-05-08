## ----------------------
##
## Author: Anna Fijarczyk
##
## Date: 23-02-2024
##
## ----------------------


rm(list=ls())
#setwd("C:/Users/aniaf/Projects/GAPP/21_mardy_paper_phylogeny/plot_phylogeny_clean_scripts")
setwd("C:/Users/aniaf/Projects/GAPP/22_mardy_paper_phylogeny")
dir()


library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggtree)
library(ape)
library(ggnewscale)


sessionInfo()
#R version 4.4.2 (2024-10-31 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 26100)

#other attached packages:
#  [1] adephylo_1.1-16    ade4_1.7-22        castor_1.8.3       Rcpp_1.0.13-1      RRphylo_2.8.1      emmeans_1.10.7    
#[7] ggnewscale_0.5.0   readxl_1.4.3       ggtree_3.15.0      RColorBrewer_1.1-3 phytools_2.4-4     maps_3.4.2.1      
#[13] ape_5.8            cowplot_1.1.3      tidyr_1.3.1        dplyr_1.1.4        ggplot2_3.5.1     



####################
###     DATA     ###
####################


dh_meta <- read.csv("metatable_for_phylogeny.tsv", sep="\t", header = TRUE)
head(dh_meta)
dim(dh_meta)

dsub_meta <- dh_meta %>% filter(!is.na(Subphylum))
head(dsub_meta)
dim(dsub_meta)
unique(dsub_meta$Subphylum)



# Reading tree
tree <- read.tree("../21_mardy_paper_phylogeny/plot_phylogeny_clean_scripts/adding_branches.tree")
tree
tree_labeles <- tree$tip.label
tree_labeles
length(tree_labeles)



outgroup <- c("Capniomyces_stellatus", "Smittium_mucronatum", "Conidiobolus_coronatus")
rooted_tree.sub <- drop.tip(tree, outgroup)
rooted_tree.sub


length(rooted_tree.sub$tip.label)
dsub_meta[!dsub_meta$label %in% rooted_tree.sub$tip.label,]


####################
###- PLOT TREE  -###
####################


########################
### Version downward ###
########################


### Tree with colored tips

dsub_meta$Subphylum <- factor(dsub_meta$Subphylum, levels = c("Pezizomycotina","Saccharomycotina", "Taphrinomycotina",
                                                              "Agaricomycotina","Pucciniomycotina","Ustilaginomycotina",
                                                              "Mucoromycotina"))

p2 <- ggtree(rooted_tree.sub) %<+% dsub_meta + 
  geom_tiplab(align=TRUE, aes(label = species_name), angle=270, fontface="italic") +
  geom_tippoint(aes(fill=Subphylum),
                pch=21,
                size=4,
                color="black",
                alpha=0.95) +
  scale_fill_manual(values=c("#DAACC6","#91CEF4","#C7B3FF","#83BBBD","#71AED4","#A67490","#C8C8C8"),
                    name="Subphylum",
                    guide = guide_legend(direction="vertical")) +
  new_scale_fill() +
  theme(legend.position = "bottom") +
  layout_dendrogram()
p2




### Adding heatmap with host information

df_mat <- dsub_meta %>% dplyr::select(Host)
rownames(df_mat) <- dsub_meta$label
head(df_mat)

p2H <- gheatmap(p2, df_mat, offset=1, width=0.1, colnames=FALSE) +
  scale_fill_manual(values = c("#DAACC6","#91CEF4","#83BBBD","white"),
                    name="Host", na.translate = F,
                    guide = guide_legend(direction="vertical")) +
  new_scale_fill() +
  theme(legend.position = "bottom")
p2H


### Adding heatmap with mutation information

df2_mat <- dsub_meta %>% dplyr::select(log10_mutations)
rownames(df2_mat) <- dsub_meta$label
head(df2_mat)

df2_mat$group <- as.numeric(cut(df2_mat$log10_mutations, breaks = c(-1,1,2,3,4)))
df2_mat$group
head(df2_mat)

df3_mat <- df2_mat %>% dplyr::select("group")
df3_mat$group <- as.factor(df3_mat$group)

p2M <- gheatmap(p2H, df3_mat, offset=1.2, width=0.1, colnames=FALSE) +
  #scale_fill_gradient(low="white", high="grey20", name="Log10 number\nof mutations") +
  scale_fill_manual(values=c("grey95","grey80","grey60","grey20"),
                    labels=c("1-10","10-100","100-1000","~4000"),
                    name="Number\nof mutations",
                    guide = guide_legend(direction="vertical")) +
  theme(legend.position = "bottom")
p2M







ggsave(file="plotting_tree.svg", plot=p2M, w=3000,h=3000, units="px", dpi=300)

