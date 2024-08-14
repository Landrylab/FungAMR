rm(list=ls())
dir()


library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggtree)
library(ape)
library(ggnewscale)

sessionInfo()
#R version 4.2.2 (2022-10-31 ucrt)
#tidyr_1.3.0
#dplyr_1.1.4
#ape_5.6-2
#ggtree_3.7.2
#ggnewscale_0.4.10
#RColorBrewer_1.1-3





####################
###     DATA     ###
####################


dh_meta <- read.csv("metatable_for_phylogeny.tsv", sep="\t", header = TRUE)
head(dh_meta)
dim(dh_meta)

dsub_meta <- dh_meta %>% filter(!is.na(Subphylum))
head(dsub_meta)
unique(dsub_meta$Subphylum)



# Reading tree
tree <- read.tree("adding_branches.tree")
tree
tree_labeles <- tree$tip.label
tree_labeles
length(tree_labeles)

outgroup <- c("Capniomyces_stellatus", "Smittium_mucronatum", "Conidiobolus_coronatus")
rooted_tree.sub <- drop.tip(tree, outgroup)
rooted_tree.sub






####################
###- PLOT TREE  -###
####################


########################
### Version downward ###
########################


### Tree with colored tips

p2 <- ggtree(rooted_tree.sub) %<+% dsub_meta + 
  geom_tiplab(align=TRUE, aes(label = Species), angle=270) +
  geom_tippoint(aes(fill=Subphylum),
                pch=21,
                size=4,
                color="black",
                #angle=90,
                alpha=0.95) +
  scale_fill_manual(values=c("#DAACC6","#91CEF4","#C7B3FF","#83BBBD","#71AED4","#C8C8C8"), name="Subphylum") +
  new_scale_fill() +
  theme(legend.position = "bottom") +
  layout_dendrogram()
p2


### Adding heatmap with host information

df_mat <- dsub_meta %>% dplyr::select(Host)
rownames(df_mat) <- dsub_meta$label
head(df_mat)

p2H <- gheatmap(p2, df_mat, offset=0.8, width=0.1, colnames=FALSE) +
  scale_fill_manual(values = c("#DAACC6","#91CEF4","#83BBBD","white"),
                    name="Host", na.translate = F) +
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

p2M <- gheatmap(p2H, df3_mat, offset=1, width=0.1, colnames=FALSE) +
  #scale_fill_gradient(low="white", high="grey20", name="Log10 number\nof mutations") +
  scale_fill_manual(values=c("grey95","grey80","grey60","grey20"), labels=c("1-10","10-100","100-1000","~4000"), name="Number\nof mutations") +
  theme(legend.position = "bottom")
p2M





#ggsave(file="plotting_tree.svg", plot=p2M, w=1600,h=1800, units="px", dpi=150)

