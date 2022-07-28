#################################################################################################################################
####### PLOTTING MIC AND CLUSTER ON THE PHYLOGENETIC TREE (CLADE A) #############################################################
#################################################################################################################################

## This script uses ggtreeExtra to plot metadata on the phylogenetic tree file

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(treeio)
library(ggtree)
library(ggplot2)
library(gdata)
library(ggtreeExtra)
library(dplyr)
library(ggstar)
library(ggnewscale)

tr_dir <- "C://Users/rache/OneDrive/Documents/Cambridge/MPHIL/proj/"
tr_pre <- "efm_dap_dataset_v1.cladeA.rmRCB.rmMGE.iqtree.contree.tree"
trfile <- paste(tr_dir, tr_pre, sep="")
tree <- read.tree(trfile)

p <- ggtree(tree, layout="fan", open.angle=15, size = 0.1) # size = 0.1 chosen to make lines narower


ring1 <- read.csv("pheno_linear_A.csv")
tips <- complete_phylo
View(tips)#phylo cluster
tips$cluster <- tips$Cluster.DSF[[1]] 
View(tips)
tips <- tips[,-2]
View(tips)
tips$cluster <- as.factor(tips$cluster)
tips_a <- tips[-640:-669,]
View(tips_a)
library(RColorBrewer)
colourCount <- 46
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)

p1 <- p %<+% tips_a +
  geom_tippoint(aes(fill = cluster, colour=cluster),
                size = 0.1,
                shape = 1) +
  scale_colour_manual(values = getPalette,
                    name = "Phylogenetic cluster",
                    breaks = c(1:45, "singleton"),
                    guide=guide_legend(keywidth=0.2,
                                       keyheight = 0.2))

p2 <- p1 +
  new_scale_fill()+
  geom_fruit(
    data = ring1,
    geom=geom_col,
    mapping=aes(y=sample_id, x=trans_MIC),
    pwidth=0.4,
    orientation="y",
    axis.params = list(
      axis="trans_MIC",
      text.angle=-45,
      hjust=0
    ),
    grid.params = list()
  ) +
  theme(
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=3),
    legend.text=element_text(size=2),
    legend.spacing = unit(0, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    legend.spacing.x = unit(0, "cm"),
    legend.key.size = unit(0.02, "cm"),
    plot.margin = margin(10, 10, 10, 100)
  )
p2