#Script inspired from https://davidcerny.github.io/post/plotting_beast/

library(phytools)
library(strap)
library(phyloch)

annot_tree <- phyloch::read.beast("/Users/yacinebenchehida/Desktop/Convergent_evolution/All_Melinaea/TreePL/New/3_calibrations/Results/treeannotator_results.txt")

if (is.null(annot_tree$`CAheight_95%_HPD_MIN`)) {
  annot_tree$min_ages <- annot_tree$`height_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`height_95%_HPD_MAX`
} else {
  annot_tree$min_ages <- annot_tree$`CAheight_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`CAheight_95%_HPD_MAX`
}

annot_tree$root.time <- max(nodeHeights(annot_tree)) 


pdf("busco_tree_treepl_3_calibrations.pdf",18,10)
geoscalePhylo(ladderize(annot_tree, right = F), x.lim = c(0, 170), cex.tip = 0.6, cex.age = 1, cex.ts = 0.9, units = c("Epoch","Age"),quat.rm=TRUE,width= 2, align.tip.label=FALSE, ftype = "b")

T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)

for(i in (Ntip(annot_tree) + 1):(annot_tree$Nnode + Ntip(annot_tree))) {
  lines(x = c(T1$root.time - annot_tree$min_ages[i - Ntip(annot_tree)],
              T1$root.time - annot_tree$max_ages[i - Ntip(annot_tree)]),
        y = rep(T1$yy[i], 2), lwd = 6, lend = 0,
        col = make.transparent("blue", 0.4))
}
dev.off()
