library(ape)
library(phytools)
setwd("~/Desktop/misse_mme_paper/missecomparison/empirical_2021")
tree <- read.tree("Eucalypts_ML1_dated_r8s.tre")

tree<- drop.tip(tree, c("Kjellbergiodendron_celebicum",                     
                                    "Lophostemon_confertus", 
                                    "Lophostemon_suaveolens_CANB_750277",
                                    "Myrtus_communis",                     
                                    "Decaspermum_humile",
                                    "Amomyrtus_luma" , 
                                    "Cloezia_floribunda", 
                                    "Tepualia_stipularis",
                                    "Acmena_smithii",
                                    "Syzygium_angophoroides",
                                    "Syncarpia_glomulifera_Jervis_Bay",
                                    "Syncarpia_hillii_612122",
                                    "Hypocalymma_linifolium",
                                    "Agonis_flexuosa", 
                                    "Kunzea_ericoides"))



options(digits=4)
tree_modified <- tree
tree_modified <- multi2di(tree_modified)
tree_modified <- force.ultrametric(tree_modified)
tree_modified$edge.length <- tree_modified$edge.length + 0.0001
tree_modified <- ladderize(tree_modified)

# solution copied from https://jonathanchang.org/blog/three-ways-to-check-and-fix-ultrametric-phylogenies/
tre <- tree_modified
tre_node_adjust <- reorder(tre, "postorder")

e1 <- tre_node_adjust$edge[, 1] # parent node
e2 <- tre_node_adjust$edge[, 2] # child node
EL <- tre_node_adjust$edge.length
ages <- numeric(N + tre_node_adjust$Nnode)

for (ii in seq_along(EL)) {
  if (ages[e1[ii]] == 0) {
    ages[e1[ii]] <- ages[e2[ii]] + EL[ii]
  } else {
    recorded_age <- ages[e1[ii]]
    new_age <- ages[e2[ii]] + EL[ii]
    if (recorded_age != new_age) {
      cat(sprintf("node %i age %.6f != %.6f\n", e1[ii], recorded_age, new_age))
   EL[ii] <- recorded_age - ages[e2[ii]]
    }
  }
}

tre_node_adjust$edge.length <- EL

is.ultrametric(tre_node_adjust) # good
min(tre_node_adjust$edge.length) # good
write.tree(tre_node_adjust, "ML1_modified.tre")

#-------
install.packages("BAMMtools")
library(BAMMtools)
setBAMMpriors(read.tree("ML1_modified.tre"))

tree <- ladderize(force.ultrametric(read.tree(file.choose())))
par(mfrow=c(1,2))
plot(tree, show.tip.label = F)
plot(read.tree("ML1_modified.tre"), show.tip.label = F)
axisPhylo()
