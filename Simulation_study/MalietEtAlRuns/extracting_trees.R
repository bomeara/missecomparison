# rm(list=ls())

# setwd("/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison")
# Extracting trees from simulated data (ClaDS2) to run other methods
library(ape)
tree_files <- list.files("original_trees/trees/ClaDS2/")
tree_files <- tree_files[grep("Rdata",tree_files)]
for(i in 1:length(tree_files)) {
  load(paste0(getwd(),"/original_trees/trees/ClaDS2/", tree_files[i]))
  if(exists("tree")){
    label <- gsub(".Rdata","",tree_files[i])
    write.tree(tree, file=paste0(getwd(),"/original_trees/trees/ClaDS2/",label,".tre"))
    rm("tree")
  }
}
