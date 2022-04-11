
setwd("~/Desktop/misse_mme_paper/missecomparison/Simulation_study/")
local_path <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/Simulation_study/" #local
lab_computer_path <- "/home/tvasconcelos/missecomparison/Simulation_study" #labcomputer

all_trees <- list.files(paste0(local_path, "MalietEtAlRuns/MalietEtAl_ClaDS2_trees/"))
all_trees <- all_trees[grep(".tre$", all_trees)]

# install.packages("BAMMtools")
library(BAMMtools)
library(phytools)


all_tree_clads <- list()
for(i in 1:length(all_trees)) {
  one_tree <- read.tree(paste0(local_path,"MalietEtAlRuns/MalietEtAl_ClaDS2_trees/", all_trees[i]))
  label <- gsub(".tre$","", all_trees[i])
  if(!is.ultrametric(one_tree) | min(one_tree$edge.length)<0) {
    options(digits=4)
    tree_modified <- one_tree
    tree_modified <- multi2di(tree_modified)
    tree_modified <- force.ultrametric(tree_modified)
    tree_modified$edge.length <- tree_modified$edge.length + 0.0001
    tree_modified <- ladderize(tree_modified)
    # solution copied from https://jonathanchang.org/blog/three-ways-to-check-and-fix-ultrametric-phylogenies/
    tre_node_adjust <- reorder(tree_modified, "postorder")
    N <- Ntip(tree_modified)
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
    write.tree(tre_node_adjust, file=paste0(local_path,"MalietEtAlRuns/MalietEtAl_ClaDS2_trees/", all_trees[i]))
  }
  priors <- setBAMMpriors(one_tree, outfile=NULL)
  generateControlFile(file = paste0("MalietEtAlRuns/BAMM_files/control_files/", label, '_divcontrol.txt'), 
                      params = list( treefile = paste0(local_path,"MalietEtAlRuns/MalietEtAl_ClaDS2_trees/", all_trees[i]),
                                     runInfoFilename = paste0(local_path, "MalietEtAlRuns/BAMM_files/all_results/", label, "_run_info.txt"),
                                     mcmcOutfile = paste0(local_path, "MalietEtAlRuns/BAMM_files/all_results/", label, "_mcmc_out.txt"),
                                     eventDataOutfile = paste0(local_path, "MalietEtAlRuns/BAMM_files/all_results/", label, "_event_data.txt"),
                                     globalSamplingFraction = '1',
                                     numberOfGenerations = '10000000',
                                     overwrite = '1',
                                     lambdaInitPrior = as.numeric(priors['lambdaInitPrior']), lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']), muInitPrior = as.numeric(priors['muInitPrior']), expectedNumberOfShifts = '1'),
                      
  )
}




#----------
all_files <- list.files("MalietEtAlRuns/BAMM_files/control_files", full.names = T)
sink("MalietEtAlRuns/BAMM_files/all_control_files.txt")
for(u in 1:length(all_files)){
  cat(paste0("bamm -c ", getwd(),"/", all_files[u]),"\n")
}
sink()
#----------



mcmcout_files <- list.files("../all_results/")[grep("_mcmc_out.txt", list.files("../all_results/"))]
library(coda)
for(i in 1:length(mcmcout_files)) {
  mcmcout <- read.csv(paste0("../all_results/", mcmcout_files[i]), header=T)
  label <- gsub("_mcmc_out.txt","", mcmcout_files[i])
  #plot(mcmcout$logLik ~ mcmcout$generation)  
  burnstart <- floor(0.1 * nrow(mcmcout))
  postburn <- mcmcout[burnstart:nrow(mcmcout), ]
  tmp1 <- effectiveSize(postburn$N_shifts)
  if(tmp1 < 200) {
    print(paste0("rerun tree ", label))
  }
  tmp2 <- effectiveSize(postburn$logLik)
  if(tmp2 < 200) {
    print(paste0("rerun tree ", label))
  }
}

# checking median height
all_tree_clads <- list()
for(i in 1:length(all_trees)) {
  one_tree <- read.tree(paste0(local_path,"MalietEtAlRuns/MalietEtAl_ClaDS2_trees/", all_trees[i]))
  label <- gsub(".tre$","", all_trees[i])
  all_tree_clads[[i]] <- one_tree
  names(all_tree_clads)[i] <- label
}

totalheight <- function(phy) { 
  return(max(branching.times(phy))) 
}
lapply(all_tree_clads, totalheight)

