
setwd("~/Desktop/misse_mme_paper/missecomparison/Simulation_study/")
local_path <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/Simulation_study/" #local
lab_computer_path <- "/home/tvasconcelos/missecomparison/ClaDScomparison" #labcomputer

all_trees <- list.files(paste0(local_path, "MalietEtAlRuns/MalietEtAl_ClaDS2_trees/"))
all_trees <- all_trees[grep(".tre$", all_trees)]

# install.packages("BAMMtools")
library(BAMMtools)
library(phytools)

all_tree_clads <- list()
for(i in 1:length(all_trees)) {
  one_tree <- read.tree(paste0(local_path,"MalietEtAlRuns/MalietEtAl_ClaDS2_trees/", all_trees[i]))
  label <- gsub(".tre$","", all_trees[i])
  if(!is.ultrametric(one_tree)) {
    #plot(one_tree, cex=0.5)
    #max(one_tree$edge.length)
    axisPhylo()
    dev.off()
    one_tree <- force.ultrametric(one_tree, file=paste0(local_path,"/original_trees/trees/ClaDS2/", all_trees[i]))
    write.tree(one_tree, file=paste0(local_path,"/original_trees/trees/ClaDS2/", all_trees[i]))
  }
  priors <- setBAMMpriors(one_tree, outfile=NULL)
  generateControlFile(file = paste0(label, '_divcontrol.txt'), 
                      params = list( treefile = paste0(local_path,"/original_trees/trees/ClaDS2/", all_trees[i]),
                                     runInfoFilename = paste0(local_path, "/BAMM_ClaDS_comparison/all_results/", label, "_run_info.txt"),
                                     mcmcOutfile = paste0(local_path, "/BAMM_ClaDS_comparison/all_results/", label, "_mcmc_out.txt"),
                                     eventDataOutfile = paste0(local_path, "/BAMM_ClaDS_comparison/all_results/", label, "_event_data.txt"),
                                                              globalSamplingFraction = '1',
                                                              numberOfGenerations = '10000000',
                                                              overwrite = '1',
                                                              lambdaInitPrior = as.numeric(priors['lambdaInitPrior']), lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']), muInitPrior = as.numeric(priors['muInitPrior']), expectedNumberOfShifts = '1'),

                      )
}



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


#get bamm priors to supply to control file priors <- setBAMMpriors(whales, outfile = NULL)
all_files <- list.files()
sink("../all_control_files.txt")
for(u in 1:length(all_files)){
  cat(paste0("bamm -c ", getwd(),"/", all_files[u]),"\n")
}
sink()




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


