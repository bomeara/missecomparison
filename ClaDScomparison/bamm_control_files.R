
# setwd("~/Desktop/misse_mme_paper/missecomparison/ClaDScomparison/BAMM_ClaDS_comparison/control_files")
local_path <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison" #local
lab_computer_path <- "/home/tvasconcelos/missecomparison/ClaDScomparison" #labcomputer

all_trees <- list.files(paste0(base.dir, "/original_trees/trees/ClaDS2"))
all_trees <- all_trees[grep(".tre$", all_trees)]

# install.packages("BAMMtools")
library(BAMMtools)

for(i in 1:length(all_trees)) {
  one_tree <- read.tree(paste0(lab_computer_path,"/original_trees/trees/ClaDS2/", all_trees[i]))
  priors <- setBAMMpriors(one_tree, outfile=NULL)
  generateControlFile(file = paste0(all_trees[i], '_divcontrol.txt'), 
                      params = list( treefile = paste0(lab_computer_path,"/original_trees/trees/ClaDS2/", all_trees[i]),
                                     runInfoFilename = paste0(lab_computer_path, "/BAMM_ClaDS_comparison/all_results/", all_trees[i], "_run_info.txt"),
                                     mcmcOutfile = paste0(lab_computer_path, "/BAMM_ClaDS_comparison/all_results/", all_trees[i], "_mcmc_out.txt"),
                                     eventDataOutfile = paste0(lab_computer_path, "/BAMM_ClaDS_comparison/all_results/", all_trees[i], "_event_data.txt"),
                                                              globalSamplingFraction = '1',
                                                              numberOfGenerations = '1000000',
                                                              overwrite = '1',
                                                              lambdaInitPrior = as.numeric(priors['lambdaInitPrior']), lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']), muInitPrior = as.numeric(priors['muInitPrior']), expectedNumberOfShifts = '1'),

                      )
}

#get bamm priors to supply to control file priors <- setBAMMpriors(whales, outfile = NULL)
