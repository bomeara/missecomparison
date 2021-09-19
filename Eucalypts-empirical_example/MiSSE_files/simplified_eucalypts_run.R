#########################################################################
# Install hisse
# The first step is to make sure the package hisse is installed in your 
# computer.You can do that easily with following commands:
#########################################################################
library(hisse)

#########################################################################
# The package hisse is under constant development. If you have problems 
# using the CRAN version, it may be wise to get the most up to date version,
# of hisse by downloading it directly from github with the following commands:
#########################################################################
# install.packages("devtools")
# devtools::install_github("thej022214/hisse")

#########################################################################
# Code to run empirical example with Eucalypts
#########################################################################
# Load the Eucalypts tree from Thornhill et al. (2019) ("A dated molecular 
# perspective of eucalypt taxonomy, evolution and diversification". 
# Australian Systematic Botany, 32(1), 29-48.)
#########################################################################

#########################################################################
# (1) Load tree files
#########################################################################
phy_ml1 <- read.tree("ML1_modified-Thornhill_et_al_2019.tre")

#########################################################################
# (2) Specify number of possible combos to run
#########################################################################
possibilities_ml1 = generateMiSSEGreedyCombinations(max.param=round(ape::Ntip(phy_ml1)/10))

#########################################################################
# (3) Specify a sampling fraction
#########################################################################
f = 0.85 # Based on the original paper

#########################################################################
# (4) Specify a name to save the progress of the greedy
#########################################################################
save.file_ml1 = "Eucalypts_fit_ml1.Rsave"

#########################################################################
# (5) Specify a delta AICc to stop MiSSEgreedy
#########################################################################
stop.deltaAICc = 10

#########################################################################
# (6) Specify number of cores 
#########################################################################
n.cores = 20

#########################################################################
# (7) Specify the chunk size of how many models to run per chunk
#########################################################################
chunk.size = 20

#########################################################################
# (8) Fit MiSSE models
#########################################################################
model.set_ml1 = MiSSEGreedy(phy=phy_ml1, # the phylogeny 
                        f=f, # sampling fraction
                        possible.combos=possibilities_ml1, # possible combinations of models
                        save.file=save.file_ml1, # the name of the file to save
                        stop.deltaAICc=stop.deltaAICc, # the deltaAIC to stop running
                        n.cores=n.cores, # number of cores
                        chunk.size=chunk.size) # size of the "chunk" of models

#########################################################################
# (8) Prune redundant
#########################################################################
load("Eucalypts_fit_ml1.Rsave")
model.set_ml1 <- misse.list
possible.combos_ml1 <- final.combos
model.set_pruned_ml1 <- PruneRedundantModels(model.set_ml1)

#########################################################################
# (9) Reconstruct rates and get tip rates
#########################################################################
n.cores=96

#########################################################################
model.recons_ml1 <- as.list(1:length(model.set_pruned_ml1))
for (model_index in 1:length(model.set_pruned_ml1)) {
  nturnover <- length(unique(model.set_pruned_ml1[[model_index]]$turnover))
  neps <- length(unique(model.set_pruned_ml1[[model_index]]$eps))
  model.recons_ml1[[model_index]] <- hisse::MarginReconMiSSE(phy = model.set_pruned_ml1[[model_index]]$phy, f = 0.85, hidden.states = max(c(nturnover, neps)), 
                                                               pars = model.set_pruned_ml1[[model_index]]$solution, fixed.eps=model.set_pruned_ml1$fixed.eps , 
                                                               AIC = model.set_pruned_ml1[[model_index]]$AIC, root.type = "madfitz",n.cores=n.cores)   
}
save(model.recons_ml1, file="recon.Rsave")

#########################################################################
tip.rates_ml1 <- GetModelAveRates(model.recons_ml1, type = "tips")
#########################################################################

#########################################################################
# (10) Save all files
#########################################################################
save(model.set_pruned_ml1, 
     model.recons_ml1, 
     tip.rates_ml1, 
     possible.combos_ml1, 
     file="Eucalypts_example_ml1_tree.Rsave")


