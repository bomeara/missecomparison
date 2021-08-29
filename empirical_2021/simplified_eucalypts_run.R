#########################################################################
# Install hisse
# The first step is to make sure the package hisse is installed in your 
# computer.You can do that easily with following commands:
#########################################################################
library(hisse)
library(phytools)

#########################################################################
# The package hisse is under constant development. If you have problems 
# using the CRAN version # (i.e. that installed with the commands above),
# it may be wise to get the most up to date version of hisse by downloading 
# it directly from github with the following commands:
#########################################################################
#install.packages("devtools")
#devtools::install_github("thej022214/hisse")

# TO BE REMOVED FOR SUBMISSION #
# empirical.wd <- "~/Desktop/misse_mme_paper/missecomparison/empirical_2021/"
# setwd(empirical.wd)
# setwd("~/missecomparison/empirical_2021") # lab computer
# TO BE REMOVED FOR SUBMISSION #

#########################################################################
# Code to run empirical example with Eucalypts
#########################################################################
# Let's now replicate the example with Eucalypts from Vasconcelos et al. 
# (in prep.). The goal here is to prepare all arguments that go in the 
# functions MiSSEGreedy and MarginReconMiSSE. First of all, let's load
# the Eucalypts tree from Thornhill et al. (2019) ("A dated molecular 
# perspective of eucalypt taxonomy, evolution and diversification". 
# Australian Systematic Botany, 32(1), 29-48.)
#########################################################################

#########################################################################
# (1) Load tree files
#########################################################################
#phy_bayes <- force.ultrametric(read.tree("trees/Eucalypts_Bayes.tre"))
#phy_ml1 <- force.ultrametric(read.tree("trees/Eucalypts_ML1.tre"))
phy_ml1 <- read.tree("ML1_modified.tre")
#phy_ml2 <- force.ultrametric(read.tree("trees/Eucalypts_ML2.tre"))

#########################################################################
# (2) Specify number of possible combos to run
#########################################################################
#possibilities_bayes = generateMiSSEGreedyCombinations(max.param=round(ape::Ntip(phy_bayes)/50), vary.both=TRUE, fixed.eps.tries = NA)
possibilities_ml1 = generateMiSSEGreedyCombinations(max.param=round(ape::Ntip(phy_ml1)/10))
#possibilities_ml2 = generateMiSSEGreedyCombinations(max.param=round(ape::Ntip(phy_ml2)/50), vary.both=TRUE, fixed.eps.tries = NA)

#########################################################################
# (3) Specify a sampling fraction
#########################################################################
f = 0.85

#########################################################################
# (4) Specify a name to save the progress of the greedy
#########################################################################
#save.file_bayes = "Eucalypts_fit_bayes.Rsave"
save.file_ml1 = "Eucalypts_fit_ml1.Rsave"
#save.file_ml2 = "Eucalypts_fit_ml2.Rsave"

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
#model.set_bayes = MiSSEGreedy(phy=phy_bayes, # the phylogeny 
#                        f=f, # sampling fraction
#                        possible.combos=possibilities_bayes, # possible combinations of models
#                        save.file=save.file_bayes, # the name of the file to save
#                        stop.deltaAICc=stop.deltaAICc, # the deltaAIC to stop running
#                        n.cores=n.cores, # number of cores
#                        chunk.size=chunk.size) # size of the "chunk" of models

model.set_ml1 = MiSSEGreedy(phy=phy_ml1, # the phylogeny 
                        f=f, # sampling fraction
                        possible.combos=possibilities_ml1, # possible combinations of models
                        save.file=save.file_ml1, # the name of the file to save
                        stop.deltaAICc=stop.deltaAICc, # the deltaAIC to stop running
                        n.cores=n.cores, # number of cores
                        chunk.size=chunk.size) # size of the "chunk" of models

#model.set_ml2 = MiSSEGreedy(phy=phy_ml2, # the phylogeny 
#                        f=f, # sampling fraction
#                        possible.combos=possibilities_ml2, # possible combinations of models
#                        save.file=save.file_ml2, # the name of the file to save
#                        stop.deltaAICc=stop.deltaAICc, # the deltaAIC to stop running
#                        n.cores=n.cores, # number of cores
#                        chunk.size=chunk.size) # size of the "chunk" of models

#########################################################################
# (8) Prune redundant
#########################################################################
#load("Eucalypts_fit_bayes.Rsave")
#model.set_bayes <- misse.list
#possible.combos_bayes <- possible.combos
#model.set_pruned_bayes <- PruneRedundantModels(model.set_bayes)

load("Eucalypts_fit_ml1.Rsave")
model.set_ml1 <- misse.list
possible.combos_ml1 <- final.combos
model.set_pruned_ml1 <- PruneRedundantModels(model.set_ml1)

#load("Eucalypts_fit_ml2.Rsave")
#model.set_ml2 <- misse.list
#possible.combos_ml2 <- possible.combos
#model.set_pruned_ml2 <- PruneRedundantModels(model.set_ml2)

#########################################################################
# (9) Reconstruct rates and get tip rates
#########################################################################
n.cores=96

#model.recons_bayes <- as.list(1:length(model.set_pruned_bayes))
#for (model_index in 1:length(model.set_pruned_bayes)) {
#  nturnover <- length(unique(model.set_pruned_bayes[[model_index]]$turnover))
#  neps <- length(unique(model.set_pruned_bayes[[model_index]]$eps))
#  model.recons_bayes[[model_index]] <- hisse::MarginReconMiSSE(phy = model.set_pruned_bayes[[model_index]]$phy, f = 0.85, hidden.states = max(c(nturnover, neps)), 
#                                                         pars = model.set_pruned_bayes[[model_index]]$solution, fixed.eps=model.set_pruned_bayes$fixed.eps , 
#                                                         AIC = model.set_pruned_bayes[[model_index]]$AIC, root.type = "madfitz",n.cores=n.cores)   
#}
#tip.rates_bayes <- GetModelAveRates(model.recons_bayes, type = "tips")
#save(model.set_pruned_bayes, 
#     model.recons_bayes, 
#     tip.rates_bayes, 
#     possible.combos_bayes, # MiSSEgreedy will overwrite this object as it runs, showing the updated AIC and lokLik it may be interesting to have it
#     file="Eucalypts_example_bayes_tree.Rsave")

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

tip.rates_ml1 <- GetModelAveRates(model.recons_ml1, type = "tips")
save(model.set_pruned_ml1, 
     model.recons_ml1, 
     tip.rates_ml1, 
     possible.combos_ml1, # MiSSEgreedy will overwrite this object as it runs, showing the updated AIC and lokLik it may be interesting to have it
     file="Eucalypts_example_ml1_tree.Rsave")

#########################################################################
#model.recons_ml2 <- as.list(1:length(model.set_pruned_ml2))
#for (model_index in 1:length(model.set_pruned_ml2)) {
#  nturnover <- length(unique(model.set_pruned_ml2[[model_index]]$turnover))
#  neps <- length(unique(model.set_pruned_ml2[[model_index]]$eps))
#  model.recons_ml2[[model_index]] <- hisse::MarginReconMiSSE(phy = model.set_pruned_ml2[[model_index]]$phy, f = 0.85, hidden.states = max(c(nturnover, neps)), 
#                                                               pars = model.set_pruned_ml2[[model_index]]$solution, fixed.eps=model.set_pruned_ml2$fixed.eps , 
#                                                               AIC = model.set_pruned_ml2[[model_index]]$AIC, root.type = "madfitz",n.cores=n.cores)   
#}
#tip.rates_ml2 <- GetModelAveRates(model.recons_ml2, type = "tips")
#save(model.set_pruned_ml2, 
#     model.recons_ml2, 
#     tip.rates_ml2, 
#     possible.combos_ml2, # MiSSEgreedy will overwrite this object as it runs, showing the updated AIC and lokLik it may be interesting to have it
#     file="Eucalypts_example_ml2_tree.Rsave")


