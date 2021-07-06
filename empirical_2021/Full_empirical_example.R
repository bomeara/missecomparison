#########################################################################
#
#
#
#
#
#
#
#
#
#
#
#
#########################################################################

#########################################################################
# Install hisse
# The first step is to make sure the package hisse is installed in your 
# computer.You can do that easily with following commands:
#########################################################################

install.packages("hisse")
library(hisse)

#########################################################################
# The package hisse is under constant development. If you have problems 
# using the CRAN version # (i.e. that installed with the commands above),
# it may be wise to get the most up to date version of hisse by downloading 
# it directly from github with the following commands:
#########################################################################
install.packages("devtools")
devtools::install_github("thej022214/hisse")

# TO BE REMOVED FOR SUBMISSION #
# empirical.wd <- "~/Desktop/misse_mme_paper/missecomparison/empirical_2021/"
# setwd(empirical.wd)
# TO BE REMOVED FOR SUBMISSION #

#########################################################################
# Code to run empirical example with Eucalypts
#########################################################################
# Let's now replicate the example with Eucalypts from Vasconcelos et al. 
# (in prep.). The goal here is to prepare all arguments that go in the 
# functions MiSSEGreedy and MarginReconMiSSE. First of all, let's load
# the Eucalypts tree from Thornhill et al. (2019) "A dated molecular 
# perspective of eucalypt taxonomy, evolution and diversification". 
# Australian Systematic Botany, 32(1), 29-48.)
#########################################################################
## Individual models can be fitted out of the greedy approach with MiSSE().

#########################################################################
# (1) Load tree file
#########################################################################

tree_wd <- "./trees"
phy <- read.tree("Eucalypts.tre")

# TO BE REMOVED FOR SUBMISSION #
# phy <- extract.clade(phy, node=1340) # for speed
# TO BE REMOVED FOR SUBMISSION #

#########################################################################
# (2) Specify number of possible combos to run
# A combo here is defined as..
#########################################################################
turnover.tries = sequence(3) 
eps.tries = sequence(3)
max.param = length(turnover.tries) + length(eps.tries)
fixed.eps.tries = NA
possible.combos = generateMiSSEGreedyCombinations(max.param=max.param, 
                                                  turnover.tries=turnover.tries, 
                                                  eps.tries=eps.tries, 
                                                  fixed.eps.tries=fixed.eps.tries)

#########################################################################
# And this is how the possible.combos data.frame should look like:
#########################################################################
head(possible.combos)

#########################################################################
# (3) Specify a sampling fraction
# Sampling fraction is ...
# In our case, the total number of species of Eucalypts is estimated at.
# So the sampling fraction can be set as xx
# for more information on why to use a global sampling fraction instead of 
# a clade specific one, see 
#########################################################################
f = 0.87

#########################################################################
# (4) Specify a name to save the progress of the greedy
# IMPORTANT NOTE: Remember that you should change the name of this object 
# for different runs or the file will be overwritten
#########################################################################
save.file = "Eucalypts_fit_fast.Rsave"

#########################################################################
# (5) Specify a delta AIC to stop MiSSEgreedy
#  
#########################################################################
stop.deltaAICc = 10

#########################################################################
# (6) Specify number of cores to run - to find out, use
# To run in my macbook, I will use 4 cores, but for big big trees it is comendable that
# you increase the number of cores when in parallel.
#########################################################################
n.cores = 4

#########################################################################
# (6) Specify the chunk size of how many models per chunk
#########################################################################
chunk.size = 5

#########################################################################
# (7) FIT MiSSE MODEL 
# Specify the bounds 
# all other parameters set to default
# NOTE: This command can take a while to run depending on te 
#########################################################################
model.set = MiSSEGreedy(phy=phy, # the phylogeny 
                        f=f, # sampling fraction
                        possible.combos=possible.combos, # possible combinations of models
                        save.file=save.file, # the name of the file to save
                        stop.deltaAICc=stop.deltaAICc, # the deltaAIC to stop running
                        n.cores=n.cores, # number of cores
                        chunk.size=chunk.size, sann=F) # size of the "chunk" of models

#########################################################################
# (8) Prune redundant
# Your set of models should be a list and look like this:
#########################################################################
class(model.set)
length(model.set)

#########################################################################
# You should then prune the redundant models with the following command:
#
#########################################################################
model.set_pruned <- PruneRedundantModels(model.set)
class(model.set_pruned)
length(model.set_pruned)

#########################################################################
# (9) Reconstruct rates
# Because we want to reconstruct all non-redundant models, let's loop over all
# possibilities with a for loop:
#########################################################################
model.recons <- as.list(1:length(model.set_pruned))
for (model_index in 1:length(model.set_pruned)) {
  nturnover <- length(unique(model.set_pruned[[model_index]]$turnover))
  neps <- length(unique(model.set_pruned[[model_index]]$eps))
  model.recons[[model_index]] <- hisse::MarginReconMiSSE(phy = model.set_pruned[[model_index]]$phy, f = 0.87, hidden.states = max(c(nturnover, neps)), 
                                                         pars = model.set_pruned[[model_index]]$solution, fixed.eps=model.set_pruned$fixed.eps , 
                                                         AIC = model.set_pruned[[model_index]]$AIC, root.type = "madfitz",n.cores=n.cores)   
}

# Calculating marginal probabilities for 40 internal nodes... 
# Finished. Calculating marginal probabilities for 41 tips... 
# Done. 

class(model.recons)
length(model.recons)

#########################################################################
# (10) Get tip rates
# Finally, we can use the reconstruct models to estimate the rates at the tips
#########################################################################

tip.rates <- GetModelAveRates(model.recons, type = "tips")

#########################################################################
# We can also set to get the tips based only on the "best" model (though this is not commendable)
# To do that mannualy, you can uncomment the following lines:
#########################################################################

# AIC.weights <- rep(0, length(model.recons))
# AIC.weights[which.max(GetAICWeights(model.recons))] <- 1
# tip.rates.best <- GetModelAveRates(model.recons, AIC.weights=AIC.weights, type = "tips")

#########################################################################
# The tip.rates should look like this:
#########################################################################

head(tip.rates)

#########################################################################
# Recommended to save all files:
save(model.set_pruned, 
     model.recons, 
     tip.rates, 
     possible.combos, # MiSSEgreedy will overwrite this object as it runs, showing the updated AIC and lokLik it may be interesting to have it
     file="Eucalypts_example.Rsave")
#save(model.set, model.recons, tip.rates.best, possible.combos, file="Eucalypts_example.best.Rsave")


#########################################################################

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#                            Post analyses
#########################################################################

# First, let's reload the results of the run that we saved in the previous part

load("Eucalypts_example.Rsave")
#

# Some aditional options for ploting MiSSE tip rates are given as follows:
# the file "plot_functions_for_MiSSE" can be found in the supplementary of
# Vasconcelos et al.

# Some plots:

# Let's see the improvement in AIC and lokLik with increasing in model complexity
plot(x=1:length(model.set_pruned), y=unlist(lapply(model.set_pruned[1:length(model.set_pruned)], "[[", "loglik")), xlab="comb", ylab="loglik")
plot(x=1:length(model.set_pruned), y=unlist(lapply(model.set_pruned[1:length(model.set_pruned)], "[[", "AICc")), xlab="comb", ylab="AICc")

# Where the X-axis corresponds to the models in the table:


#boxplot(tip.rates_1$turnover, tip.rates_1$net.div, 
#        tip.rates_1$speciation,tip.rates_1$extinct.frac, 
#        tip.rates_1$extinction, main="tip rates run 1", 
#        names=c("turnover", "netdiv", "speciation", "eps", "extinction"))

painted.tree <- hisse::plot.misse.states(x = model.recons, 
                                          rate.param = "speciation", type = "fan",  
                                          show.tip.label = F, fsize=0.9) 



# get confidence intervals with dentist?


