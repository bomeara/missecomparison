#########################################################################
# Quick example on how to run MiSSE 
# [The whole script should take approx. 30min to run]
#########################################################################

#########################################################################
# (1) Install package 'hisse'.
# The first step is to make sure hisse is installed in your computer. You 
# can do that easily with following commands:
#########################################################################
#install.packages("hisse")
library(hisse)

#########################################################################
# Package hisse is under constant development. If you have problems 
# using the CRAN version (i.e. the one installed with the commands above),
# it may be wise to try the most up to date version of hisse by downloading 
# it directly from github with the following commands:
#########################################################################
#install.packages("devtools")
#devtools::install_github("thej022214/hisse", force=T)

#########################################################################
# (2) Load tree file.
# Let's now run a quick MiSSE analysis. The goal here is to prepare all 
# arguments that will go in the functions MiSSEGreedy (which will fit a 
# series of MiSSE models of increasing complexities) and MarginReconMiSSE 
# (which will reconstruct the rates at the tips using a marginal reconstruction 
# algorithm). First of all, let's load the tree file that will be used 
# in this example. Here we analyze the Lupinus tree from Drummond 
# et al. (2012) [Systematic Biology, 61(3), 443-460].
#########################################################################
phy <- read.tree("Lupinus.tre")

#########################################################################
# (3) Specify number of possible combinations of rates classes to run.
# MiSSE can run models with up to 26 different rate classes. In this step,
# we will set the possible combinations of how many different net-turnover 
# and extinction-fractions we can have in our set of MiSSE models. A model 
# with 2 turnover rates and 2 extinction fraction, for example, is a model 
# with two rate classes of variable rates of net-turnover and extinction 
# fractions, whereas a model with 2 turnover rates and 1 extinction fraction 
# is also a model with two rate classes, but the same extinction fraction 
# is shared between rate classes. 
#########################################################################
max.param = max(4, round(ape::Ntip(phy)/10))
set.seed(42)
possible.combos = generateMiSSEGreedyCombinations(max.param=max.param, vary.both=TRUE, fixed.eps.tries=NA)

#########################################################################
# And this is how the possible.combos data.frame should look like:
#####################################################q####################
head(possible.combos)
# turnover eps fixed.eps
# 1         7   1        NA
# 2         4   3        NA
# 3         1   1        NA
# 4         4   2        NA
# 5         4   1        NA
# 6         1   2        NA

#########################################################################
# (4) Specify a sampling fraction.
# Sampling fraction is the proportion of species that were sampled in a 
# tree relative to how many species we believe that exist in that clade In 
# our example, the total number Lupinus species is estimated in 250, so 
# let's use a sampling fraction of 0.45 for our tree of 120 species. 
# For more information on why to use a global sampling fraction instead of 
# a clade specific one, see [link to Jeremy's vignette]
#########################################################################
f = 0.45

#########################################################################
# (5) Specify a file name to save the progress of MiSSEGreedy.
# Here we give a name for the .Rsave file that will save MiSSEGreedy's 
# progress in our working directory. IMPORTANT: Remember that you should 
# change the name of this object between runs or the file will be 
# overwritten.
#########################################################################
save.file = "Lupinus_fit_fast.Rsave"

#########################################################################
# (6) Specify a deltaAICc to stop MiSSEGreedy.
# MiSSEGreedy will fit models of increasing complexities following the order
# presented in the possible.combos object. However, very complex models, 
# for example, will likely be a poor fit for most small trees. MiSSEGreedy 
# will then not necessarily fit all models of the possible.combos object 
# and will stop running when the likelihood and AICc of new models stop 
# improving significantly. We then want to inform MiSSEGreedy when we think 
# the AICc is not improving anymore with the argument stop.deltaAICc. 
# Because we are using a relatively small tree in our example (120 tips), 
# we can also use a small value for the stop.deltaAICc (but this number should 
# be larger for larger trees). 
#########################################################################
stop.deltaAICc = 2

#########################################################################
# (7) Specify number of cores to run MiSSE.
# To find out how many cores you have in your computer, you can use:
#########################################################################
#install.packages("parallel")
parallel::detectCores()
# [1] 8 

#########################################################################
# My macbook has 8 cores, but I will use 4 cores to run MiSSE. You may 
# want to use a more powerful computer and more cores if you have a larger
# tree.
#########################################################################
n.cores = 20

#########################################################################
# (8) Specify the chunk size.
# This number refers to how many different models MiSSEGreedy will fit at
# each time before checking if the lowest AICc of the set has a higher 
# deltaAICc than the number defined as stop.deltaAICc. Here we are going 
# to set the chunk.size to 4, which means that MiSSEGreedy will fit 4 models 
# before checking if the fit stopped improving.
#########################################################################
chunk.size = 5

#########################################################################
# (9) Fit MiSSE models using MiSSEGreedy.
# The following command will take all the arguments set in the steps above
# to fit a series of MiSSE models using a greedy algorithm. In this run, 
# all other arguments in the function are set to default, but you should 
# look at the help page of MiSSEGreedy to set the other arguments according 
# to your needs.
#########################################################################
#?MiSSEGreedy()

#########################################################################
# Note that this command can take a while to run depending on the size of 
# the tree and the number of cores available in your computer. In my computer 
# (a MacBook Pro 2.4 GHz Quad-Core Intel Core i5), it takes approximately 
# 20 minutes using 4 cores. IMPORTANT NOTE: Note that we set sann=FALSE 
# in this run for speed since sann=TRUE (the default) tends to make runs 
# slower. However, it is strongly advised to use sann=TRUE for real world 
# optimizations, as it tends to give more accurate parameter estimates.
#########################################################################
model.set = MiSSEGreedy(phy=phy, # the phylogeny 
                        f=f, # sampling fraction
                        possible.combos=possible.combos, # possible combinations of models
                        save.file=save.file, # the name of the file to save
                        stop.deltaAICc=stop.deltaAICc, # the deltaAIC to stop running
                        n.cores=n.cores, # number of cores
                        chunk.size=chunk.size, # number of models to run at each time
                        sann=FALSE) # IMPORTANT IMPORTANT IMPORTANT - This argument is here set to F for speed, 
                                # but it should be set to TRUE for your actual runs. See above. 

#########################################################################
# You should see this printed on your console:
#########################################################################
# Starting at 2021-08-30 16:33:06
# running on 4 cores.
# turnover eps fixed.eps
# 1        1   1        NA
# 2        2   1        NA
# 3        1   4        NA
# 4        2   2        NA

#########################################################################
# As you will see, even though we gave MiSSEGreedy several possible 
# rate-class combinations, it stopped to run as soon as the AICc reached 
# the specified deltaAICc in contrast to the "best model". In our example,
# this happened after the second "chunk" of models.
# In the next steps, we will work with this set of models that have 
# relatively low AICc.
#########################################################################

#########################################################################
# (10) Prune redundant
# Your set of models should be a list and look like this:
#########################################################################
class(model.set) 
# [1] "list"
length(model.set) 
# [1] 8

#########################################################################
# You should then prune the redundant models with the following command. 
# More details about this step, use: 
#########################################################################
# ?PruneRedundantModels()
#########################################################################
model.set_pruned <- PruneRedundantModels(model.set)
length(model.set_pruned) 
# [1] 8

#########################################################################
# In this case, we kept our same 8 models after prunning, though this number
# may decrease in cases where there are redundant models in the set.
#########################################################################

#########################################################################
# (11) Reconstruct rates.
# Because we want to estimate tip-rates using all non-redundant models, 
# let's loop over all possibilities (this should also take a few minutes to run):
#########################################################################
model.recons <- as.list(1:length(model.set_pruned))
for (model_index in 1:length(model.set_pruned)) {
  nturnover <- length(unique(model.set_pruned[[model_index]]$turnover))
  neps <- length(unique(model.set_pruned[[model_index]]$eps))
  model.recons[[model_index]] <- hisse::MarginReconMiSSE(phy = model.set_pruned[[model_index]]$phy, f = 0.45, hidden.states = max(c(nturnover, neps)), 
                                                         pars = model.set_pruned[[model_index]]$solution, fixed.eps=model.set_pruned$fixed.eps , 
                                                         AIC = model.set_pruned[[model_index]]$AIC, root.type = "madfitz",n.cores=n.cores)   
}

#########################################################################
# Something like this should be printed on your console:
#########################################################################
# Calculating marginal probabilities for 119 internal nodes... 
# Finished. Calculating marginal probabilities for 120 tips... 
# Done. 

#########################################################################
# (12) Get tip-rates.
# Finally, we can use the reconstructed models to estimate rates of speciation, 
# extinction, net-diversification, net-turnover and extinction fraction at 
# the tips with:
#########################################################################
tip.rates <- GetModelAveRates(model.recons, type = "tips")

#########################################################################
# We can also get the tip-rates based only on the "best" model 
# (though this is not always commendable). To do that 'mannualy', you can 
# uncomment the following lines:
#########################################################################
# AIC.weights <- rep(0, length(model.recons))
# AIC.weights[which.max(GetAICWeights(model.recons))] <- 1
# tip.rates.best <- GetModelAveRates(model.recons, AIC.weights=AIC.weights, type = "tips")

#########################################################################
# We believe that these tip rate values can be used in an analogous way to 
# continuous traits evolving in a tree (e.g. as we did with Eucalypts in 
# Vasconcelos et al. in prep). The tip.rates should look like this:
#########################################################################
head(tip.rates)
#                  taxon state  turnover   net.div speciation extinct.frac extinction
# 1      Lupinus_diffusus     0 0.8371675 0.2901501  0.5636588    0.4678924  0.2735087
# 2    Lupinus_cumulicola     0 0.8404156 0.2924021  0.5664089    0.4672625  0.2740068
# 3      Lupinus_villosus     0 0.8404156 0.2924021  0.5664089    0.4672625  0.2740068
# 4    Lupinus_micranthus     0 0.8132777 0.2757810  0.5445294    0.4705892  0.2687484
# 5 Lupinus_angustifolius     0 0.8211815 0.2795363  0.5503589    0.4704875  0.2708226
# 6        Lupinus_luteus     0 0.8258186 0.2819912  0.5539049    0.4700129  0.2719137

#########################################################################
# Recommended to save all files:
save(model.set_pruned, 
     model.recons, 
     tip.rates, 
     possible.combos,
     file="Lupinus_example.Rsave")

#########################################################################
# You can also visualize the rates evolving in the tree with the following command, 
# though rates through time *should not* be interpreted literally. The painting
# is just to get a "feeling" for the model.
#painted.tree <- hisse::plot.misse.states(x = model.recons, 
#                                         rate.param = "speciation", type = "phylo", show.tip.label = F) 

