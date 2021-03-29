library(hisse)
#empirical.wd <- "~/Desktop/missecomparison/empirical_2021/"
#setwd(empirical.wd)

# Code to run empirical example with Conifers
# Load tree file
tree <- read.tree("Cupressophytes.tre")

# Specify number of possible combos to run
turnover.tries = sequence(12) 
eps.tries = sequence(3)
max.param = length(turnover.tries) + length(eps.tries)
possible.combos = hisse::generateMiSSEGreedyCombinations(max.param = max.param, turnover.tries = turnover.tries, 
                                                         eps.tries = eps.tries, fixed.eps.tries=c(0, 0.9, NA), vary.both=T)

# Fit MiSSE models
model.set = hisse::MiSSEGreedy(tree, f = 0.9, possible.combos = possible.combos, save.file=paste0(empirical.wd,"/conifer_example.Rsave"), 
                               root.type="madfitz", stop.deltaAICc=Inf, n.cores=20, chunk.size=20, turnover.upper=20, trans.upper=10, sann=TRUE, sann.its=1000) #

# Reconstruct rates
model.recons <- as.list(1:length(model.set))
for (model.index in 1:length(model.set)) {
  nturnover <- length(unique(model.set[[model_index]]$turnover))
  neps <- length(unique(model.set[[model_index]]$eps))
  model.recons[[model.index]] <- hisse::MarginReconMiSSE(phy = model.set[[model_index]]$phy, f = ape::Ntip(phy)/total_clade, hidden.states = max(c(nturnover, neps)), 
                                                pars = model.set[[model_index]]$solution, fixed.eps=model.set$fixed.eps , AIC = model.set[[model_index]]$AIC, root.type = "madfitz")
}

tip.rates <- hisse::GetModelAveRates(model.recons, type = "tips")
save(model.set, model.recons, tip.rates, file="Conifer_example.Rsave")


