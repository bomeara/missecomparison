
#-----------

MarginReconMiSSE_full <- function(model.set, possible.combos, f=1, root.type = "madfitz", models.to.recon=c("all","best"), prune.redundant=TRUE) {
  misse_fit <- model.set
  reference_table <- possible.combos[!is.na(possible.combos$lnL),]
  
  # categorize by class of models
  eps0 <- eps0.9 <- eps_var <- c()
  for(model_index in sequence(nrow(reference_table))) {
    if(is.na(reference_table[model_index,"fixed.eps"])) { eps_var <- c(eps_var, model_index)
    } else {
      if(reference_table[model_index,"fixed.eps"] == 0) { eps0 <- c(eps0, model_index)
      }
      if(reference_table[model_index,"fixed.eps"] == 0.9) { eps0.9 <- c(eps0.9, model_index) }
    }
  }
  models_eps0 <- misse_fit[eps0]; models_eps0.9 <- misse_fit[eps0.9]; models_eps_var <- misse_fit[eps_var]
  if(prune.redundant){
    models_eps0 <- hisse::PruneRedundantModels(models_eps0)
    models_eps0.9 <- hisse::PruneRedundantModels(models_eps0.9)
    models_eps_var <- hisse::PruneRedundantModels(models_eps_var)
  } 
  
  if(models.to.recon=="best"){
  best_models <- c(models_eps0[get.95.aicw(models_eps0)],models_eps0.9[get.95.aicw(models_eps0.9)],models_eps_var[get.95.aicw(models_eps_var)])
  } else {
    best_models <-  c(models_eps0 ,models_eps0.9, models_eps_var)
  }

  # Reconstruct rates
  model.recons <- as.list(1:length(best_models))
  cat("Running marginal reconstruction for",length(best_models),"models.")
  for (model_index in 1:length(best_models)) {
    nturnover <- length(unique(best_models[[model_index]]$turnover))
    neps <- length(unique(best_models[[model_index]]$eps))
    model.recons[[model_index]] <- hisse::MarginReconMiSSE(phy = best_models[[model_index]]$phy, f = 0.87, hidden.states = max(c(nturnover, neps)), 
                                                           pars = best_models[[model_index]]$solution, fixed.eps=best_models$fixed.eps , AIC = best_models[[model_index]]$AIC, root.type = "madfitz")
  }
  return(model.recons)
}

get.95.aicw <- function(misse.models) {
  aicw = 0
  models_to_weigh <- NULL
  tmp_aicw <- hisse::GetAICWeights(misse.models,criterion="AIC")
  while(aicw < 0.95) {
      maxweight <- max(tmp_aicw)
      models_to_weigh <- c(models_to_weigh, which(tmp_aicw==maxweight))
      tmp_aicw[models_to_weigh] <- 0
      aicw <- sum(aicw, maxweight) 
};return(models_to_weigh);}


#----------------------------------------------
#----------------------------------------------
#----------------------------------------------

library(hisse)
#empirical.wd <- "~/Desktop/missecomparison/empirical_2021/"
#setwd(empirical.wd)

# Code to run empirical example with Conifers
# Load tree file
tree <- read.tree("Eucalypts.tre")
#plot(tree, cex=0.2)
#write.csv(tree$tip.label, file="Eucalypt_heights.csv")


# Specify number of possible combos to run
turnover.tries = sequence(12) 
eps.tries = sequence(3)
max.param = length(turnover.tries) + length(eps.tries)

possible.combos = hisse::generateMiSSEGreedyCombinations(max.param = max.param, turnover.tries = turnover.tries, 
                                                         eps.tries = eps.tries, fixed.eps.tries=c(0, 0.9, NA), vary.both=T)

# Fit MiSSE models
model.set = hisse::MiSSEGreedy(tree, f = 0.87, possible.combos = possible.combos, save.file="Eucalypts_fit.Rsave", 
                               root.type="madfitz", stop.deltaAICc=10, n.cores=5, chunk.size=5, turnover.upper=20, trans.upper=10, sann=TRUE, sann.its=1000) #

# load("Eucalypts_fit.Rsave")
# model.set <- misse.list

recon_models <- MarginReconMiSSE_full(model.set, possible.combos, f=1, root.type = "madfitz", models.to.recon=c("all"), prune.redundant=F)

tip.rates <- hisse::GetModelAveRates(model.recons, type = "tips")
save(model.set, model.recons, tip.rates, file="Eucalypts_example.Rsave")


# 
# PLOTS

heights <- read.csv("Eucalypt_heights.csv")[,c(2,3)]
heights <- heights[heights$max_height_m!="no_info_yet",]
heights$max_height_m <- as.numeric(heights$max_height_m)

hist(heights$max_height_m, breaks = 50)





