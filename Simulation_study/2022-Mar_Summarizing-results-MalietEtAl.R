############################################################
# Code to summarize results from simulation studies
############################################################
# rm(list=ls())
############################################################
# MalietEtAl results (TB, ND, DR, BAMM, MiSSE and ClaDS)
############################################################

# setwd("~/Desktop/misse_mme_paper/missecomparison/Simulation_study/")
library(Metrics)
library(tidyverse)
library(magrittr)
library(stats)
library(progress)
library(PMCMR)
library(BAMMtools)
require(ape)
require(diversitree)
require(phangorn)

# Get the true rates from the ClaDS2 model
tree_files <- list.files("MalietEtAlRuns/MalietEtAl_ClaDS2_trees/")
tree_files <- tree_files[grep("Rdata",tree_files)]

# Function to get true rates from ClaDS2 trees
get.true.ratesClaDS2 <- function(path="/MalietEtAlRuns/MalietEtAl_ClaDS2_trees/") {
  true_clads_tiprates <- list()
  for(i in 1:length(tree_files)) {
    load(paste0(getwd(),path, tree_files[i]))
    if(exists("tree")){
      label <- gsub(".Rdata","",tree_files[i])
      n_tip <- Ntip(tree)
      tmp_tiprates <- as.data.frame(matrix(nrow=n_tip, ncol=7))
      tip_speciation_rates <- c()
      tip_extinction_rates <- c()
      for(tip_index in 1:n_tip){
        tip_speciation_rates[tip_index] <- speciation_rates[which(tree$edge[,2]==tip_index)]
        tip_extinction_rates[tip_index] <- extinction_rates[which(tree$edge[,2]==tip_index)]
      }
      tmp_tiprates[,1] <- label
      tmp_tiprates[,2] <- tree$tip.label
      #library(RPANDA)
      #plot_ClaDS_phylo(tree, speciation_rates) # have to make sure the order is right
      tmp_tiprates[,3] <- tip_speciation_rates
      tmp_tiprates[,4] <- tip_extinction_rates
      tmp_tiprates[,5] <- tip_speciation_rates + tip_extinction_rates 
      tmp_tiprates[,6] <- tip_speciation_rates - tip_extinction_rates 
      tmp_tiprates[,7] <- tip_extinction_rates / tip_speciation_rates 
      true_clads_tiprates[[i]] <- tmp_tiprates
      rm("tree")
    }
  }
  true_clads_tiprates <- do.call(rbind, true_clads_tiprates)
  colnames(true_clads_tiprates) <- c("treeName","tipName","lambdaTRUE","muTRUE",
                                     "netDivTRUE","turnoverTRUE","extinctionFractionTRUE")
  return(true_clads_tiprates)
}


# Functions to summarize rates from MiSSE runs
treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}

true_clads_tiprates <- get.true.ratesClaDS2(path="/MalietEtAlRuns/MalietEtAl_ClaDS2_trees/")
  
all_results <- data.frame()
model_average_results <- data.frame()
best_results <- data.frame()
dones <- c(list.files("MalietEtAlRuns/MiSSE_files/misse_results", pattern="done_recon", full.names=TRUE))

for (i in seq_along(dones)) {
  print(paste0("Loading ", i, " of ", length(dones)))
  try(load(paste0(getwd(),"/", dones[i])))
  if(exists("summary_df")) {
    
    all_results <- plyr::rbind.fill(all_results, summary_df)
    local_weighted <- summary_df %>% group_by(treeName, taxon_id_in_phy, tipName) %>% summarise(
      turnover = weighted.mean(turnover, AICc_weight),
      extinction.fraction = weighted.mean(extinction.fraction, AICc_weight),
      net.div = weighted.mean(net.div, AICc_weight), 
      speciation = weighted.mean(speciation, AICc_weight),
      extinction = weighted.mean(extinction, AICc_weight)
    )
    model_average_results <- plyr::rbind.fill(model_average_results, local_weighted)
    best_results <- plyr::rbind.fill(best_results, subset(summary_df, deltaAICc==0))
    
    rm(summary_df)
    rm(delta_AICc)
    rm(local_weighted)
  }
}


rates.ours <- treeTipMerge(all_results)
rates.ours.model.average <- treeTipMerge(model_average_results)
rates.ours$unique_string <- paste(rates.ours$treeName, "nturnover", rates.ours$nturnover, "neps", rates.ours$neps, "root", rates.ours$root_type, sep="_")
rates.ours.best <- best_results

colnames(rates.ours.best)[which(colnames(rates.ours.best)=="turnover")] <- "turnoverMiSSEbest"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="extinction.fraction")] <- "extinctionFractionMiSSEbest"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="net.div")] <- "netDivMiSSEbest"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="speciation")] <- "lambdaMiSSEbest"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="extinction")] <- "muMiSSEbest"

colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="turnover")] <- "turnoverMiSSEavg"
colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="extinction.fraction")] <- "extinctionFractionMiSSEavg"
colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="net.div")] <- "netDivMiSSEavg"
colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="speciation")] <- "lambdaMiSSEavg"
colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="extinction")] <- "muMiSSEavg"

rates.combined <- merge(rates.ours.best, rates.ours.model.average)
rates.combined <- merge(rates.combined, true_clads_tiprates)

# Retriving ClaDS results:
clads.results <- list.files("MalietEtAlRuns/ClaDS_files/clads", full.names = T)

clads_tip_rates <- list()
for(result_index in 1:length(clads.results)) {
  load(clads.results[result_index])
  tmp_tiprates <- data.frame(CladsOutput$tree$tip.label, CladsOutput$lambdatip_map, CladsOutput$eps_map)
  colnames(tmp_tiprates) <- c("tipName","lambdaClaDS", "extinctionFractionClaDS")
  tmp_tiprates$muClaDS <- tmp_tiprates$lambdaClaDS * tmp_tiprates$extinctionFractionClaDS  
  tmp_tiprates$turnoverClaDS <- tmp_tiprates$lambdaClaDS + tmp_tiprates$muClaDS
  tmp_tiprates$netDivClaDS <- tmp_tiprates$lambdaClaDS - tmp_tiprates$muClaDS
  label <- gsub("_clads_results.Rdata","", clads.results[result_index])
  label <- gsub("MalietEtAlRuns/ClaDS_files/clads/","", label)
  tmp_tiprates$treeName <- label
  clads_tip_rates[[result_index]] <- tmp_tiprates
  print(result_index)
}
all_clads_tip_rates <- do.call(rbind, clads_tip_rates)

#subset_results <- subset(rates.cleaned, rates.cleaned$treeName %in% tree_labels)
# tree_labels <- gsub("_clads_results.Rdata","", clads.results)
# subset_results <- subset(rates.combined, rates.combined$treeName %in% tree_labels)
subset_results <- merge(rates.combined, all_clads_tip_rates, by=c("treeName","tipName"))

# Get BAMM rates:
# setwd("~/Desktop/misse_mme_paper/missecomparison/ClaDScomparison/BAMM_ClaDS_comparison/control_files")

tree.dir <-  "MalietEtAlRuns/MalietEtAl_ClaDS2_trees/"
result.dir <- "MalietEtAlRuns/BAMM_files/all_results/"
tree_files <- list.files(tree.dir)
tree.files <- tree_files[grep(".tre$", tree_files)]

results_files  <- list.files(result.dir)
results_files  <- results_files[grep("_event_data.txt", results_files)]

all_results <- list()
for(tree_index in 1:length(tree.files)){
  label <- gsub(".tre$","", tree.files[tree_index])
  one_result <- results_files[grep(label, results_files)]
  if(length(one_result)!=0) {
    tree_bamm <- read.tree(paste0(tree.dir, tree.files[tree_index]))
    event_data <- getEventData(tree_bamm, paste0(result.dir, one_result), nsamples = 9000) # increase when running final version
    bamm_tiprates <- getTipRates(event_data)    
    # Getting BAMM rates
    lambda <- bamm_tiprates$lambda.avg
    mu <- bamm_tiprates$mu.avg
    bamm_tiprates_final <- data.frame(treeName = label, species=names(lambda), lambdaBAMM=unname(lambda), muBAMM=unname(mu))
    colnames(bamm_tiprates_final) <- c("treeName","tipName","lambdaBAMM","muBAMM")
    bamm_tiprates_final$netDivBAMM <- bamm_tiprates_final$lambdaBAMM - bamm_tiprates_final$muBAMM
    bamm_tiprates_final$turnoverBAMM <- bamm_tiprates_final$lambdaBAMM + bamm_tiprates_final$muBAMM
    bamm_tiprates_final$extinctionFractionBAMM <- bamm_tiprates_final$muBAMM / bamm_tiprates_final$lambdaBAMM 
    all_results[[tree_index]] <- bamm_tiprates_final
  }
}

all_bamm_results <- do.call(rbind, all_results)
subset_results <- merge(subset_results, all_bamm_results, by=c("treeName","tipName"))




# Got function from Title&Rabosky 2019
# DR metric / inverse equal splits
DRstat <- function(tree) {
  
  spRate <- function(sp, tree) {
    #get branch lengths from root to tip
    edges <- vector()
    daughterNode <- match(sp, tree$tip.label)
    while (daughterNode != (length(tree$tip.label) + 1)) {
      parentNode <- tree$edge[which(tree$edge[,2] == daughterNode), 1]
      edges <- c(edges, tree$edge.length[which(tree$edge[,1] == parentNode & tree$edge[,2] == daughterNode)])
      daughterNode <- parentNode
    }
    
    res <- sum(sapply(1:length(edges), function(x) edges[x] * (1/(2 ^ (x-1)))))
    res <- res ^ (-1)
    
    return(res)
  }
  
  rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
  names(rates) <- tree$tip.label
  
  return(rates)
}




# ratio of number of speciation events / age of clade
nodeDensity <- function(tree) {
  maxBT <- max(branching.times(tree))
  nodeCounts <- sapply(1:length(tree$tip.label), function(x) {
    n <- 0
    childnode <- x
    parentNode <- tree$edge[which(tree$edge[,2] == childnode), 1]
    while(parentNode != (length(tree$tip.label) + 1)) {
      n <- n + 1
      childnode <- parentNode
      parentNode <- tree$edge[which(tree$edge[,2] == childnode), 1]
    }
    return(n)
  })
  
  # avoid rates of 0 by adding the root node
  nodeCounts <- nodeCounts + 1
  
  return(setNames(nodeCounts / maxBT, tree$tip.label))
}




# inverse of terminal branch lengths
inverseTerminalBranches <- function(tree) {
  tb <- sapply(1:length(tree$tip.label), function(x) {
    tree$edge.length[which(tree$edge[,2] == x)]	
  })
  return(setNames((1 / (2*tb)), tree$tip.label))	
}




# Fit constant-rate birth-death process function
#    note bounds on lambda optimization; may want to increase these or check for boundary fails
fitCRBD <- function(phy, nopt=5, lmin=0.001, lmax=5.0, MAXBAD = 200){
  
  if (length(phy$tip.label) < 3){
    pars <- c(0.0001,0)
    names(pars) <- c("lambda", "mu")
    return(pars)
  }
  
  fx <- make.bd(phy)
  
  for (i in 1:nopt){
    
    lam <- runif(1, 0, 0.5)	
    mu <- lam * runif(1, 0, 1)
    
    badcount <- 0
    
    resx <- try(optim(c(lam, mu) ,fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T)
    while (class(resx) == 'try-error'){
      
      lam <- runif(1, 0, 0.5)	
      mu <- lam * runif(1, 0, 1)
      
      resx <- try(optim(c(lam, mu) , fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T);
      
      badcount <- badcount + 1;
      if (badcount > MAXBAD){
        stop("Too many fails in fitDiversitree\n")
      }
    }
    
    if (i == 1){
      best <- resx
    }else{
      if (best$value < resx$value){
        best <- resx
      }
    }
    
  }
  
  fres <- list(pars=best$par, loglik=best$value)
  fres$AIC <- -2*fres$loglik + 2*length(argnames(fx))
  fres$counts <- best$counts
  #fres$like_function <- fx
  fres$convergence <- best$convergence
  fres$message <- best$message
  
  pars <- fres$pars
  names(pars) <- c("lambda", "mu")
  
  return(pars)
}


tree_files <- list.files("MalietEtAlRuns/MalietEtAl_ClaDS2_trees", full.names = T)
tree_files <- tree_files[grep("tre$", tree_files)]
all_trees <- lapply(tree_files, read.tree)
names(all_trees) <- gsub(paste0(c("MalietEtAlRuns/MalietEtAl_ClaDS2_trees/",".tre$"), collapse="|"),"", tree_files)

DR_results <- lapply(all_trees, DRstat)
ND_results <- lapply(all_trees, nodeDensity)
TB_results <- lapply(all_trees, inverseTerminalBranches)

# colnames: 
# "treeName" "tipName" "lambdaDR" "lambdaND"  "lambdaTB"

summarize_results <- function(results, type=c("DR","ND","TB")) {
  all_results_tmp <- list()
  for (i in 1:length(results)) {
    one_result <- results[[i]]
    all_results <- data.frame(matrix(nrow=length(one_result), ncol=3))
    all_results[,1] <- names(results)[i]
    all_results[,2] <- names(one_result)
    all_results[,3] <- one_result
    colnames(all_results) <- c("treeName","tipName", paste0("lambda",type))
    all_results_tmp[[i]] <- all_results
  }
  return(do.call(rbind, all_results_tmp))
}

DR <- summarize_results(DR_results, type="DR")
ND <- summarize_results(ND_results, type="ND")
TB <- summarize_results(TB_results, type="TB")

partial_results <- merge(DR, ND, by=c("treeName","tipName"))
nonparametric_results <- merge(partial_results, TB, by=c("treeName","tipName"))


final_results <- merge(subset_results, nonparametric_results, by=c("treeName","tipName"))
write.csv(final_results, file=paste0(getwd(),"/", Sys.Date(), "_preliminary_results_malietetal.csv"), row.names=F)


