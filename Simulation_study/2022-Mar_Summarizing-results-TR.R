############################################################
# Code to summarize results from simulation studies
############################################################
# rm(list=ls())
############################################################
# Title&Rabosky results (original results + MiSSE and ClaDS)
############################################################

# setwd("~/Desktop/misse_mme_paper/missecomparison/Simulation_study/")
library(Metrics)
library(tidyverse)
library(magrittr)
library(stats)
library(ggplot2)
library(progress)
library(viridis)

treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}

all_results <- data.frame()
model_average_results <- data.frame()
best_results <- data.frame()
dones <- c(list.files("TitleRaboskyRuns/results", pattern="done_recon", full.names=TRUE), list.files("TitleRaboskyRuns/new_results", pattern="done_recon", full.names=TRUE))

for (i in seq_along(dones)) {
  print(paste0("Loading ", i, " of ", length(dones)))
  try(load(dones[i]))
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
rates.theirs <- treeTipMerge(read.csv("TitleRaboskyRuns/original_data/title_rabosky_dryad/tipRates_dryad/dataFiles/estimatedTipRates.csv", stringsAsFactors=FALSE))
rates.true <- treeTipMerge(read.csv("TitleRaboskyRuns/original_data/title_rabosky_dryad/tipRates_dryad/dataFiles/trueTipRates.csv", stringsAsFactors=FALSE))
rates.ours$unique_string <- paste(rates.ours$treeName, "nturnover", rates.ours$nturnover, "neps", rates.ours$neps, "root", rates.ours$root_type, sep="_")
rates.ours.best <- best_results

rates.theirs <- rates.theirs[rates.theirs$treeName %in% rates.ours.best$treeName,]
rates.true <- rates.true[rates.true$treeName %in% rates.ours.best$treeName,]

colnames(rates.true)[which(colnames(rates.true)=="tipLambda")] <- "lambdaTRUE"
colnames(rates.true)[which(colnames(rates.true)=="tipMu")] <- "muTRUE"
colnames(rates.true)[which(colnames(rates.true)=="tipNetDiv")] <- "netDivTRUE"
rates.true$turnoverTRUE <- rates.true$lambdaTRUE + rates.true$muTRUE
rates.true$extinctionFractionTRUE <- rates.true$muTRUE / rates.true$lambdaTRUE

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
rates.combined <- merge(rates.combined, rates.theirs)
rates.combined <- merge(rates.combined, rates.true)

rates.combined$muBAMM <- rates.combined$lambdaBAMM - rates.combined$netDivBAMM
rates.combined$turnoverBAMM <- rates.combined$muBAMM + rates.combined$lambdaBAMM
rates.combined$extinctionFractionBAMM <- rates.combined$muBAMM / rates.combined$lambdaBAMM

# get rid of NAs and weird values for now
rates.cleaned <- rates.combined
rates.cleaned <- rates.cleaned[!is.na(rates.cleaned$netDivBAMM),]
rates.cleaned <- rates.cleaned[(rates.cleaned$lambdaTRUE>0),] #yeah, not sure why there'd be no speciation in reality for a tree with >2 taxa

#unique(rates.cleaned$treeName)

############################################################
# Include ClaDS results
############################################################
clads.dir <- "/Volumes/NO NAME/ClaDS_results" # they are very heavy
clads.results <- list.files(clads.dir)

clads_tip_rates <- list()
for(result_index in 1:length(clads.results)) {
  load(paste0(clads.dir, "/", clads.results[result_index]))
  tmp_tiprates <- data.frame(CladsOutput$tree$tip.label, CladsOutput$lambdatip_map, CladsOutput$eps_map)
  colnames(tmp_tiprates) <- c("tipName","lambdaClaDS", "extinctionFractionClaDS")
  tmp_tiprates$muClaDS <- tmp_tiprates$lambdaClaDS * tmp_tiprates$extinctionFractionClaDS  
  tmp_tiprates$turnoverClaDS <- tmp_tiprates$lambdaClaDS + tmp_tiprates$muClaDS
  tmp_tiprates$netDivClaDS <- tmp_tiprates$lambdaClaDS - tmp_tiprates$muClaDS
  tmp_tiprates$treeName <- gsub("_clads_results.Rdata","", clads.results[result_index])
  clads_tip_rates[[result_index]] <- tmp_tiprates
  print(result_index)
}
all_clads_tip_rates <- do.call(rbind, clads_tip_rates)

############################################################
# Merge all Title&Rabosky results
############################################################
tree_labels <- gsub("_clads_results.Rdata","", clads.results)
subset_results <- subset(rates.cleaned, rates.cleaned$treeName %in% tree_labels)
subset_results <- merge(subset_results, all_clads_tip_rates,  by=c("treeName","tipName"))

write.csv(subset_results, file=paste0(getwd(), Sys.Date(), "_preliminary_results_titlerabosky.csv"), row.names=F)


