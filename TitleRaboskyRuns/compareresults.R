# setwd("~/Desktop/missecomparison/TitleRaboskyRuns") # TV local
rm(list=ls())
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
dones <- list.files("results", pattern="done", full.names=TRUE)
for (i in seq_along(dones)) {
	print(paste0("Loading ", i, " of ", length(dones)))
	try(load(dones[i]))
	if(exists("summary_df")) {
		
		# summary_df$loglik <- NA
		# summary_df$nparam <- NA
		# good_enough_recon <- which(delta_AICc<20)
		# model_likelihoods <- sapply(hisse_result_all, "[[", "loglik")
		# current_row <- 0
		# for(good_enough_index in seq_along(good_enough_recon)) {
		# 	for (taxon_index in sequence(length(unique(summary_df$taxon_id_in_phy)))) {
		# 		current_row <- current_row+1
		# 		summary_df$loglik[current_row] <- model_likelihoods[good_enough_recon[good_enough_index]]
		# 		summary_df$nparam[current_row] <- 0.5*(summary_df$AIC[current_row] + 2*summary_df$loglik[current_row])
		# 	}
		# }


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
#rates.ours <- treeTipMerge(read.csv("result.csv", stringsAsFactors=FALSE))
rates.ours <- treeTipMerge(all_results)
rates.ours.model.average <- treeTipMerge(model_average_results)

rates.theirs <- treeTipMerge(read.csv("data/title_rabosky_dryad/tipRates_dryad/dataFiles/estimatedTipRates.csv", stringsAsFactors=FALSE))
rates.true <- treeTipMerge(read.csv("data/title_rabosky_dryad/tipRates_dryad/dataFiles/trueTipRates.csv", stringsAsFactors=FALSE))

rates.ours$unique_string <- paste(rates.ours$treeName, "nturnover", rates.ours$nturnover, "neps", rates.ours$neps, "root", rates.ours$root_type, sep="_")

rates.ours.best <- best_results
# rates.ours.best <- data.frame()
# unique.trees <- unique(rates.ours$treeName)
# for (i in seq_along(unique.trees)) {
#   rates.local.df <- rates.ours[which(rates.ours$treeName == unique.trees[i]),]
#   rates.local <- split(rates.local.df, rates.local.df$unique_string)
#   AICc <- sapply(lapply(rates.local, "[[", "AICc"),min)
#   deltaAICcToNext <- NA
#   if(length(rates.local)>1) {
#     AICc.sorted <- sort(AICc, decreasing=FALSE)
#     deltaAICcToNext <- AICc.sorted[2] - AICc.sorted[1]
#   }
#   rates.local.best <- rates.local[[which.min(AICc)[1]]]
#   rates.local.best$numberAlternatives <- length(rates.local)
#   rates.local.best$deltaAICcToNext <- deltaAICcToNext
#   rates.ours.best <- rbind(rates.ours.best, rates.local.best)
# }

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

approaches <- c("TB", "ND", "DR", "BAMM", "MiSSEbest", "MiSSEavg")
parameters <- c("mu", "lambda", "netDiv", "turnover", "extinctionFraction")
RMSE.results <- data.frame(matrix(nrow=length(approaches), ncol=length(parameters)))
rownames(RMSE.results) <- approaches
colnames(RMSE.results) <- parameters

absoluteError.mean.results <- RMSE.results
absoluteError.median.results <- RMSE.results
cv.results <- RMSE.results


# get rid of NAs and weird values for now
rates.cleaned <- rates.combined
rates.cleaned <- rates.cleaned[!is.na(rates.cleaned$netDivBAMM),]
rates.cleaned <- rates.cleaned[(rates.cleaned$lambdaTRUE>0),] #yeah, not sure why there'd be no speciation in reality for a tree with >2 taxa

coef.var <- function(x){
  cv <- sd(x) / mean(x) * 100
  return(cv)
}

for (approach.index in seq_along(approaches)) {
  for (parameter.index in seq_along(parameters)) {
    truth <- rates.cleaned[,paste0(parameters[parameter.index],"TRUE")]
    estimate_name <- paste0(parameters[parameter.index], approaches[approach.index])
    if(estimate_name %in% colnames(rates.cleaned)) {
      estimate <- rates.cleaned[,estimate_name]
      RMSE.results[approach.index, parameter.index] <- Metrics::rmse(estimate, truth)
      absoluteError.mean.results[approach.index, parameter.index] <- mean(abs(estimate-truth))
      absoluteError.median.results[approach.index, parameter.index] <- median(abs(estimate-truth))
      cv.results[approach.index, parameter.index] <- coef.var(abs(estimate-truth))
    }
  }
}

print("RMSE")
print(round(RMSE.results,3))

print("absolute error, mean")
print(round(absoluteError.mean.results,3))

print("absolute error, median")
print(round(absoluteError.median.results,3))

pivot_names <- colnames(rates.combined)[grepl("(lambda)|(mu)|(^turnover)|(extinctionFraction)|(netDiv)", colnames(rates.combined))]
pivot_names <- pivot_names[!grepl("TRUE", pivot_names)]
rates.cleaned.longer <- pivot_longer(rates.cleaned, cols=all_of(pivot_names), names_to="parameter_name", values_to="parameter_value")

rates.cleaned.longer$estimator <- gsub("(lambda)|(mu)|(^turnover)|(extinctionFraction)|(netDiv)", "", rates.cleaned.longer$parameter_name)
rates.cleaned.longer$actual_parameter <- gsub("(MiSSEbest)|(MiSSEavg)|(TB)|(ND)|(DR)|(BAMM)", "", rates.cleaned.longer$parameter_name)

rates.cleaned.longer$true_value <- NA

save(rates.cleaned.longer, file="results/rates.cleaned.longer.rda")

print(dim(rates.cleaned.longer))
unique_params <- unique(rates.cleaned.longer$actual_parameter)
for(i in sequence(length(unique_params))) {
	matching_rows <- which(rates.cleaned.longer$actual_parameter==unique_params[i])
	matching_col <- which(colnames(rates.cleaned.longer)==paste0(unique_params[i], "TRUE"))
	rates.cleaned.longer$true_value[matching_rows] <- as.data.frame(rates.cleaned.longer[matching_rows, matching_col])[,1]
	print(dim(rates.cleaned.longer))
}

rates.cleaned.good <- subset(rates.cleaned.longer, estimator!="ND")
rates.cleaned.good <- subset(rates.cleaned.good, estimator!="TB")
rates.cleaned.good <- subset(rates.cleaned.good, estimator!="DR")
rates.cleaned.good <- subset(rates.cleaned.good, actual_parameter!="extinctionFraction") #wrong scale for this

g <- ggplot(rates.cleaned.good, aes(true_value, parameter_value)) + geom_point(alpha=0.05) + facet_grid(vars(estimator), vars(actual_parameter)) + geom_abline(slope=1, intercept=0) + xlim(min(rates.cleaned.good$true_value), 2*max(rates.cleaned.good$true_value)) + ylim(min(rates.cleaned.good$true_value), 2*max(rates.cleaned.good$true_value))
pdf(file="results/rateplot.pdf", width=20, height=20)
print(g)
dev.off()
system("open results/rateplot.pdf")

unique_params <- unique(rates.cleaned.good$actual_parameter)
for(i in sequence(length(unique_params))) {
	focal <- subset(rates.cleaned.good, actual_parameter==unique_params[i])
	g <- ggplot(focal, aes(true_value, parameter_value)) + facet_wrap(vars(estimator)) + geom_abline(slope=1, intercept=0) + geom_hex(bins=100) + scale_fill_viridis() + theme_bw() 
	pdf(file=paste0("results/rates_", unique_params[i], ".pdf"), width=20, height=7)
	print(g)
	dev.off()
	system(paste0("open results/rates_", unique_params[i], ".pdf"))

}