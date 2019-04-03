rm(list=ls())
library(Metrics)

treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}

rates.ours <- treeTipMerge(read.csv("result.csv", stringsAsFactors=FALSE))
rates.theirs <- treeTipMerge(read.csv("data/title_rabosky_dryad/tipRates_dryad/dataFiles/estimatedTipRates.csv", stringsAsFactors=FALSE))
rates.true <- treeTipMerge(read.csv("data/title_rabosky_dryad/tipRates_dryad/dataFiles/trueTipRates.csv", stringsAsFactors=FALSE))

rates.ours$unique_string <- paste(rates.ours$treeName, "nturnover", rates.ours$nturnover, "neps", rates.ours$neps, "root", rates.ours$root_type, sep="_")

rates.ours.best <- data.frame()
unique.trees <- unique(rates.ours$treeName)
for (i in seq_along(unique.trees)) {
  rates.local.df <- rates.ours[which(rates.ours$treeName == unique.trees[i]),]
  rates.local <- split(rates.local.df, rates.local.df$unique_string)
  AICc <- sapply(lapply(rates.local, "[[", "AICc"),min)
  deltaAICcToNext <- NA
  if(length(rates.local)>1) {
    AICc.sorted <- sort(AICc, decreasing=FALSE)
    deltaAICcToNext <- AICc.sorted[2] - AICc.sorted[1]
  }
  rates.local.best <- rates.local[[which.min(AICc)[1]]]
  rates.local.best$numberAlternatives <- length(rates.local)
  rates.local.best$deltaAICcToNext <- deltaAICcToNext
  rates.ours.best <- rbind(rates.ours.best, rates.local.best)
}

rates.theirs <- rates.theirs[rates.theirs$treeName %in% rates.ours.best$treeName,]
rates.true <- rates.true[rates.true$treeName %in% rates.ours.best$treeName,]

colnames(rates.true)[which(colnames(rates.true)=="tipLambda")] <- "lambdaTRUE"
colnames(rates.true)[which(colnames(rates.true)=="tipMu")] <- "muTRUE"
colnames(rates.true)[which(colnames(rates.true)=="tipNetDiv")] <- "netDivTRUE"
rates.true$turnoverTRUE <- rates.true$lambdaTRUE + rates.true$muTRUE
rates.true$extinctionFractionTRUE <- rates.true$muTRUE / rates.true$lambdaTRUE

colnames(rates.ours.best)[which(colnames(rates.ours.best)=="turnover")] <- "turnoverMiSSE"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="extinction.fraction")] <- "extinctionFractionMiSSE"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="net.div")] <- "netDivMiSSE"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="speciation")] <- "lambdaMiSSE"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="extinction")] <- "muMiSSE"

rates.combined <- merge(rates.ours.best, rates.theirs, by="treeTipString")
rates.combined <- merge(rates.combined, rates.true, by="treeTipString")

rates.combined$muBAMM <- rates.combined$lambdaBAMM - rates.combined$netDivBAMM
rates.combined$turnoverBAMM <- rates.combined$muBAMM + rates.combined$lambdaBAMM
rates.combined$extinctionFractionBAMM <- rates.combined$muBAMM / rates.combined$lambdaBAMM

approaches <- c("TB", "ND", "DR", "BAMM", "MiSSE")
parameters <- c("mu", "lambda", "netDiv", "turnover", "extinctionFraction")
RMSE.results <- data.frame(matrix(nrow=length(approaches), ncol=length(parameters)))
rownames(RMSE.results) <- approaches
colnames(RMSE.results) <- parameters

absoluteError.mean.results <- RMSE.results
absoluteError.median.results <- RMSE.results


# get rid of NAs and weird values for now
rates.cleaned <- rates.combined
rates.cleaned <- rates.cleaned[!is.na(rates.cleaned$netDivBAMM),]
rates.cleaned <- rates.cleaned[(rates.cleaned$lambdaTRUE>0),] #yeah, not sure why there'd be no speciation in realisty for a tree with >2 taxa


for (approach.index in seq_along(approaches)) {
  for (parameter.index in seq_along(parameters)) {
    truth <- rates.cleaned[,paste0(parameters[parameter.index],"TRUE")]
    estimate_name <- paste0(parameters[parameter.index], approaches[approach.index])
    if(estimate_name %in% colnames(rates.cleaned)) {
      estimate <- rates.cleaned[,estimate_name]
      RMSE.results[approach.index, parameter.index] <- Metrics::rmse(estimate, truth)
      absoluteError.mean.results[approach.index, parameter.index] <- mean(abs(estimate-truth))
      absoluteError.median.results[approach.index, parameter.index] <- median(abs(estimate-truth))

    }
  }
}

print(round(RMSE.results,3))
print(round(absoluteError.mean.results,3))
print(round(absoluteError.median.results,3))
