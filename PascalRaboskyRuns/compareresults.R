rm(list=ls())
library(Metrics)

treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}

rates.ours <- treeTipMerge(read.csv("result.csv", stringsAsFactors=FALSE))
rates.theirs <- treeTipMerge(read.csv("data/pascal_rabosky_dryad/tipRates_dryad/dataFiles/estimatedTipRates.csv", stringsAsFactors=FALSE))
rates.true <- treeTipMerge(read.csv("data/pascal_rabosky_dryad/tipRates_dryad/dataFiles/trueTipRates.csv", stringsAsFactors=FALSE))

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
    }
  }
}


lambda.to.compare <- c("lambdaTB", "lambdaND", "lambdaDR", "lambdaBAMM", "netDivBAMM", "tipLambdaMiSSE",  "tipNetDivMiSSE")
lambdasRMSE.all <- rep(NA, length(lambda.to.compare))
rates.combined.no.na <- rates.combined[!is.na(rates.combined$lambdaBAMM),]
names(lambdasRMSE.all) <- lambda.to.compare
for (i in seq_along(lambda.to.compare)) {
  lambdasRMSE.all[i] <- Metrics::rmse(rates.combined.no.na[,lambda.to.compare[i]], rates.combined.no.na$tipLambdaTRUE)
}

netdivRMSE.all <- rep(NA, length(lambda.to.compare))
rates.combined.no.na <- rates.combined[!is.na(rates.combined$lambdaBAMM),]
names(netdivRMSE.all) <- lambda.to.compare
for (i in seq_along(lambda.to.compare)) {
  netdivRMSE.all[i] <- Metrics::rmse(rates.combined.no.na[,lambda.to.compare[i]], rates.combined.no.na$tipNetDivTRUE)
}

print(t(t(round(netdivRMSE.all,3))))



lambdasRMSE <- rep(NA, length(lambda.to.compare))
names(lambdasRMSE) <- lambda.to.compare
for (i in seq_along(lambda.to.compare)) {
 lambdasRMSE[i] <- Metrics::rmse(rates.combined[,lambda.to.compare[i]], rates.combined$tipLambdaTRUE)
}
