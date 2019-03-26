rm(list=ls())
library(Metrics)

treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}

rates.ours <- treeTipMerge(read.csv("result.csv", stringsAsFactors=FALSE))
rates.theirs <- treeTipMerge(read.csv("data/pascal_rabosky_dryad/tipRates_dryad/dataFiles/estimatedTipRates.csv", stringsAsFactors=FALSE))
rates.true <- treeTipMerge(read.csv("data/pascal_rabosky_dryad/tipRates_dryad/dataFiles/trueTipRates.csv", stringsAsFactors=FALSE))

rates.theirs <- rates.theirs[rates.theirs$treeName %in% rates.ours$treeName,]
rates.true <- rates.true[rates.true$treeName %in% rates.ours$treeName,]

colnames(rates.true)[which(colnames(rates.true)=="tipLambda")] <- "tipLambdaTRUE"
colnames(rates.true)[which(colnames(rates.true)=="tipMu")] <- "tipMuTRUE"
colnames(rates.true)[which(colnames(rates.true)=="tipNetDiv")] <- "tipNetDivTRUE"
rates.true$tipTurnoverTRUE <- rates.true$tipLambdaTRUE + rates.true$tipMuTRUE
rates.true$tipExtinctionFractionTrue <- rates.true$tipMuTRUE / rates.true$TipLambdaTRUE

colnames(rates.ours)[which(colnames(rates.ours)=="turnover")] <- "tipTurnoverMiSSE"
colnames(rates.ours)[which(colnames(rates.ours)=="extinction.fraction")] <- "tipExtinctionFractionMiSSE"
colnames(rates.ours)[which(colnames(rates.ours)=="net.div")] <- "tipNetDivMiSSE"
colnames(rates.ours)[which(colnames(rates.ours)=="speciation")] <- "tipLambdaMiSSE"
colnames(rates.ours)[which(colnames(rates.ours)=="extinction")] <- "tipMuMiSSE"

rates.combined <- merge(rates.ours, rates.theirs, by="treeTipString")
rates.combined <- merge(rates.combined, rates.true, by="treeTipString")

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
