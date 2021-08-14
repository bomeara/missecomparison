# Summarize results by true regime
## For multi-regime trees, summarize tip rates (true and inferred) by known regimes

outfile <- 'dataFiles/regimeSummary.csv'

trueRates <- read.csv('dataFiles/trueTipRates.csv', stringsAsFactors=FALSE)
estimatedRates <- read.csv('dataFiles/estimatedTipRates.csv', stringsAsFactors=FALSE)
treeDat <- read.csv('dataFiles/treeSummary.csv', stringsAsFactors=FALSE)
crbd <- read.csv('dataFiles/estimatedCRBD.csv', stringsAsFactors=FALSE)

trueRatesList <- split(trueRates, trueRates$treeName)
estimatedRatesList <- split(estimatedRates, estimatedRates$treeName)

# make sure lists are in the same order, and tips are in the same order in each list element
identical(names(trueRatesList), names(estimatedRatesList))

for (i in 1:length(trueRatesList)) {
	xx <- trueRatesList[[i]]$tipName
	xx <- sort(xx)
	trueRatesList[[i]] <- trueRatesList[[i]][match(xx, trueRatesList[[i]]$tipName),]
	estimatedRatesList[[i]] <- estimatedRatesList[[i]][match(xx, estimatedRatesList[[i]]$tipName),]
}

# put CRBD table in same order
ordering <- sapply(names(trueRatesList), function(x) which(crbd$treeName == x))
crbd <- crbd[ordering,]

# merge the two, and add CRBD lambda and mu (identical values across tip names)
allRatesList <- vector('list', length = length(trueRatesList))
names(allRatesList) <- names(trueRatesList)
for (i in 1:length(trueRatesList)) {
	
	allRatesList[[i]] <- cbind.data.frame(trueRatesList[[i]], estimatedRatesList[[i]][, c('lambdaTB', 'lambdaND', 'lambdaDR', 'lambdaBAMM', 'netDivBAMM')], CRBDlambda = crbd[i, 'lambda'], CRBDmu = crbd[i, 'mu'])
	
}


# use regime ID to create table of regimes, and true/inferred rate averages
# single-regime trees will be a single regime
# also record size of regimes

setvec <- c('MitchellRabosky2016','MeyerWiens2017','fossilBAMM','MooreEtAl2016','RaboskyEtAl2017','lambdaConstantVariance', 'netDivConstantVariance')

setInd <- unique(treeDat[which(treeDat$setname %in% setvec), 'treeName'])

regimeList <- list()

counter <- 1
for (i in 1:length(setInd)) {
	
	cat(i, 'of', length(setInd), '\n')
	xx <- allRatesList[[setInd[i]]]
	
	# take average of all tip metrics, grouped by regimeID
	tmp <- aggregate(xx[, colnames(xx)[which(colnames(xx) == 'tipLambda'):ncol(xx)]], list(xx$regimeID), mean)
	regimeList[[counter]] <- cbind(treeName = xx[1, 'treeName'], setname = xx[1, 'setname'], nTips = aggregate(xx$regimeID, list(xx$regimeID), length)$x, tmp)
		
	counter <- counter + 1	
	
}

regimeSummary <- do.call(rbind, regimeList)
regimeSummary <- regimeSummary[, -which(colnames(regimeSummary) == 'Group.1')]
regimeSummary <- as.data.frame(regimeSummary, stringsAsFactors=FALSE)

for (i in 3:ncol(regimeSummary)) {
	regimeSummary[, i] <- as.numeric(regimeSummary[, i])
}
regimeSummary[,1] <- as.character(regimeSummary[,1])
regimeSummary[,2] <- as.character(regimeSummary[,2])

write.csv(regimeSummary, outfile, row.names = FALSE)






## Diversity-dependent trees
outfile <- 'dataFiles/regimeSummary_dd.csv'

setvec <- c('Rabosky2014_DD_k1', 'Rabosky2014_DD_k2', 'Rabosky2014_DD_k3', 'Rabosky2014_DD_k4')

setInd <- unique(treeDat[which(treeDat$setname %in% setvec), 'treeName'])

regimeList <- list()

counter <- 1
for (i in 1:length(setInd)) {
	
	cat(i, 'of', length(setInd), '\n')
	xx <- allRatesList[[setInd[i]]]
	
	# take average of all tip metrics, grouped by regimeID
	tmp <- aggregate(xx[, colnames(xx)[which(colnames(xx) == 'tipLambda'):ncol(xx)]], list(xx$regimeID), mean)
	regimeList[[counter]] <- cbind(treeName = xx[1, 'treeName'], setname = xx[1, 'setname'], nTips = aggregate(xx$regimeID, list(xx$regimeID), length)$x, tmp)
		
	counter <- counter + 1	
	
}

regimeSummary <- do.call(rbind, regimeList)
regimeSummary <- regimeSummary[, -which(colnames(regimeSummary) == 'Group.1')]
regimeSummary <- as.data.frame(regimeSummary, stringsAsFactors=FALSE)

for (i in 3:ncol(regimeSummary)) {
	regimeSummary[, i] <- as.numeric(regimeSummary[, i])
}
regimeSummary[,1] <- as.character(regimeSummary[,1])
regimeSummary[,2] <- as.character(regimeSummary[,2])

write.csv(regimeSummary, outfile, row.names = FALSE)




## Evolving Rates trees

outfile <- 'dataFiles/regimeSummary_evolvingRates.csv'

setvec <- c('evolvingRates')

setInd <- unique(treeDat[which(treeDat$setname %in% setvec), 'treeName'])

regimeList <- list()

counter <- 1
for (i in 1:length(setInd)) {
	
	cat(i, 'of', length(setInd), '\n')
	xx <- allRatesList[[setInd[i]]]
	
	# take average of all tip metrics, grouped by regimeID
	tmp <- aggregate(xx[, colnames(xx)[which(colnames(xx) == 'tipLambda'):ncol(xx)]], list(xx$regimeID), mean)
	regimeList[[counter]] <- cbind(treeName = xx[1, 'treeName'], setname = xx[1, 'setname'], nTips = aggregate(xx$regimeID, list(xx$regimeID), length)$x, tmp)
		
	counter <- counter + 1	
	
}

regimeSummary <- do.call(rbind, regimeList)
regimeSummary <- regimeSummary[, -which(colnames(regimeSummary) == 'Group.1')]
regimeSummary <- as.data.frame(regimeSummary, stringsAsFactors=FALSE)

for (i in 3:ncol(regimeSummary)) {
	regimeSummary[, i] <- as.numeric(regimeSummary[, i])
}
regimeSummary[,1] <- as.character(regimeSummary[,1])
regimeSummary[,2] <- as.character(regimeSummary[,2])

write.csv(regimeSummary, outfile, row.names = FALSE)

