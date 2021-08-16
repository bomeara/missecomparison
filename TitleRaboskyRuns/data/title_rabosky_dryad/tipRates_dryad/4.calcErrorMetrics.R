## This script generates rateErrorMetrics.csv, which is already present in the Dryad repository.

# read in true and estimated rates, and calculate error metrics

## Read in datasets
trueRates <- read.csv('dataFiles/trueTipRates.csv', stringsAsFactors=FALSE)
estimatedRates <- read.csv('dataFiles/estimatedTipRates.csv', stringsAsFactors=FALSE)
treeDat <- read.csv('dataFiles/treeSummary.csv', stringsAsFactors=FALSE)
crbd <- read.csv('dataFiles/estimatedCRBD.csv', stringsAsFactors=FALSE)

outfile <- 'dataFiles/rateErrorMetrics.csv'

options(warn=1)
# ----------------------------------------------------
# Combine true and estimated datasets into one object

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

# this object is convenient - save for later use
saveRDS(allRatesList, 'dataFiles/allRatesList.rds', compress = 'xz')

### -----------------------------

# error metric definitions

rmse <- function(estimated, true) {
	sqrt(mean((estimated - true) ^ 2))
}

absoluteError <- function(estimated, true) {
	mean(abs(estimated - true))
}

# PROPORTIONAL ERROR METRICS (from Rabosky et al. 2014 appendix)

# estimated = vector of estimated tip values
# true = true tip values

# overall mean proportional error
propError1 <- function(estimated, true) {
	sum((estimated - true) / true) / length(true)
}

# overall mean proportional error
propError1abs <- function(estimated, true) {
	sum(abs(estimated - true) / true) / length(true)
}

# overall mean proportional error, weighted by the relative size of the regime
# regimeVec is a categorical vector of regimeID of the same length as vec estimated and vec true.
# it should represent the relative size of the regime, to the whole tree
propError2 <- function(estimated, true, regimeVec) {
	regimeVecProp <- setNames(as.numeric(table(regimeVec) / length(regimeVec)), names(table(regimeVec)))
	regimeVec <- regimeVecProp[regimeVec]
	sum(((estimated - true) / true) * regimeVec)
}

# mean error over all tips, using log-transformed rates
propError3 <- function(estimated, true) {
	exp(sum(log(estimated) - log(true)) / length(true))
}

# mean error over all tips, using log-transformed rates, where rates are weighted by regime size.
propError4 <- function(estimated, true, regimeVec) {
	regimeVecProp <- setNames(as.numeric(table(regimeVec) / length(regimeVec)), names(table(regimeVec)))
	regimeVec <- regimeVecProp[regimeVec]
	exp(sum((log(estimated) - log(true)) * regimeVec))
}


# ols regression R2 and slope
olsRegression <- function(estimated, true) {
	if (!all(is.na(estimated))) {
		ols <- lm(estimated ~ true)
		if (nrow(summary(ols)$coefficients) > 1) {
			setNames(c(summary(ols)$coefficients[2,1], summary(ols)$adj.r.squared), c('slope','r2'))
		} else {
			setNames(c(NA, NA), c('slope','r2'))
		}
	} else {
		setNames(c(NA, NA), c('slope','r2'))
	}
}

# -----------------------------------------------

# calculate different error metrics between true and estimated rates for each tree

dat <- as.data.frame(matrix(nrow = nrow(crbd), ncol = 6))
colnames(dat) <- c('treeName', 'setname', 'nTips', 'nRegimes', 'processType')

dat[, c('treeName','setname', 'CRBD_lambda', 'CRBD_mu')] <- crbd[,c('treeName','setname','lambda','mu')]

for (i in 1:nrow(dat)) {
	
	tmp <- treeDat[which(treeDat$treeName == dat[i, 'treeName']),]
	
	dat[i, 'nTips'] <- tmp[1, 'nTips_all']
	dat[i, 'nRegimes'] <- nrow(tmp)
	
	if (grepl('Rabosky2014', tmp[1, 'setname'])) {
		dat[i, 'processType'] <- 'diversity-dependent'
	} else if (tmp[1, 'setname'] == 'evolvingRates') {
		dat[i, 'processType'] <- 'evolvingRates'
	} else {
		dat[i, 'processType'] <- 'birth-death'
	}
}

# check ordering
identical(names(allRatesList), dat$treeName)

# calculate proportional error of estimated to true lambda tip rates for each different tip rate metric

tipMetricHeaders <- c('CRBDlambda','lambdaTB','lambdaND','lambdaDR','lambdaBAMM')
tipMetricHeadersNetDiv <- c('CRBDlambda','lambdaTB','lambdaND','lambdaDR','netDivBAMM')

for (i in 1:nrow(dat)) {
	
	cat(i, '\n')
	
	for (j in 1:length(tipMetricHeaders)) {
		
		# speciation rate
		dat[i, paste0('tipLambda_', tipMetricHeaders[j], '_PE1')] <- propError1(allRatesList[[i]][, tipMetricHeaders[j]], allRatesList[[i]]$tipLambda)

		dat[i, paste0('tipLambda_', tipMetricHeaders[j], '_PE2')] <- propError2(allRatesList[[i]][, tipMetricHeaders[j]], allRatesList[[i]]$tipLambda, allRatesList[[i]]$regimeID)

		dat[i, paste0('tipLambda_', tipMetricHeaders[j], '_PE3')] <- propError3(allRatesList[[i]][, tipMetricHeaders[j]], allRatesList[[i]]$tipLambda)

		dat[i, paste0('tipLambda_', tipMetricHeaders[j], '_PE4')] <- propError4(allRatesList[[i]][, tipMetricHeaders[j]], allRatesList[[i]]$tipLambda, allRatesList[[i]]$regimeID)

		dat[i, paste0('tipLambda_', tipMetricHeaders[j], '_slope')] <- olsRegression(allRatesList[[i]][, tipMetricHeaders[j]], allRatesList[[i]]$tipLambda)[1]

		dat[i, paste0('tipLambda_', tipMetricHeaders[j], '_r2')] <- olsRegression(allRatesList[[i]][, tipMetricHeaders[j]], allRatesList[[i]]$tipLambda)[2]

		dat[i, paste0('tipLambda_', tipMetricHeaders[j], '_absoluteError')] <- absoluteError(allRatesList[[i]][, tipMetricHeaders[j]], allRatesList[[i]]$tipLambda)

		dat[i, paste0('tipLambda_', tipMetricHeaders[j], '_rmse')] <- rmse(allRatesList[[i]][, tipMetricHeaders[j]], allRatesList[[i]]$tipLambda)
		

		# Net diversification rate
		dat[i, paste0('tipNetDiv_', tipMetricHeadersNetDiv[j], '_PE1')] <- propError1(allRatesList[[i]][, tipMetricHeadersNetDiv[j]], allRatesList[[i]]$tipNetDiv)

		dat[i, paste0('tipNetDiv_', tipMetricHeadersNetDiv[j], '_PE2')] <- propError2(allRatesList[[i]][, tipMetricHeadersNetDiv[j]], allRatesList[[i]]$tipNetDiv, allRatesList[[i]]$regimeID)

		dat[i, paste0('tipNetDiv_', tipMetricHeadersNetDiv[j], '_PE3')] <- propError3(allRatesList[[i]][, tipMetricHeadersNetDiv[j]], allRatesList[[i]]$tipNetDiv)

		dat[i, paste0('tipNetDiv_', tipMetricHeadersNetDiv[j], '_PE4')] <- propError4(allRatesList[[i]][, tipMetricHeadersNetDiv[j]], allRatesList[[i]]$tipNetDiv, allRatesList[[i]]$regimeID)

		dat[i, paste0('tipNetDiv_', tipMetricHeadersNetDiv[j], '_slope')] <- olsRegression(allRatesList[[i]][, tipMetricHeadersNetDiv[j]], allRatesList[[i]]$tipNetDiv)[1]

		dat[i, paste0('tipNetDiv_', tipMetricHeadersNetDiv[j], '_r2')] <- olsRegression(allRatesList[[i]][, tipMetricHeadersNetDiv[j]], allRatesList[[i]]$tipNetDiv)[2]

		dat[i, paste0('tipNetDiv_', tipMetricHeadersNetDiv[j], '_absoluteError')] <- absoluteError(allRatesList[[i]][, tipMetricHeadersNetDiv[j]], allRatesList[[i]]$tipNetDiv)

		dat[i, paste0('tipNetDiv_', tipMetricHeadersNetDiv[j], '_rmse')] <- rmse(allRatesList[[i]][, tipMetricHeadersNetDiv[j]], allRatesList[[i]]$tipNetDiv)


	}
		
}

head(dat)

write.csv(dat, outfile, row.names=FALSE)
