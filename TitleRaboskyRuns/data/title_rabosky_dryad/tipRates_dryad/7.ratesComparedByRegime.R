# Analyze estimated rates as a function of true rates, when averaged by true regimes

regimeSummary <- read.csv('dataFiles/regimeSummary.csv', stringsAsFactors=FALSE)

regimeSummary[, 'CRBDnetdiv'] <- regimeSummary$CRBDlambda - regimeSummary$CRBDmu
# ---------------------------------------------------------
# ---------------------------------------------------------

# # Examine strength of estimated-to-true rate correlation across regimes, as a function of regime size. 
# # Use only multi-regime time-constant trees
# # Regime size along x-axis and Pearson correlation coefficient and slope as y-axis

regimeSummary2 <- regimeSummary[which(regimeSummary$setname %in% c('MeyerWiens2017','fossilBAMM','MooreEtAl2016','RaboskyEtAl2017')),]

maxRegimeCount <- 50

regimeSize <- sort(unique(regimeSummary2$nTips))

regimeSize <- 1:maxRegimeCount

rateMetrics <- c('CRBDlambda','lambdaTB','lambdaND','lambdaDR','lambdaBAMM')
rateMetricsNetDiv <- c('CRBDnetdiv','lambdaTB','lambdaND','lambdaDR','netDivBAMM')


regimeSizeTable <- as.data.frame(matrix(nrow = length(regimeSize), ncol = 21))
colnames(regimeSizeTable) <- c('regimeSize', paste0(rateMetrics, 'lambda_cor'), paste0(rateMetrics, 'lambda_slope'), paste0(rateMetricsNetDiv, 'netdiv_cor'), paste0(rateMetricsNetDiv, 'netdiv_slope'))
regimeSizeTable$regimeSize <- regimeSize

for (i in 1:length(regimeSize)) {
	
	xx <- regimeSummary2[which(regimeSummary2$nTips >= regimeSize[i]),]
	
	# speciation rate
	if (nrow(xx) > 1) {
		for (j in 1:length(rateMetrics)) {
			regimeSizeTable[i, paste0(rateMetrics[j], 'lambda_cor')] <- cor(xx$tipLambda, xx[, rateMetrics[j]])
			regimeSizeTable[i, paste0(rateMetrics[j], 'lambda_slope')] <- summary(lm(xx[, rateMetrics[j]] ~ xx$tipLambda))$coef[2,1]
		}
	}
	
	# net diversification
	if (nrow(xx) > 1) {
		for (j in 1:length(rateMetrics)) {
			regimeSizeTable[i,paste0(rateMetricsNetDiv[j], 'netdiv_cor')] <- cor(xx$tipNetDiv, xx[, rateMetricsNetDiv[j]])
			regimeSizeTable[i,paste0(rateMetricsNetDiv[j], 'netdiv_slope')] <- summary(lm(xx[, rateMetricsNetDiv[j]] ~ xx$tipNetDiv))$coef[2,1]
		}
	}
	
}



# Speciation rate, Pearson correlation
tipMetrics <- c('CRBDlambda','lambdaTB','lambdaND','lambdaDR','lambdaBAMM')

metricColors <- c('black','dark green','dark orange','red','blue')

pdf('fig4.pdf', width=10, height=5)

par(mfrow=c(1,2), mar=c(4,4,0,0))

# speciation rate, correlation
plot.new()
plot.window(xlim = range(regimeSize), ylim = c(0,1))
axis(1, at = c(-1, axTicks(1)))
axis(2, at = c(-1, axTicks(2)))
mtext('minimum regime size', side = 1, line = 2.5, cex=1)
mtext('Pearson correlation', side = 2, line = 2.5, cex=1)
for (i in 1:length(rateMetrics)) {
	lines(regimeSize, regimeSizeTable[, paste0(rateMetrics[i], 'lambda_cor')], col= metricColors[i], lwd=2)
}
abline(h=1, lty=3)
mtext('(a)', outer = TRUE, line = -1.5, at = 0.02, font=2)

legend(35, 0.55, legend = c(expression(lambda['CRBD']), expression(lambda['TB']), expression(lambda['ND']), expression(lambda['DR']), expression(lambda['BAMM'])), col = metricColors, lwd=2, bty='n', xpd=NA, cex=1.2)

# speciation rate, slope
plot.new()
plot.window(xlim = range(regimeSize), ylim = c(0,1.5))
axis(1, at = c(-1, axTicks(1)))
axis(2, at = c(-1, axTicks(2)))
mtext('minimum regime size', side = 1, line = 2.5, cex=1)
mtext('OLS slope', side = 2, line = 2.5, cex=1)
for (i in 1:length(rateMetrics)) {
	lines(regimeSize, regimeSizeTable[, paste0(rateMetrics[i], 'lambda_slope')], col= metricColors[i], lwd=2)
}
abline(h=1, lty=3)
mtext('(b)', outer = TRUE, line = -1.5, at = 0.52, font=2)

dev.off()


# Now, net div

pdf('figS10.pdf', width=10, height=5)

par(mfrow=c(1,2), mar=c(4,4,0,0))

# net div rate, correlation
plot.new()
plot.window(xlim = range(regimeSize), ylim = c(0,1))
axis(1, at = c(-1, axTicks(1)))
axis(2, at = c(-1, axTicks(2)))
mtext('minimum regime size', side = 1, line = 2.5, cex=1)
mtext('Pearson correlation', side = 2, line = 2.5, cex=1)
for (i in 1:length(rateMetrics)) {
	lines(regimeSize, regimeSizeTable[, paste0(rateMetricsNetDiv[i], 'netdiv_cor')], col= metricColors[i], lwd=2)
}
abline(h=1, lty=3)
mtext('(a)', outer = TRUE, line = -1.5, at = 0.02, font=2)

legend(35, 0.97, legend = c(expression(italic(r)['CRBD']), expression(lambda['TB']), expression(lambda['ND']), expression(lambda['DR']), expression(italic(r)['BAMM'])), col = metricColors, lwd=2, bty='n', xpd=NA, cex=1.2)


# net div rate, slope
plot.new()
plot.window(xlim = range(regimeSize), ylim = c(0,1.5))
axis(1, at = c(-1, axTicks(1)))
axis(2, at = c(-1, axTicks(2)))
mtext('minimum regime size', side = 1, line = 2.5, cex=1)
mtext('OLS slope', side = 2, line = 2.5, cex=1)

for (i in 1:length(rateMetrics)) {
	lines(regimeSize, regimeSizeTable[, paste0(rateMetricsNetDiv[i], 'netdiv_slope')], col= metricColors[i], lwd=2)
}
abline(h=1, lty=3)
mtext('(b)', outer = TRUE, line = -1.5, at = 0.52, font=2)

dev.off()









# Calculate absolute and proportional error as a function of regime size

allRatesList <- readRDS('dataFiles/allRatesList.rds')
listnames <- names(allRatesList)
listnames <- sapply(listnames, function(x) strsplit(x, '_')[[1]][1], USE.NAMES=FALSE)
allRatesList <- allRatesList[which(listnames %in% c('MeyerWiens2017','fossilBAMM','MooreEtAl2016','RaboskyEtAl2017'))]

# combine treeName and regimeID to create one big lookup table
for (i in 1:length(allRatesList)) {
	allRatesList[[i]]$regimeID <- paste0(allRatesList[[i]]$treeName, '_', allRatesList[[i]]$regimeID)
}
allRates <- do.call(rbind, allRatesList)
allRates$CRBDnetdiv <- allRates$CRBDlambda - allRates$CRBDmu

# average each regime 
regimeAvg <- matrix(nrow = length(unique(allRates$regimeID)), ncol = 11)
colnames(regimeAvg) <- c('treeName', 'regimeSize', 'tipLambda', 'tipNetDiv', 'lambdaTB', 'lambdaND', 'lambdaDR', 'lambdaBAMM', 'netDivBAMM', 'CRBDlambda', 'CRBDnetdiv')
regimeAvg <- as.data.frame(regimeAvg, stringsAsFactors=FALSE)
for (i in 1:length(unique(allRates$regimeID))) {
	
	qq <- allRates[which(allRates$regimeID == unique(allRates$regimeID)[i]),]
	regimeAvg[i, 'treeName'] <- qq[1, 'treeName']
	regimeAvg[i, 'regimeSize'] <- nrow(qq)
	regimeAvg[i, 'tipLambda'] <- mean(qq[, 'tipLambda'])
	regimeAvg[i, 'tipNetDiv'] <- mean(qq[, 'tipNetDiv'])
	regimeAvg[i, 'lambdaTB'] <- mean(qq[, 'lambdaTB'])
	regimeAvg[i, 'lambdaND'] <- mean(qq[, 'lambdaND'])
	regimeAvg[i, 'lambdaDR'] <- mean(qq[, 'lambdaDR'])
	regimeAvg[i, 'lambdaBAMM'] <- mean(qq[, 'lambdaBAMM'])
	regimeAvg[i, 'netDivBAMM'] <- mean(qq[, 'netDivBAMM'])
	regimeAvg[i, 'CRBDlambda'] <- mean(qq[, 'CRBDlambda'])
	regimeAvg[i, 'CRBDnetdiv'] <- mean(qq[, 'CRBDnetdiv'])
}

regimeAvgError <- matrix(nrow = nrow(regimeAvg), ncol = 12)
colnames(regimeAvgError) <- c('treeName','regimeSize','lambdaError_lambdaTB', 'lambdaError_lambdaND', 'lambdaError_lambdaDR', 'lambdaError_lambdaBAMM', 'lambdaError_CRBDlambda', 'netdivError_lambdaTB', 'netdivError_lambdaND', 'netdivError_lambdaDR', 'netdivError_netDivBAMM', 'netdivError_CRBDnetdiv')
regimeAvgError <- as.data.frame(regimeAvgError, stringsAsFactors=FALSE)

for (i in 1:nrow(regimeAvg)) {
	regimeAvgError[i, 'treeName'] <- regimeAvg[i, 'treeName']
	regimeAvgError[i, 'regimeSize'] <- regimeAvg[i, 'regimeSize']
	
	regimeAvgError[i, 'lambdaError_lambdaTB'] <- abs(regimeAvg[i, 'tipLambda'] - regimeAvg[i, 'lambdaTB'])
	regimeAvgError[i, 'lambdaError_lambdaND'] <- abs(regimeAvg[i, 'tipLambda'] - regimeAvg[i, 'lambdaND'])
	regimeAvgError[i, 'lambdaError_lambdaDR'] <- abs(regimeAvg[i, 'tipLambda'] - regimeAvg[i, 'lambdaDR'])
	regimeAvgError[i, 'lambdaError_lambdaBAMM'] <- abs(regimeAvg[i, 'tipLambda'] - regimeAvg[i, 'lambdaBAMM'])
	regimeAvgError[i, 'lambdaError_CRBDlambda'] <- abs(regimeAvg[i, 'tipLambda'] - regimeAvg[i, 'CRBDlambda'])

	regimeAvgError[i, 'netdivError_lambdaTB'] <- abs(regimeAvg[i, 'tipNetDiv'] - regimeAvg[i, 'lambdaTB'])
	regimeAvgError[i, 'netdivError_lambdaND'] <- abs(regimeAvg[i, 'tipNetDiv'] - regimeAvg[i, 'lambdaND'])
	regimeAvgError[i, 'netdivError_lambdaDR'] <- abs(regimeAvg[i, 'tipNetDiv'] - regimeAvg[i, 'lambdaDR'])
	regimeAvgError[i, 'netdivError_netDivBAMM'] <- abs(regimeAvg[i, 'tipNetDiv'] - regimeAvg[i, 'netDivBAMM'])
	regimeAvgError[i, 'netdivError_CRBDnetdiv'] <- abs(regimeAvg[i, 'tipNetDiv'] - regimeAvg[i, 'CRBDnetdiv'])
	
}

rateMetrics <- c('CRBDlambda','lambdaND','lambdaDR','lambdaBAMM')
rateMetricsNetDiv <- c('CRBDnetdiv','lambdaND','lambdaDR','netDivBAMM')

# bin into 10 tip groupings

brks <- c(1,10,20,30,40,50)

yrange <- c(0, 0.5)
width <- 0.4

rateMetricLabels <- c(expression(bold(lambda['CRBD'])), expression(bold(lambda['ND'])), expression(bold(lambda['DR'])),expression(bold(lambda['BAMM'])))

rateMetricLabelsNetDiv <- c(expression(bold(italic(r)['CRBD'])), expression(bold(lambda['TB'])), expression(bold(lambda['ND'])), expression(bold(lambda['DR'])),expression(bold(italic(r)['BAMM'])))

pdf('fig5.pdf', width=10, height=3)
par(mfrow = c(1, 4), mar = c(5,3,0,0), oma = c(0,1,1,0))


for (i in 1:length(rateMetrics)) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- intersect(grep('lambdaError', colnames(regimeAvgError)), grep(tipMetric, colnames(regimeAvgError)))
	plot.new()
	plot.window(xlim = c(0.5,6.5), ylim = yrange)
	axis(1, at = c(-0.5, 1:6), labels = FALSE)
	text(x = (1:6) + 0.35, y = -0.05, labels = c('1-10', '11-20','21-30','31-40','41-50','> 50'), srt=35, pos=2, xpd=NA)
	axis(2, at = c(-10, axTicks(2)))
	
	for (j in 1:length(brks)) {

		if (j < 6) {
			qStats <- quantile(regimeAvgError[which(regimeAvgError[, 'regimeSize'] > brks[j] & regimeAvgError[, 'regimeSize'] <= brks[j+1]), errorCol], c(0.05, 0.25, 0.5, 0.75, 0.95))
		} else {
			qStats <- quantile(regimeAvgError[which(regimeAvgError[, 'regimeSize'] > brks[j]), errorCol], c(0.05, 0.25, 0.5, 0.75, 0.95))
		}
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1, xpd =NA)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1, xpd =NA)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1, xpd =NA)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1, xpd =NA)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1, xpd =NA)
	}
	abline(h=0, lty=3)
	mtext('regime size', side = 1, cex=0.7, line = 3)
	if (i == 1) mtext('absolute error', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.6, cex.main=1.2)
}

dev.off()





