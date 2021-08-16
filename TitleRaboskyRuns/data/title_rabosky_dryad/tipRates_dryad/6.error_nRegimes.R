# Analyze error in rates as a function of number of regimes

# only for speciation rate, as at this point, these metrics have been shown to track speciation rate. 

dat <- read.csv('dataFiles/rateErrorMetrics.csv', stringsAsFactors=FALSE)

# set infinite values to NA
for (i in 6:ncol(dat)) {
	if (length(which(!is.finite(dat[,i]))) > 0) {
		dat[which(!is.finite(dat[,i])), i] <- NA
	}
}
# ---------------------------------------------------------
# ---------------------------------------------------------

# combined figure of all 5 tip rate metrics, under 3 scenarios
## leaving lambdaTB out for a separate figure

rateMetrics <- c('CRBDlambda', 'lambdaND', 'lambdaDR', 'lambdaBAMM')
errorMetricNames <- c(PE1='Prop. Error 1', PE2='Prop. Error 2', PE3='Prop. Error 3', PE4='Prop. 4', r2='OLS R2', slope='OLS slope', absoluteError='mean absolute error')

rateMetricLabels <- c(expression(bold(paste(lambda['CRBD']))), expression(bold(paste(lambda['ND']))), expression(bold(paste(lambda['DR']))),expression(bold(paste(lambda['BAMM']))))

pdf('fig3.pdf', width=10, height=7)
par(mfrow = c(3, 4), mar = c(4,3,3,1), oma = c(0,1,0,0))

width <- 0.4
rateType <- 'Lambda'
qrange <- c(0.05, 0.95)
defaultYrange <- c(0, 0.4)
xAxisCex <- 1
yAxisCex <- 1
errorMetric <- 'absoluteError'


## time-constant multi-regime
dat2 <- dat[which(dat$setname %in% c('MooreEtAl2016','RaboskyEtAl2017','fossilBAMM','MeyerWiens2017')),]

datSplit <- split(dat2, dat2$nRegimes)

# collapse shift number >= 10 into same category
datSplit2 <- datSplit
datSplit <- datSplit2[as.character(1:10)]
datSplit[['10']] <- do.call(rbind, datSplit2[which(names(datSplit2) == '10'):length(datSplit2)])
names(datSplit)[10] <- '>10'

for (i in 1:4) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- Reduce(intersect, list(grep(rateType, colnames(dat)), grep(rateMetrics[i], colnames(dat)), grep(errorMetric, colnames(dat))))
	if (tipMetric != 'lambdaTB') {
		yrange <- defaultYrange
	} else {
		yrange <- c(0, quantile(unlist(sapply(datSplit, function(x) x[, errorCol])), 0.975))
	}
	plot.new()
	plot.window(xlim = c(1, length(datSplit)), ylim = yrange)
	axis(1, at = c(0, 1:length(datSplit)), labels = NA)
	axis(1, at = c(0, 1:length(datSplit)), tick=FALSE, labels = c(NA, names(datSplit)), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(1, at = c(0, 1:length(datSplit)), tick=FALSE, labels = c(rep(NA,10), names(datSplit)[10]), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(2, at = c(-10, axTicks(2)), cex.axis = yAxisCex, mgp = c(3, 0.7, 0))
	for (j in 1:length(datSplit)) {
		qStats <- quantile(datSplit[[j]][, errorCol], c(qrange[1], 0.25, 0.5, 0.75, qrange[2]), na.rm=TRUE)
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1)
	}
	abline(h=0, lty=3)
	mtext('n regimes', side = 1, cex=0.7, line = 2)
	if (i == 1) mtext('mean absolute error', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.5)
}
text(x=grconvertX(0.5, "ndc", "user"), y=grconvertY(0.97, "ndc", "user"), 'multi-regime', font=2, cex=1.5, xpd=NA)


# multi-regime diversity-dependent
dat2 <- dat[which(dat$setname %in% c('Rabosky2014_DD_k1', 'Rabosky2014_DD_k2', 'Rabosky2014_DD_k3', 'Rabosky2014_DD_k4')),]

datSplit <- split(dat2, dat2$nRegimes)

defaultYrange <- c(0, 0.2)

for (i in 1:4) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- Reduce(intersect, list(grep(rateType, colnames(dat)), grep(rateMetrics[i], colnames(dat)), grep(errorMetric, colnames(dat))))
	if (tipMetric != 'lambdaTB') {
		yrange <- defaultYrange
	} else {
		yrange <- c(0, quantile(unlist(sapply(datSplit, function(x) x[, errorCol])), 0.975))
	}
	plot.new()
	plot.window(xlim = c(1.5, 5.5), ylim = yrange)
	axis(1, at = c(1, as.numeric(names(datSplit))), labels = NA)
	axis(1, at = c(1, as.numeric(names(datSplit))), tick=FALSE, labels = c(NA, names(datSplit)), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(2, at = c(-0.25, axTicks(2)), cex.axis = yAxisCex, mgp = c(3, 0.5, 0))
	for (j in as.numeric(names(datSplit))) {
		qStats <- quantile(datSplit[[as.character(j)]][, errorCol], c(qrange[1], 0.25, 0.5, 0.75, qrange[2]), na.rm=TRUE)
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1)
	}
	abline(h=0, lty=3)
	mtext('n regimes', side = 1, cex=0.7, line = 2)
	if (i == 1) mtext('mean absolute error', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.5)
}
text(x=grconvertX(0.5, "ndc", "user"), y=grconvertY(0.64, "ndc", "user"), 'diversity-dependent', font=2, cex=1.5, xpd=NA)


# Evolving Rates

treeDat <- read.csv('dataFiles/treeSummary.csv', stringsAsFactors=FALSE)

dat2 <- dat[which(dat$setname == 'evolvingRates'),]
dat2$sigma <- sapply(dat2$treeName, function(x) treeDat[which(treeDat$treeName == x), 'sigma'])
table(dat2$sigma)

datSplit <- split(dat2, dat2$sigma)

defaultYrange <- c(0, 0.4)

for (i in 1:4) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- Reduce(intersect, list(grep(rateType, colnames(dat)), grep(rateMetrics[i], colnames(dat)), grep(errorMetric, colnames(dat))))
	if (tipMetric != 'lambdaTB') {
		yrange <- defaultYrange
	} else {
		yrange <- c(0, quantile(unlist(sapply(datSplit, function(x) x[, errorCol])), 0.975))
	}
	plot.new()
	plot.window(xlim = c(0.5,3.5), ylim = yrange)
	axis(1, at = c(-0.5, 1:3), labels = NA)
	axis(1, at = c(-0.5, 1:3), tick=FALSE, labels = c(NA, names(datSplit)), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(2, at = c(-10, axTicks(2)), cex.axis = yAxisCex, mgp = c(3, 0.5, 0))
	for (j in 1:length(datSplit)) {
		qStats <- quantile(datSplit[[j]][, errorCol], c(qrange[1], 0.25, 0.5, 0.75, qrange[2]), na.rm=TRUE)
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1)
	}
	abline(h=0, lty=3)
	mtext(expression(sigma), side = 1, cex=0.7, line = 2)
	if (i == 1) mtext('mean absolute error', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.5)
}
text(x=grconvertX(0.5, "ndc", "user"), y=grconvertY(0.3, "ndc", "user"), 'evolving rates', font=2, cex=1.5, xpd=NA)
dev.off()



#########################
# Supplemental figure using RMSE (root mean square error)

pdf('figS8.pdf', width=10, height=7)
par(mfrow = c(3, 4), mar = c(4,3,3,1), oma = c(0,1,0,0))

width <- 0.4
rateType <- 'Lambda'
qrange <- c(0.05, 0.95)
defaultYrange <- c(0, 0.8)
errorMetric <- 'rmse'


## time-constant multi-regime
dat2 <- dat[which(dat$setname %in% c('MooreEtAl2016','RaboskyEtAl2017','fossilBAMM','MeyerWiens2017')),]

datSplit <- split(dat2, dat2$nRegimes)

# collapse shift number >= 10 into same category
datSplit2 <- datSplit
datSplit <- datSplit2[as.character(1:10)]
datSplit[['10']] <- do.call(rbind, datSplit2[which(names(datSplit2) == '10'):length(datSplit2)])
names(datSplit)[10] <- '>10'

for (i in 1:4) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- Reduce(intersect, list(grep(rateType, colnames(dat)), grep(rateMetrics[i], colnames(dat)), grep(errorMetric, colnames(dat))))
	if (tipMetric != 'lambdaTB') {
		yrange <- defaultYrange
	} else {
		yrange <- c(0, quantile(unlist(sapply(datSplit, function(x) x[, errorCol])), 0.975))
	}
	plot.new()
	plot.window(xlim = c(1, length(datSplit)), ylim = yrange)
	axis(1, at = c(0, 1:length(datSplit)), labels = NA)
	axis(1, at = c(0, 1:length(datSplit)), tick=FALSE, labels = c(NA, names(datSplit)), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(1, at = c(0, 1:length(datSplit)), tick=FALSE, labels = c(rep(NA,10), names(datSplit)[10]), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(2, at = c(-10, axTicks(2)), cex.axis = yAxisCex, mgp = c(3, 0.7, 0))
	for (j in 1:length(datSplit)) {
		qStats <- quantile(datSplit[[j]][, errorCol], c(qrange[1], 0.25, 0.5, 0.75, qrange[2]), na.rm=TRUE)
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1)
	}
	abline(h=0, lty=3)
	mtext('n regimes', side = 1, cex=0.7, line = 2)
	if (i == 1) mtext('RMSE', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.5)
}
text(x=grconvertX(0.5, "ndc", "user"), y=grconvertY(0.97, "ndc", "user"), 'multi-regime', font=2, cex=1.5, xpd=NA)


# multi-regime diversity-dependent
dat2 <- dat[which(dat$setname %in% c('Rabosky2014_DD_k1', 'Rabosky2014_DD_k2', 'Rabosky2014_DD_k3', 'Rabosky2014_DD_k4')),]

datSplit <- split(dat2, dat2$nRegimes)

defaultYrange <- c(0, 0.2)

for (i in 1:4) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- Reduce(intersect, list(grep(rateType, colnames(dat)), grep(rateMetrics[i], colnames(dat)), grep(errorMetric, colnames(dat))))
	if (tipMetric != 'lambdaTB') {
		yrange <- defaultYrange
	} else {
		yrange <- c(0, quantile(unlist(sapply(datSplit, function(x) x[, errorCol])), 0.975))
	}
	plot.new()
	plot.window(xlim = c(1.5, 5.5), ylim = yrange)
	axis(1, at = c(1, as.numeric(names(datSplit))), labels = NA)
	axis(1, at = c(1, as.numeric(names(datSplit))), tick=FALSE, labels = c(NA, names(datSplit)), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(2, at = c(-0.10, axTicks(2)), cex.axis = yAxisCex, mgp = c(3, 0.5, 0))
	for (j in as.numeric(names(datSplit))) {
		qStats <- quantile(datSplit[[as.character(j)]][, errorCol], c(qrange[1], 0.25, 0.5, 0.75, qrange[2]), na.rm=TRUE)
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1)
	}
	abline(h=0, lty=3)
	mtext('n regimes', side = 1, cex=0.7, line = 2)
	if (i == 1) mtext('RMSE', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.5)
}
text(x=grconvertX(0.5, "ndc", "user"), y=grconvertY(0.64, "ndc", "user"), 'diversity-dependent', font=2, cex=1.5, xpd=NA)


# Evolving Rates

treeDat <- read.csv('dataFiles/treeSummary.csv', stringsAsFactors=FALSE)

dat2 <- dat[which(dat$setname == 'evolvingRates'),]
dat2$sigma <- sapply(dat2$treeName, function(x) treeDat[which(treeDat$treeName == x), 'sigma'])
table(dat2$sigma)

datSplit <- split(dat2, dat2$sigma)

defaultYrange <- c(0, 0.8)

for (i in 1:4) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- Reduce(intersect, list(grep(rateType, colnames(dat)), grep(rateMetrics[i], colnames(dat)), grep(errorMetric, colnames(dat))))
	if (tipMetric != 'lambdaTB') {
		yrange <- defaultYrange
	} else {
		yrange <- c(0, quantile(unlist(sapply(datSplit, function(x) x[, errorCol])), 0.975))
	}
	plot.new()
	plot.window(xlim = c(0.5,3.5), ylim = yrange)
	axis(1, at = c(-0.5, 1:3), labels = NA)
	axis(1, at = c(-0.5, 1:3), tick=FALSE, labels = c(NA, names(datSplit)), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(2, at = c(-10, axTicks(2)), cex.axis = yAxisCex, mgp = c(3, 0.5, 0))
	for (j in 1:length(datSplit)) {
		qStats <- quantile(datSplit[[j]][, errorCol], c(qrange[1], 0.25, 0.5, 0.75, qrange[2]), na.rm=TRUE)
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1)
	}
	abline(h=0, lty=3)
	mtext(expression(sigma), side = 1, cex=0.7, line = 2)
	if (i == 1) mtext('RMSE', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.5)
}
text(x=grconvertX(0.5, "ndc", "user"), y=grconvertY(0.3, "ndc", "user"), 'evolving rates', font=2, cex=1.5, xpd=NA)
dev.off()



#########################
# Supplemental version with lambdaTB, compared to lambdaDR and lambdaBAMM


rateMetrics <- c('lambdaTB', 'lambdaDR', 'lambdaBAMM')
errorMetricNames <- c(PE1='Prop. Error 1', PE2='Prop. Error 2', PE3='Prop. Error 3', PE4='Prop. 4', r2='OLS R2', slope='OLS slope', absoluteError='mean absolute error')

rateMetricLabels <- c(expression(bold(paste(lambda['TB']))), expression(bold(paste(lambda['DR']))),expression(bold(paste(lambda['BAMM']))))

pdf('figS9.pdf', width=7.5, height=7)
par(mfrow = c(3, 3), mar = c(4,3,3,1), oma = c(0,1,0,0))

width <- 0.4
rateType <- 'Lambda'
qrange <- c(0.05, 0.95)
defaultYrange <- c(0, 8)
errorMetric <- 'absoluteError'


## time-constant multi-regime
dat2 <- dat[which(dat$setname %in% c('MooreEtAl2016','RaboskyEtAl2017','fossilBAMM','MeyerWiens2017')),]

datSplit <- split(dat2, dat2$nRegimes)

# collapse shift number >= 10 into same category
datSplit2 <- datSplit
datSplit <- datSplit2[as.character(1:10)]
datSplit[['10']] <- do.call(rbind, datSplit2[which(names(datSplit2) == '10'):length(datSplit2)])
names(datSplit)[10] <- '>10'

for (i in 1:3) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- Reduce(intersect, list(grep(rateType, colnames(dat)), grep(rateMetrics[i], colnames(dat)), grep(errorMetric, colnames(dat))))
	plot.new()
	plot.window(xlim = c(1, length(datSplit)), ylim = defaultYrange)
	axis(1, at = c(0, 1:length(datSplit)), labels = NA)
	axis(1, at = c(0, 1:length(datSplit)), tick=FALSE, labels = c(NA, names(datSplit)), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(1, at = c(0, 1:length(datSplit)), tick=FALSE, labels = c(rep(NA,10), names(datSplit)[10]), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(2, at = c(-10, axTicks(2)))
	for (j in 1:length(datSplit)) {
		qStats <- quantile(datSplit[[j]][, errorCol], c(qrange[1], 0.25, 0.5, 0.75, qrange[2]), na.rm=TRUE)
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1)
	}
	abline(h=0, lty=3)
	mtext('n regimes', side = 1, cex=0.7, line = 2)
	if (i == 1) mtext('mean absolute error', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.5)
}
text(x=grconvertX(0.5, "ndc", "user"), y=grconvertY(0.97, "ndc", "user"), 'multi-regime', font=2, cex=1.5, xpd=NA)


# multi-regime diversity-dependent
dat2 <- dat[which(dat$setname %in% c('Rabosky2014_DD_k1', 'Rabosky2014_DD_k2', 'Rabosky2014_DD_k3', 'Rabosky2014_DD_k4')),]

datSplit <- split(dat2, dat2$nRegimes)

defaultYrange <- c(0, 1.5)

for (i in 1:3) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- Reduce(intersect, list(grep(rateType, colnames(dat)), grep(rateMetrics[i], colnames(dat)), grep(errorMetric, colnames(dat))))
	if (tipMetric != 'lambdaTB') {
		yrange <- defaultYrange
	} else {
		yrange <- c(0, quantile(unlist(sapply(datSplit, function(x) x[, errorCol])), 0.975))
	}
	plot.new()
	plot.window(xlim = c(1.5, 5.5), ylim = defaultYrange)
	axis(1, at = c(1, as.numeric(names(datSplit))), labels = NA)
	axis(1, at = c(1, as.numeric(names(datSplit))), tick=FALSE, labels = c(NA, names(datSplit)), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(2, at = c(-1, axTicks(2)))
	for (j in as.numeric(names(datSplit))) {
		qStats <- quantile(datSplit[[as.character(j)]][, errorCol], c(qrange[1], 0.25, 0.5, 0.75, qrange[2]), na.rm=TRUE)
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1)
	}
	abline(h=0, lty=3)
	mtext('n regimes', side = 1, cex=0.7, line = 2)
	if (i == 1) mtext('mean absolute error', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.5)
}
text(x=grconvertX(0.5, "ndc", "user"), y=grconvertY(0.64, "ndc", "user"), 'diversity-dependent', font=2, cex=1.5, xpd=NA)


# Evolving Rates

treeDat <- read.csv('dataFiles/treeSummary.csv', stringsAsFactors=FALSE)

dat2 <- dat[which(dat$setname == 'evolvingRates'),]
dat2$sigma <- sapply(dat2$treeName, function(x) treeDat[which(treeDat$treeName == x), 'sigma'])
table(dat2$sigma)

datSplit <- split(dat2, dat2$sigma)

defaultYrange <- c(0, 8)

for (i in 1:3) {
	
	tipMetric <- rateMetrics[i]
	errorCol <- Reduce(intersect, list(grep(rateType, colnames(dat)), grep(rateMetrics[i], colnames(dat)), grep(errorMetric, colnames(dat))))
	if (tipMetric != 'lambdaTB') {
		yrange <- defaultYrange
	} else {
		yrange <- c(0, quantile(unlist(sapply(datSplit, function(x) x[, errorCol])), 0.975))
	}
	plot.new()
	plot.window(xlim = c(0.5,3.5), ylim = defaultYrange)
	axis(1, at = c(-0.5, 1:3), labels = NA)
	axis(1, at = c(-0.5, 1:3), tick=FALSE, labels = c(NA, names(datSplit)), lwd=0, cex.axis = xAxisCex, mgp = c(3, 0.7, 0), xpd=NA)
	axis(2, at = c(-10, axTicks(2)))
	for (j in 1:length(datSplit)) {
		qStats <- quantile(datSplit[[j]][, errorCol], c(qrange[1], 0.25, 0.5, 0.75, qrange[2]), na.rm=TRUE)
		rect(j - width/2, qStats[2], j + width/2, qStats[4], col=gray(0.75))
		segments(j, qStats[1], j, qStats[2], lty=2, lend=1, xpd=NA)
		segments(j, qStats[4], j, qStats[5], lty=2, lend=1, xpd=NA)
		segments(j - width/3, qStats[1], j + width/3, qStats[1], lend=1, xpd=NA)
		segments(j - width/3, qStats[5], j + width/3, qStats[5], lend=1, xpd=NA)
		segments(j - width/3, qStats[3], j + width/3, qStats[3], lwd=2, lend=1, xpd=NA)
	}
	abline(h=0, lty=3)
	mtext(expression(sigma), side = 1, cex=0.7, line = 2)
	if (i == 1) mtext('mean absolute error', side = 2, cex=0.7, line = 2.5)
	title(main = rateMetricLabels[i], line = -0.5)
}
text(x=grconvertX(0.5, "ndc", "user"), y=grconvertY(0.3, "ndc", "user"), 'evolving rates', font=2, cex=1.5, xpd=NA)
dev.off()



# Explore accuracy in DR and BAMM as compared to CRBD, when all trees are combined:


## time-constant multi-regime
dat2 <- dat[which(dat$setname %in% c('MooreEtAl2016','RaboskyEtAl2017','fossilBAMM','MeyerWiens2017')),]

crbd <- dat2[, 'tipLambda_CRBDlambda_absoluteError']
dr <- dat2[, 'tipLambda_lambdaDR_absoluteError']
bamm <- dat2[, 'tipLambda_lambdaBAMM_absoluteError']

cat('time-constant multi-regime\n')

cat('\tCRBD more accurate than DR', round(length(which((crbd < dr) == TRUE)) / length(crbd), 2), '\n')
cat('\tCRBD more accurate than BAMM', round(length(which((crbd < bamm) == TRUE)) / length(crbd), 2), '\n\n')


# multi-regime diversity-dependent
dat2 <- dat[which(dat$setname %in% c('Rabosky2014_DD_k1', 'Rabosky2014_DD_k2', 'Rabosky2014_DD_k3', 'Rabosky2014_DD_k4')),]

crbd <- dat2[, 'tipLambda_CRBDlambda_absoluteError']
dr <- dat2[, 'tipLambda_lambdaDR_absoluteError']
bamm <- dat2[, 'tipLambda_lambdaBAMM_absoluteError']


cat('multi-regime diversity-dependent\n')

cat('\tCRBD more accurate than DR', round(length(which((crbd < dr) == TRUE)) / length(crbd), 2), '\n')
cat('\tCRBD more accurate than BAMM', round(length(which((crbd < bamm) == TRUE)) / length(crbd), 2), '\n\n')



treeDat <- read.csv('dataFiles/treeSummary.csv', stringsAsFactors=FALSE)

dat2 <- dat[which(dat$setname == 'evolvingRates'),]
dat2$sigma <- sapply(dat2$treeName, function(x) treeDat[which(treeDat$treeName == x), 'sigma'])
table(dat2$sigma)

crbd <- dat2[, 'tipLambda_CRBDlambda_absoluteError']
dr <- dat2[, 'tipLambda_lambdaDR_absoluteError']
bamm <- dat2[, 'tipLambda_lambdaBAMM_absoluteError']

cat('evolving rates\n')

cat('\tCRBD more accurate than DR', round(length(which((crbd < dr) == TRUE)) / length(crbd), 2), '\n')
cat('\tCRBD more accurate than BAMM', round(length(which((crbd < bamm) == TRUE)) / length(crbd), 2), '\n\n')



