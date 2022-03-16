# plot tip rate values against each other, no by-tree or by-regime summarizing

require(LSD)

allRatesList <- readRDS('dataFiles/allRatesList.rds')
dat <- do.call(rbind, allRatesList)


table(dat$setname)


resolution <- 300
# ------------------------------------------------------------------
# Scatter figure for only ND, DR and BAMM, for speciation rate only. 

setlist <- list('multi-regime'= c('fossilBAMM', 'lambdaConstantVariance','MeyerWiens2017','MooreEtAl2016','RaboskyEtAl2017'), 'diversity-dependent'=c('Rabosky2014_DD_k1', 'Rabosky2014_DD_k2', 'Rabosky2014_DD_k3', 'Rabosky2014_DD_k4'), 'evolving rates'='evolvingRates')

tipMetrics <- c('lambdaND','lambdaDR','lambdaBAMM')

#gridVal <- 1000 # used in paper, but quite computationally intensive
gridVal <- 100

tipMetricLabel <- c(expression(lambda['ND']), expression(lambda['DR']),expression(lambda['BAMM']))

m <- rbind( c(10,10,10),
			c(1,2,3),
			c(11,11,11),
			c(4,5,6),
			c(12,12,12),
			c(7,8,9))

png('fig2.png', width=9, height=9, units='in', res=resolution)
			
layout(m, heights=rep(c(15/3, 85/3), times=3))
#layout.show(max(m))

par(mar=c(2,3,1,1), oma = c(2,2,0,0))

for (i in 1:length(setlist)) {
	qq <- dat[which(dat$setname %in% setlist[[i]]),]
	axisRange <- quantile(qq$tipLambda, c(0, 0.95))
	axisRange <- quantile(qq$tipLambda, c(0, 0.99))

	for (j in 1:length(tipMetrics)) {
		cat(i, '--', j, '\n')
		plot.new()
		plot.window(xlim = axisRange, ylim = axisRange)
		axis(1, at = c(-0.1, axTicks(1)))
		axis(2, at = c(-0.1, axTicks(2)))
		mtext(expression(paste(lambda['TRUE'])), side=1, line = 2.5)
		mtext(tipMetricLabel[j], side=2, line = 2.5)
		heatscatterpoints(qq$tipLambda, qq[, tipMetrics[j]], grid= gridVal, xlim = axisRange, ylim = axisRange)
		segments(axisRange[1], axisRange[1], axisRange[2], axisRange[2], lty=3)
		qqSub <- qq[which(!is.na(qq[, tipMetrics[j]])),]
		r <- cor(qqSub$tipLambda, qqSub[, tipMetrics[j]], method='spearman')
		mtext(paste('r = ', round(r,2)), side=1, line = -1.5, adj = 0.95, cex=0.7)
	}
}

titleSize <- 1.5
par(mar=c(1,1,1,1))
plot.new()
text(0.5, -1.5, 'multi-regime', font=2, xpd=NA, cex= titleSize)

plot.new()
text(0.5, -1.5, 'diversity-dependent', font=2, xpd=NA, cex= titleSize)

plot.new()
text(0.5, -1.5, 'evolving rates', font=2, xpd=NA, cex= titleSize)


dev.off()


#### Same for net diversification

tipMetricsNetDiv <- c('lambdaND','lambdaDR','netDivBAMM')

tipMetricLabel_netdiv <- c(expression(lambda['ND']), expression(lambda['DR']),expression(italic(r)['BAMM']))


png('figS7.png', width=9, height=9, units='in', res= resolution)
			
layout(m, heights=rep(c(15/3, 85/3), times=3))
#layout.show(max(m))

par(mar=c(2,3,1,1), oma = c(2,2,0,0))

for (i in 1:length(setlist)) {
	qq <- dat[which(dat$setname %in% setlist[[i]]),]
	axisRange <- quantile(qq$tipNetDiv, c(0, 0.99))
	axisRange <- list(c(-0.75, 1), c(-0.05, 0.15), c(0, 0.7))[[i]]

	for (j in 1:length(tipMetricsNetDiv)) {
		cat(i, '--', j, '\n')
		plot.new()
		plot.window(xlim = axisRange, ylim = axisRange)
		axis(1, at = c(-1, axTicks(1)), labels = c(NA, round(axTicks(1), 2)))
		axis(2, at = c(-1, axTicks(2)), labels = c(NA, round(axTicks(2), 2)))
		mtext(expression(italic(r)['TRUE']), side=1, line = 2.5)
		mtext(tipMetricLabel_netdiv[j], side=2, line = 2.5)
		heatscatterpoints(qq$tipNetDiv, qq[, tipMetricsNetDiv[j]], grid= gridVal, xlim = axisRange, ylim = axisRange)
		segments(axisRange[1], axisRange[1], axisRange[2], axisRange[2], lty=3)
		qqSub <- qq[which(!is.na(qq[, tipMetricsNetDiv[j]])),]
		r <- cor(qqSub$tipNetDiv, qqSub[, tipMetricsNetDiv[j]], method='spearman')
		mtext(paste('r = ', round(r,2)), side=1, line = -1.5, adj = 0.95, cex=0.7)
	}
}

titleSize <- 1.5
par(mar=c(1,1,1,1))
plot.new()
text(0.5, -1, 'multi-regime', font=2, xpd=NA, cex= titleSize)

plot.new()
text(0.5, -1, 'diversity-dependent', font=2, xpd=NA, cex= titleSize)

plot.new()
text(0.5, -1, 'evolving rates', font=2, xpd=NA, cex= titleSize)


dev.off()


##############

# For lambdaTB only

png('figS5.png', width=9, height=6, units='in', res= resolution)

par(mfrow = c(2,3), mar=c(2,3,1,1), oma = c(2,2,1,0))

for (i in 1:length(setlist)) {
	qq <- dat[which(dat$setname %in% setlist[[i]]),]
	axisRange <- quantile(qq$tipLambda, c(0, 0.99))

	cat(i, '\n')
	plot.new()
	plot.window(xlim = axisRange, ylim = axisRange)
	axis(1, at = c(-0.1, axTicks(1)))
	axis(2, at = c(-0.1, axTicks(2)))
	mtext(expression(lambda['TRUE']), side=1, line = 2.5)
	mtext(expression(lambda['TB']), side=2, line = 2.5)
	heatscatterpoints(qq$tipLambda, qq$lambdaTB, grid= gridVal, xlim = axisRange, ylim = axisRange)
	segments(axisRange[1], axisRange[1], axisRange[2], axisRange[2], lty=3)
	qqSub <- qq[which(!is.na(qq$lambdaTB)),]
	r <- cor(qqSub$tipLambda, qqSub$lambdaTB, method='spearman')
	mtext(paste('r = ', round(r,2)), side=1, line = -1.5, adj = 0.95, cex=0.7)
	title(main = names(setlist)[i], cex.main = 1.4)
}

for (i in 1:length(setlist)) {
	qq <- dat[which(dat$setname %in% setlist[[i]]),]
	axisRange <- list(c(-0.75, 1), c(-0.05, 0.15), c(0, 0.7))[[i]]

	cat(i, '\n')
	plot.new()
	plot.window(xlim = axisRange, ylim = axisRange)
	axis(1, at = c(-1, axTicks(1)), labels = c(NA, round(axTicks(1), 2)))
	axis(2, at = c(-1, axTicks(2)), labels = c(NA, round(axTicks(2), 2)))	
	mtext(expression(italic(r)['TRUE']), side=1, line = 2.5)
	mtext(expression(lambda['TB']), side=2, line = 2.5)
	heatscatterpoints(qq$tipNetDiv, qq$lambdaTB, grid= gridVal, xlim = axisRange, ylim = axisRange)
	segments(axisRange[1], axisRange[1], axisRange[2], axisRange[2], lty=3)
	qqSub <- qq[which(!is.na(qq$lambdaTB)),]
	r <- cor(qqSub$tipNetDiv, qqSub$lambdaTB, method='spearman')
	mtext(paste('r = ', round(r,2)), side=1, line = -1.5, adj = 0.95, cex=0.7)
}
dev.off()


# # 




# # By dataset
setlist <-c('lambdaConstantVariance', 'netDivConstantVariance', 'MooreEtAl2016', 'RaboskyEtAl2017', 'MeyerWiens2017', 'fossilBAMM')

sourceName <- c(expression(bold(constant~variance~lambda)), expression(bold(constant~variance~NetDiv)), expression(bold(Moore~et~al.~2016)), expression(bold(Rabosky~et~al.~2017)), expression(bold(Meyer~and~Wiens~2017)), expression(bold(Mitchell~et~al.2018)))
sourceName <- c(expression(lambda~uniform), expression(italic(r)~uniform), 'Moore et al. 2016', 'Rabosky et al. 2017', 'Meyer and Wiens 2017', 'Mitchell et al. 2018')

m <- rbind( c(1,2,7,8),
			c(3,4,9,10),
			c(5,6,11,12))


png('figS6.png', width=10, height=7, units='in', res= resolution)

layout(m)

axisRange <- c(0, 0.6)
tipMetrics <- c('lambdaDR','lambdaBAMM')
tipMetricLabel <- c(expression(lambda['DR']),expression(lambda['BAMM']))

for (i in 1:length(setlist)) {
	cat(i, '\n')
	qq <- dat[which(dat$setname %in% setlist[i]),]

	for (j in 1:length(tipMetrics)) {
		plot.new()
		if (j == 1) par(mar = c(4,5,4,1))
		if (j == 2) par(mar = c(4,3,4,3))
		plot.window(xlim = axisRange, ylim = axisRange)
		axis(1, c(-1, axTicks(1)))
		axis(2, c(-1, axTicks(2)), las=2)
		mtext(expression(lambda['TRUE']), side=1, line = 2.5)
		mtext(tipMetricLabel[j], side=2, line = 2.5)
		if (!all(is.na(qq[, tipMetrics[j]]))) {
			heatscatterpoints(qq$tipLambda, qq[, tipMetrics[j]], grid= gridVal, xlim = axisRange, ylim = axisRange)
			r <- cor(qq$tipLambda, qq[, tipMetrics[j]], method='spearman')
			mtext(paste('r = ', round(r,3)), side=1, line = -1.5, adj = 0.95, cex=0.8)
	
		}
		abline(a=0, b=1, lty=3)
	}
}


mtext(sourceName[1], line = -3, font=2, cex=0.8, at = 0.25, outer = TRUE)

mtext(sourceName[2], line = -21, font=2, cex=0.8, at = 0.25, outer = TRUE)

mtext(sourceName[3], line = -38, font=2, cex=0.8, at = 0.25, outer = TRUE)

mtext(sourceName[4], line = -3, font=2, cex=0.8, at = 0.75, outer = TRUE)

mtext(sourceName[5], line = -21, font=2, cex=0.8, at = 0.75, outer = TRUE)

mtext(sourceName[6], line = -38, font=2, cex=0.8, at = 0.75, outer = TRUE)

dev.off()
