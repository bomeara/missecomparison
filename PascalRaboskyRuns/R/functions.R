DoSingleRun <- function(phy, nturnover=2, neps=2) {
	max_length <- max(nturnover, neps)
	turnover <- sequence(nturnover)
	eps <- sequence(neps)
	if(length(turnover)<max_length) {
		turnover <- rep(1, max_length)
	}
	if(length(eps)<max_length) {
		eps <- rep(1, max_length)
	}
	hisse_result <- hisse::MiSSE(phy, f=1, turnover=turnover, eps=eps)
	hisse_recon <- MarginReconMiSSE(phy=phy, f=1,  hidden.states=max_length, pars=hisse_result$solution, aic=hisse_result$AIC)
	return(list(hisse_result=hisse_result, hisse_recon=hisse_recon))
}
