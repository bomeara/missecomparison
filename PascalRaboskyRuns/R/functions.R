#' Function to do a single MiSSE run
#'
#' @param dir Where the phylo object resides
#' @param nturnover How many turnover categories to use
#' @param neps_same Whether to have same number of eps categories (TRUE) or just one (FALSE)

DoSingleRun <- function(dir, nturnover=2, neps_same=TRUE, root_type="madfitz") {
#	print(system(paste0("ls data/pascal_rabosky_dryad/trees/", dir)))
#	print(paste0("data/pascal_rabosky_dryad/trees/", dir, "/", dir, ".tre"))
#	phy <- NULL
try(phy <- ape::read.tree(paste0("data/pascal_rabosky_dryad/trees/", dir, "/", dir, ".tre")))
	if(!is.null(phy)) {
	#print(phy)
	  #dir <- name(phy)
		turnover <- sequence(nturnover)
		eps <- ifelse(neps_same, turnover, rep(1, nturnover))
		hisse_result <- hisse::MiSSE(phy, f=1, turnover=turnover, eps=eps, root.type=root_type)
		hisse_recon <- hisse::MarginReconMiSSE(phy=phy, f=1,  hidden.states=nturnover, pars=hisse_result$solution, aic=hisse_result$AIC, root.type=root_type)
		summary_df <- t(hisse_recon$rates.mat %*% t(hisse_recon$tip.mat[,-1]))
		summary_df$treeName <- dir
		summary_df$taxon_id_in_phy <- sequence(nrow(summary_df))
		summary_df$tipName <- phy$tip.label
		summary_df$nturnover <- nturnover
		summary_df$neps <- length(unique(eps))
		summary_df$AIC <- hisse_result$AIC
		summary_df$AICc <- hisse_result$AICc
		summary_df$root_type <- hisse_result$root.type
		return(list(hisse_result=hisse_result, hisse_recon=hisse_recon, summary_df=summary_df))
	} else {
		return(list(failure=dir))
	}
}
