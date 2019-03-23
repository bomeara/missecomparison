#' Function to do a single MiSSE run
#'
#' @param dir Where the phylo object resides
#' @param phy The phylo object
#' @param nturnover How many turnover categories to use
#' @param neps_same Whether to have same number of eps categories (TRUE) or just one (FALSE)

DoSingleRun <- function(dir, phy, nturnover=2, neps_same=TRUE, root_type="madfitz") {
#	print(system(paste0("ls data/pascal_rabosky_dryad/trees/", dir)))
#	print(paste0("data/pascal_rabosky_dryad/trees/", dir, "/", dir, ".tre"))
#	phy <- NULL
#try(phy <- ape::read.tree(paste0("data/pascal_rabosky_dryad/trees/", dir, "/", dir, ".tre")))
	start_time <- Sys.time()
	if(!is.null(phy)) {
	#print(phy)
	  #dir <- name(phy)
		turnover <- sequence(nturnover)
		eps <- turnover
		if(neps_same) {
			eps <- rep(1, nturnover)
		}
		#eps <- ifelse(neps_same, turnover, rep(1, nturnover))
		hisse_result <- hisse::MiSSE(phy, f=1, turnover=turnover, eps=eps, root.type=root_type)
		hisse_recon <- hisse::MarginReconMiSSE(phy=phy, f=1,  hidden.states=nturnover, pars=hisse_result$solution, aic=hisse_result$AIC, root.type=root_type, n.cores=1)
		tip_mat_transformed <- hisse_recon$tip.mat[,-1]
		if(max(tip_mat_transformed) == 0) {
			tip_mat_transformed[,1] <- 1 #deal with misse bug of no weight if no hidden
		}
		summary_df <- t(hisse_recon$rates.mat %*% t(tip_mat_transformed))
		summary_df <- data.frame(summary_df)
		summary_df$treeName <- dir
		summary_df$taxon_id_in_phy <- sequence(nrow(summary_df))
		summary_df$tipName <- phy$tip.label
		summary_df$nturnover <- nturnover
		summary_df$neps <- length(unique(eps))
		summary_df$AIC <- hisse_result$AIC
		summary_df$AICc <- hisse_result$AICc
		summary_df$root_type <- hisse_result$root.type
		summary_df$elapsed_mins <- as.numeric(difftime(Sys.time, start_time, units="mins"))
		#return(list(hisse_result=hisse_result, hisse_recon=hisse_recon, summary_df=summary_df))
		return(summary_df)
	} else {
		return(list(failure=dir))
	}
}
