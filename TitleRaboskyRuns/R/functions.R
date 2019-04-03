#' Function to do a single MiSSE run
#'
#' @param dir Where the phylo object resides
#' @param phy The phylo object
#' @param root_type The type to use at the root

DoSingleRun <- function(dir, phy, root_type="madfitz") {
#	print(system(paste0("ls data/title_rabosky_dryad/trees/", dir)))
#	print(paste0("data/title_rabosky_dryad/trees/", dir, "/", dir, ".tre"))
#	phy <- NULL
#try(phy <- ape::read.tree(paste0("data/title_rabosky_dryad/trees/", dir, "/", dir, ".tre")))
	start_time <- Sys.time()
	if(!is.null(phy)) {
		summary_df <- data.frame()
	#print(phy)
	  #dir <- name(phy)
		#eps <- ifelse(neps_same, turnover, rep(1, nturnover))
		hisse_result_all <- hisse::MiSSEGreedy(phy, f=1, root.type=root_type)
		model_fit_time <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
		for(model_index in sequence(length(hisse_result_all))) {
			start_time <- Sys.time()
			nturnover <- length(unique(hisse_result_all[[model_index]]$turnover))
			neps <- length(unique(hisse_result_all[[model_index]]$eps))

			hisse_recon <- hisse::MarginReconMiSSE(phy=hisse_result_all[[model_index]]$phy, f=1,  hidden.states=nturnover, pars=hisse_result_all[[model_index]]$solution, aic=hisse_result_all[[model_index]]$AIC, root.type=root_type, get.tips.only=TRUE, n.cores=parallel::detectCores())
			tip_mat_transformed <- hisse_recon$tip.mat[,-1]
			if(max(tip_mat_transformed) == 0) {
				tip_mat_transformed[,1] <- 1 #deal with misse bug of no weight if no hidden
			}
			summary_df_local <- t(hisse_recon$rates.mat %*% t(tip_mat_transformed))
			summary_df_local <- data.frame(summary_df_local)
			summary_df_local$treeName <- dir
			summary_df_local$taxon_id_in_phy <- sequence(nrow(summary_df_local))
			summary_df_local$tipName <- phy$tip.label
			summary_df_local$nturnover <- nturnover
			summary_df_local$neps <- neps
			summary_df_local$AIC <- hisse_result_all[[model_index]]$AIC
			summary_df_local$AICc <- hisse_result_all[[model_index]]$AICc
			summary_df_local$root_type <- hisse_result_all[[model_index]]$root.type
			summary_df_local$elapsed_mins_recon <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
			summary_df_local$elapsed_mins_params_all_models_together <- model_fit_time
			#return(list(hisse_result=hisse_result, hisse_recon=hisse_recon, summary_df=summary_df))
			if(nrow(summary_df)==0) {
				summary_df <- summary_df_local
			} else {
				summary_df <- rbind(summary_df, summary_df_local)
			}
		}
		return(summary_df)
	} else {
		return(list(failure=dir))
	}
}
