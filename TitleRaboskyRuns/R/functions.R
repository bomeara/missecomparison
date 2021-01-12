#' Function to do a single MiSSE run
#'
#' @param dir Where the phylo object resides
#' @param phy The phylo object
#' @param root_type The type to use at the root
#' @param possibilities The possible_combos object
#' @param tree_index which tree is being used

DoSingleRun <- function(dir, phy, root_type="madfitz", possibilities, tree_index) {
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
		hisse_result_all <- hisse::MiSSEGreedy(phy, f=1, root.type=root_type, possible.combos=possibilities, chunk.size=6, n.cores=parallel::detectCores(), save.file=paste0(unname(Sys.info()["nodename"]), "_",tree_index, ".rda"))
		AIC_weights <- hisse::GetAICWeights(hisse_result_all, criterion="AIC")
		delta_AIC <- sapply(hisse_result_all, "[[", "AIC") - min(sapply(hisse_result_all, "[[", "AIC"))

		AICc_weights <- hisse::GetAICWeights(hisse_result_all, "AICc")
		delta_AICc <- sapply(hisse_result_all, "[[", "AICc") - min(sapply(hisse_result_all, "[[", "AICc"))

		model_fit_time <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
		for(model_index in sequence(length(hisse_result_all))) {
			start_time <- Sys.time()
			nturnover <- length(unique(hisse_result_all[[model_index]]$turnover))
			neps <- length(unique(hisse_result_all[[model_index]]$eps))

			hisse_recon <- hisse::MarginReconMiSSE(phy=hisse_result_all[[model_index]]$phy, f=1, hidden.states=nturnover, pars=hisse_result_all[[model_index]]$solution, AIC=hisse_result_all[[model_index]]$AIC, root.type=root_type, get.tips.only=TRUE)
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
			summary_df_local$AIC_weight <- AIC_weights[model_index]
			summary_df_local$AICc_weight <- AICc_weights[model_index]
			summary_df_local$deltaAIC <- delta_AIC[model_index]
			summary_df_local$deltaAICc <- delta_AICc[model_index]
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

WrapSingleRunCondor <- function(dir, phy, root_type="madfitz") {
	file_root <- paste0(gsub("/","",dir), "_", gsub("/","",phy))
	executable <- paste0("#! /bin/sh", "\n", "R CMD BATCH batch_", file_root, ".R")
	cat(executable, file=paste0("exec_",file_root,".sh"))
	system(paste0("chmod u+x ", paste0("exec_",file_root,".sh")))

	submit_script <- paste0("executable=", paste0("exec_",file_root,".sh"),"\n",
		"universe=vanilla\n",
		"output=results.output.", file_root, "\n",
		"error=results.error.", file_root, "\n",
		"transfer_input_files=batch_", , file_root, ".R", "\n",
		"notification=never\nshould_transfer_files=YES\nwhen_to_transfer_output=ON_EXIT\n"
	)
	cat(submit_script, file=paste0("job_", file_root))

	R_script <- paste0("library(hisse)\n",
		# call in DoSingleRun, have it return stuff
	)


}
