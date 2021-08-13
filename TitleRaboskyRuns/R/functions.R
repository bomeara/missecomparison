#' Function to do a single MiSSE run
#'
#' @param dir Where the phylo object resides
#' @param phy The phylo object
#' @param root_type The type to use at the root
#' @param possibilities The possible_combos object
#' @param tree_index which tree is being used

DoSingleRun <- function(dir, phy, root_type="madfitz", possibilities, tree_index, n.cores=1, chunk.size=10) {
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
		#cat(paste0("Starting tree ", tree_index, " ntip is ", ape::Ntip(phy)), file=paste0("results/",unname(Sys.info()["nodename"]), "_",tree_index, "_newrun.log"), append=FALSE)

		output_files <- list.files("results", pattern=paste0("_",tree_index, ".rda"), full.names=TRUE)
		output_files <- output_files[!grepl("recon", output_files)]
		if(length(output_files)>=1) {
			load(output_files[1])
			hisse_result_all<- misse.list
		} else {
			hisse_result_all <- hisse::MiSSEGreedy(phy, f=1, root.type=root_type, possible.combos=possibilities, chunk.size=chunk.size, n.cores=1, save.file=paste0("results/",unname(Sys.info()["nodename"]), "_",tree_index, ".rda"), stop.deltaAICc=1000, sann=TRUE)
			#cat(paste0("Finished fit to tree ", tree_index), file=paste0("results/",unname(Sys.info()["nodename"]), "_",tree_index, "_newrun.log"), append=TRUE)
		}
		hisse_result_nonredundant <- PruneRedundantModels(hisse_result_all)
		AIC_weights <- hisse::GetAICWeights(hisse_result_nonredundant, criterion="AIC")
		delta_AIC <- sapply(hisse_result_nonredundant, "[[", "AIC") - min(sapply(hisse_result_nonredundant, "[[", "AIC"))
		AICc_weights <- hisse::GetAICWeights(hisse_result_nonredundant, criterion="AICc")
		delta_AICc <- sapply(hisse_result_nonredundant, "[[", "AICc") - min(sapply(hisse_result_nonredundant, "[[", "AICc"))
		save(list=ls(), file=paste0("results/",unname(Sys.info()["nodename"]), "_",tree_index, "_donefitting.rda"))
		model_fit_time <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
		for(model_index in sequence(length(hisse_result_nonredundant))) {
			if(delta_AICc[model_index]<20) {
				cat(paste0("Doing recon on model ", model_index, " at ", Sys.time(), "\n"), file=paste0("results/",unname(Sys.info()["nodename"]), "_",tree_index, ".log"), append=TRUE)

				start_time <- Sys.time()
				nturnover <- length(unique(hisse_result_nonredundant[[model_index]]$turnover))
				neps <- length(unique(hisse_result_nonredundant[[model_index]]$eps)) - ifelse(is.null(hisse_result_nonredundant[[model_index]]$fixed.eps), 0,1)


				hisse_recon <- hisse::MarginReconMiSSE(phy=hisse_result_nonredundant[[model_index]]$phy, f=1, hidden.states=nturnover, fixed.eps=hisse_result_nonredundant[[model_index]]$fixed.eps, pars=hisse_result_nonredundant[[model_index]]$solution, AIC=hisse_result_nonredundant[[model_index]]$AIC, root.type=root_type, get.tips.only=TRUE, n.cores=1)

				save(hisse_recon, hisse_result_nonredundant, hisse_result_all, file=paste0("results/", unname(Sys.info()["nodename"]), "_pre_summarizing_recon_",tree_index, "_model_", model_index, "_raw_.rda"))


				tip_mat_transformed <- hisse_recon$tip.mat[,-1]
				if(max(tip_mat_transformed) == 0) {
					tip_mat_transformed[,1] <- 1 #deal with misse bug of no weight if no hidden
				}
				summary_df_local <- t(hisse_recon$rates.mat %*% t(tip_mat_transformed))
				summary_df_local <- data.frame(summary_df_local)
				summary_df_local$treeName <- dir
				summary_df_local$tree_index <- tree_index
				summary_df_local$model_index <- model_index
				summary_df_local$length_all_models <- length(hisse_result_all)
				summary_df_local$length_nonredundant_models <- length(hisse_result_nonredundant)
				summary_df_local$length_nonredundant_delta_below_20 <- length(which(delta_AICc<20))

				summary_df_local$taxon_id_in_phy <- sequence(nrow(summary_df_local))
				summary_df_local$tipName <- phy$tip.label
				summary_df_local$nturnover <- nturnover
				summary_df_local$neps <- neps
				summary_df_local$nparam <- neps + nturnover
				summary_df_local$fixed_eps <- hisse_result_nonredundant[[model_index]]$fixed.eps
				summary_df_local$loglik <- hisse_result_nonredundant[[model_index]]$loglik

				summary_df_local$AIC <- hisse_result_nonredundant[[model_index]]$AIC
				summary_df_local$AICc <- hisse_result_nonredundant[[model_index]]$AICc
				summary_df_local$AIC_weight <- AIC_weights[model_index]
				summary_df_local$AICc_weight <- AICc_weights[model_index]
				summary_df_local$deltaAIC <- delta_AIC[model_index]
				summary_df_local$deltaAICc <- delta_AICc[model_index]
				summary_df_local$root_type <- hisse_result_nonredundant[[model_index]]$root.type
				summary_df_local$elapsed_mins_recon <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
				summary_df_local$elapsed_mins_params_all_models_together <- model_fit_time
				summary_df_local$ntip <- ape::Ntip(hisse_result_nonredundant[[model_index]]$phy)
				summary_df_local$treeHeight <- max(branching.times(hisse_result_nonredundant[[model_index]]$phy))
				summary_df_local$min_brlen <- min(hisse_result_nonredundant[[model_index]]$phy$edge.length)
				#return(list(hisse_result=hisse_result, hisse_recon=hisse_recon, summary_df=summary_df))
				#save(summary_df_local, hisse_recon, file=paste0("results/", unname(Sys.info()["nodename"]), "_model_", model_index, "_recon_",tree_index, "_newrun.rda"))
				if(nrow(summary_df)==0) {
					try(summary_df <- summary_df_local)
				} else {
					try(summary_df <- plyr::rbind.fill(summary_df, summary_df_local))
				}
				save(summary_df, hisse_result_nonredundant, AICc_weights, delta_AICc, file=paste0("results/", unname(Sys.info()["nodename"]), "_post_summarizing_recon_",tree_index,"_model_", model_index, "_newrun.rda"))
			}
		}
		save(summary_df, hisse_result_nonredundant, hisse_result_all, AICc_weights, delta_AICc, model_fit_time, file=paste0("results/", unname(Sys.info()["nodename"]), "_done_recon_",tree_index, "_newrun.rda"))
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
