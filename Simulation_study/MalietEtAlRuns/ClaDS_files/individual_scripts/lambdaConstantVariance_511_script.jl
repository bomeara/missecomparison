# mode: julia 
	 using PANDA 
# mode: julia 
	 my_tree = load_tree("/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison/trees/lambdaConstantVariance_511.tre") 
# mode: julia 
	 output = infer_ClaDS(my_tree, print_state = 100) 
# mode: julia 
	 using JLD2 
# mode: julia 
	 save_ClaDS_in_R(output, "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison/results/lambdaConstantVariance_511_clads_results.Rdata") 