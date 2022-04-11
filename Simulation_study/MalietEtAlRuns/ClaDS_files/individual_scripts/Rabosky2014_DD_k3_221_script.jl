# mode: julia 
	 using PANDA 
# mode: julia 
	 my_tree = load_tree("/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison/trees/Rabosky2014_DD_k3_221.tre") 
# mode: julia 
	 output = infer_ClaDS(my_tree, print_state = 100) 
# mode: julia 
	 using JLD2 
# mode: julia 
	 save_ClaDS_in_R(output, "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison/results/Rabosky2014_DD_k3_221_clads_results.Rdata") 