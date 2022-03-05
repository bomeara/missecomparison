# mode: julia 
	 using PANDA 
# mode: julia 
	 using JLD2 
# mode: julia 
	 my_tree = load_tree("/home/tvasconcelos/missecomparison/Empirical_example-Eucalypts/ClaDS_files/ML1_modified-Thornhill_et_al_2019.tre") 
# mode: julia 
	 output = infer_ClaDS(my_tree, print_state = 100) 
# mode: julia 
	 save_ClaDS_in_R(output, ""/home/tvasconcelos/missecomparison/Empirical_example-Eucalypts/ClaDS_files/Eucalypts_clads_results.Rdata") 