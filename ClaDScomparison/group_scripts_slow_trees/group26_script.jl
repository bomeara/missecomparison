# mode: julia 
	 using PANDA 
# mode: julia 
	 using JLD2 
# mode: julia 
	 my_tree = load_tree("/home/tvasconcelos/missecomparison/ClaDScomparison/trees/RaboskyEtAl2017_211.tre") 
# mode: julia 
	 output = infer_ClaDS(my_tree, print_state = 100) 
# mode: julia 
	 save_ClaDS_in_R(output, "/home/tvasconcelos/missecomparison/ClaDScomparison/results/RaboskyEtAl2017_211_clads_results.Rdata") 
# mode: julia 
	 my_tree = load_tree("/home/tvasconcelos/missecomparison/ClaDScomparison/trees/RaboskyEtAl2017_221.tre") 
# mode: julia 
	 output = infer_ClaDS(my_tree, print_state = 100) 
# mode: julia 
	 save_ClaDS_in_R(output, "/home/tvasconcelos/missecomparison/ClaDScomparison/results/RaboskyEtAl2017_221_clads_results.Rdata") 
