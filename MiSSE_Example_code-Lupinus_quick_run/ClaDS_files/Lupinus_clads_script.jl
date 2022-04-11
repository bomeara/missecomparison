# mode: julia 
	 using PANDA 
# mode: julia 
	 using JLD2 
# mode: julia 
	 my_tree = load_tree("/home/tvasconcelos/missecomparison/MiSSE_Example_code-Lupinus_quick_run/Lupinus.tre") 
# mode: julia 
	 output = infer_ClaDS(my_tree, print_state = 100) 
# mode: julia 
	 save_ClaDS_in_R(output, "/home/tvasconcelos/missecomparison/MiSSE_Example_code-Lupinus_quick_run/Lupinus_clads_results.Rdata") 
