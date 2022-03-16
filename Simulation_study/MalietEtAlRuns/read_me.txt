This folders contain the simulated phylogenies used to check the statistical behavior of ClaDS 
(O. Maliet, F. Hartig and H. Morlon, Nature Ecology and Evolution 2019). More specifically :

	/ClaDS0/ contains the phylogenies used to test the model without extinction, 
ClaDS0, on phylogenies simulated with the ClaDS0 model. The .Rdata files in this folder contain the following objects :
	sigma				the sigma parameter used to simulate the tree
	alpha				the alpha parameter used to simulate the tree
	lambda_0			the initial speciation rate used to simulate the tree  
	speciation_rates	a vector with the simulated branch-specific speciation rates 
	tree				the simulated phylogeny
	
	/ClaDS2/ contains the phylogenies used to test the model with constant turnover rate, 
ClaDS2, on phylogenies simulated with the ClaDS0 model. The .Rdata files in this folder contain the following objects :
	sigma				the sigma parameter used to simulate the tree
	alpha				the alpha parameter used to simulate the tree
	epsilon				the turnover rate used to simulate the tree
	lambda_0			the initial speciation rate used to simulate the tree  
	speciation_rates	a vector with the simulated branch-specific speciation rates 
	extinction_rates	a vector with the simulated branch-specific extinction rates (speciation_rates * epsilon)
	tree				the simulated phylogeny
	
	/OneShift/ contains the phylogenies used to test the model without extinction, 
ClaDS0, on phylogenies simulated with one big rate shift. The .Rdata files in this folder contain the following objects : 
	speciation_rates	a vector with the simulated branch-specific speciation rates 
	tree				the simulated phylogeny