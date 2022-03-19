source("MiSSE_functions.R") 
base.dir = "/home/tvasconcelos/missecomparison/Simulation_study/TitleRaboskyRuns"
subset_trees1 = c("ClaDS2tree_100_3_2","ClaDS2tree_100_4_1","ClaDS2tree_100_6_2","ClaDS2tree_100_7_4")
trees_subset = load.subset(subset_trees1)
DoSingleRun_new(dir=subset_trees1[3], phy=trees_subset, root_type="madfitz", n.cores=NULL)
