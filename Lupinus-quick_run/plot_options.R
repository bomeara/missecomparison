
#########################################################################
# Plot options
#########################################################################

# First, let's reload the results of the run that we saved in the previous part
load("Lupinus_example.Rsave")
#

# Some aditional options for ploting MiSSE tip rates are given as follows:
# the file "plot_functions_for_MiSSE" can be found in the supplementary of
# Vasconcelos et al.

# Some plots:

# Let's see the improvement in AIC and lokLik with increasing in model complexity
plot(x=1:length(model.set_pruned), 
     y=unlist(lapply(model.set_pruned[1:length(model.set_pruned)], "[[", "loglik")), 
     xlab="comb", ylab="loglik", pch=19, col="blue")

plot(x=1:length(model.set_pruned), 
     y=unlist(lapply(model.set_pruned[1:length(model.set_pruned)], "[[", "AICc")), 
     xlab="comb", ylab="AICc", pch=19, col="red")

# Where the X-axis corresponds to the models in the table:


boxplot(tip.rates$turnover, tip.rates$net.div, 
        tip.rates$speciation,tip.rates$extinct.frac, 
        tip.rates$extinction, main="tip rates run 1", 
        names=c("turnover", "netdiv", "speciation", "eps", "extinction"))

painted.tree <- hisse::plot.misse.states(x = model.recons, 
                                         rate.param = "speciation", type = "phylo") 

