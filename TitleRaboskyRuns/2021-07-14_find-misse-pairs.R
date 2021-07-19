# Finding adjacent models (hardcode)

find.adjacent.misse <- function(model1, npar_prox=1) {
  if(nrow(model1)==1) {
    rownames(model1) <- NULL
    model1_plus = model2_plus = model3_plus = model1_minus = model2_minus = model3_minus = model1
    if(npar_prox == 1) {  # How many parameters up or down from the set we have now
      model1_plus[1,c(1:1)] <- model1[1,c(1:1)]+1
      model2_plus[1,c(2:2)] <- model1[1,c(2:2)]+1
      #model3_plus[1,c(1:2)] <- model1[1,c(1:2)]+1
      model1_minus[1,c(1:1)] <- model1[1,c(1:1)]-1
      model2_minus[1,c(2:2)] <- model1[1,c(2:2)]-1
      #model3_minus[1,c(1:2)] <- model1[1,c(1:2)]-1
      all_close_models <- rbind(model1_plus, model2_plus, #model3_plus,
                              model1_minus, model2_minus)# , model3_minus)
      all_close_models <- subset(all_close_models, all_close_models$turnover!=0)
      all_close_models <- subset(all_close_models, all_close_models$eps!=0)
    }
  } else { stop("more than one model selected to pair") }
return(all_close_models)
}

# Example 
max.param = 14
possibilities = hisse::generateMiSSEGreedyCombinations(max.param=max.param, 
                                                       vary.both=TRUE, fixed.eps.tries = NA)
test <- possibilities[1,]

find.adjacent.misse(test)


