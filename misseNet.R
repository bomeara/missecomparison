

library(hisse)
library(igraph)


######################################################################################################################################
######################################################################################################################################
### Step 1: Get some data to play around with
######################################################################################################################################
######################################################################################################################################

#Takes a vector of free parameters and find the distance between them and gives back an
#edge matrix for models 1 freepar away from each other
GetEdges <- function(free.pars, nodes){
    distances <- dist(free.pars)
    dist.mat <- as.matrix(distances)
    dist.mat[lower.tri(dist.mat)] <- max(dist.mat)
    one.par.diff <- which(dist.mat == 1, arr.ind=TRUE)
    tmp.edges.mat <- c()
    for(index in 1:dim(one.par.diff)[1]){
        tmp.edges.mat <- rbind(tmp.edges.mat, c(nodes[one.par.diff[index,1]], nodes[one.par.diff[index,2]]))
    }
    final.edge.mat <- data.frame(x=tmp.edges.mat[,1], y=tmp.edges.mat[,2], weight=rep(1,dim(tmp.edges.mat)[1]))
    return(final.edge.mat)
}


#Start with a fake model space
set.seed(42)
model.space <- generateMiSSEGreedyCombinations(turnover.tries=1:4, eps.tries=1:4)

#What follows is sort of what I would imagine the model results table to look like. The weights are not real:
fake.weights <- rlnorm(dim(model.space)[1])
fake.weights <- fake.weights/sum(fake.weights)
model.space <- cbind(model.space[,c(1,2)], fake.weights)


######################################################################################################################################
######################################################################################################################################
### Step 2: Make a list of nodes and edges to make our network
######################################################################################################################################
######################################################################################################################################

#Very basic igraph just requires that you specify two things: 1) nodes and 2) edges.
#We have all the information we need in the table above. Note this is a very generic
#set of code and plot. But its a good starting point.

nodes <- paste("T", model.space[,1], "_", "E", model.space[,2], sep="")
free.pars <- model.space[,1] + model.space[,2]
edges <- GetEdges(free.pars, nodes)
node.size <- setNames(model.space[,3], nodes)

#Now we can put the pieces together for igraph:
df <- graph.data.frame(edges)
V(df)$names <- nodes
plot(df,vertex.label=V(df)$names, edge.arrow.size=0.01, vertex.size=node.size*100)

#Some alternatives:
tkplot(df,vertex.label=V(df)$names, edge.arrow.size=0.01, vertex.size=node.size*100)







