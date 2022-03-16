# This function can be found on Data Dryad, associated with the following publication:
# Beaulieu JM, O'Meara BC (2015) Extinction can be estimated from moderately sized molecular phylogenies. Evolution, 69, 1036–1043.
# https://datadryad.org/resource/doi:10.5061/dryad.33p91


#epsvec <- c(0, .25, .5, .75);
#lamvec <- c(.06437, 0.07834, 0.1018213, 0.1512680);
#sdvec <- c(0.001, seq(0.01, .06, by=0.01));
#REPS <- 2000;

evolveRateTree <- function (b, stdev, time.stop = 0, mintax = 3, maxtax = 5000) 
{
	
	if (time.stop == 0) 
	stop("Must have stopping criterion\n");
    
	breaker <- 0;
	while (1) {
 		edge <- rbind(c(1, 2), c(1, 3));
		
		birth.time <- rep(0, maxtax);
		end.time <- rep(-1, maxtax);
  		lambda <- rep(b, maxtax);
		edge.status <- rep(FALSE, maxtax);
		
		current <- 2;
		next.node <- 4;
		
		currvec <- current;
		
		while (sum(edge.status[1:nrow(edge)]) < nrow(edge)){
			
			current <- edge[,2][edge.status[1:nrow(edge)] == FALSE][1]; #new current
			t <- birth.time[1:nrow(edge)][edge[,2] == current];
			
			repeat{
				
				currvec <- c(currvec, current);
				
				dt <- rexp(1, lambda[1:nrow(edge)][edge[,2] == current]);
				t <- t + dt;
				print(t)
				if (t >= time.stop) {
					end.time[1:nrow(edge)][edge[,2] == current] <- time.stop;
					edge.status[1:nrow(edge)][edge[,2] == current] <- TRUE;
#cat('Tmax exceeded:', t, 'curr:', current, '\n', sep='\t');
#print((1:nrow(edge))[edge[,2] == current]);
#print(edge.status[edge[,2] == current]); 
#print(edge.status[(1:nrow(edge))[edge[,2] == current]]);
#print(edge.status[1:nrow(edge)][edge[,2] == current]) 
					break;
				}
##
				
				edge.status[1:nrow(edge)][edge[,2] == current] <- TRUE;
				
				edge <- rbind(edge, c(current, next.node), c(current, next.node + 1));
				
				birth.time[1:nrow(edge)][edge[,1] == current] <- t; 
				
				end.time[1:nrow(edge)][edge[,2] == current] <- t;
				
## evolving lambda:
				
				lambda[1:nrow(edge)][edge[,1] == current] <- rlnorm(1, log(lambda[1:nrow(edge)][edge[,2] == current]), sqrt(dt)*stdev);
				
				
				
#print(end.time[1:nrow(edge)][edge[,2 ] == current]);
				
#cat('curr:', current, 'next:', next.node, 'time:', t, '\n', sep='\t');
#z <- data.frame(v1=edge[,1], v2=edge[,2], es=edge.status[1:nrow(edge)], bt=birth.time[1:nrow(edge)], et=end.time[1:nrow(edge)]);
#print(z);
				
				
				current <- next.node;
				next.node <- next.node + 2;
				
				if (nrow(edge) >= maxtax)
				break;
				
#breaker <- breaker + 1;
#if (breaker > 5000)
#	stop('breaker exceeded\n');
				
			}#repeat
			
			if (nrow(edge) >= maxtax){
				print('maxtax exceeded');
				break;	
			}
			
			
		}#while (sum(edge.status))
		
		birth.time <- birth.time[1:nrow(edge)];
		end.time <- end.time[1:nrow(edge)];
		lambda <- lambda[1:nrow(edge)];
        
		edge.length <- end.time - birth.time;
		
		if (nrow(edge) >= 2*mintax & nrow(edge) < maxtax)
		break;
		
    }# while (1)
    
	
	n <- -1
	for (i in 1:max(edge)) {
		if (any(edge[, 1] == i)) {
			edge[which(edge[, 1] == i), 1] <- n;
			edge[which(edge[, 2] == i), 2] <- n;
			n <- n - 1;
		}
	}
	edge[edge > 0] <- 1:sum(edge > 0);
	tip.label <- 1:sum(edge > 0);
	mode(edge) <- "character";
	mode(tip.label) <- "character";
	obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label, birth.time= birth.time, end.time = end.time, lambda=lambda, currvec = currvec);
    class(obj) <- "phylo";
    obj <- old2new.phylo(obj);
    obj;
    
}




# This evolves a tree with extinction, holding the ratio mu/lambda constant, but allowing the 
# speciation rate to evolve. 
# parameter b is the speciation rate, not the net diversification rate.

evolveRateTree.eps <- function (b, eps, stdev, time.stop=50, mintax = 3, maxtax = 100000, return.all.extinct=FALSE) 
{
	if (time.stop == 0) 
	stop("Must have stopping criterion\n");
    too.many<-c()
	min.lambda<-b
	breaker <- 0;
	while (1) {
 		edge <- rbind(c(1, 2), c(1, 3));
		
		birth.time <- rep(0, maxtax);
		end.time <- rep(-1, maxtax);
		max.time=time.stop
		lambda <- rep(b, maxtax);
  		mu<- rep(b*eps, maxtax);
  		dt.check<-0
		#cat(lambda[1], mu[1], '\n');
  		
		edge.status <- rep(FALSE, maxtax);
		
		current <- 2;
		next.node <- 4;
		
		currvec <- current;
		
		while (sum(edge.status[1:nrow(edge)]) < nrow(edge)){
			current <- edge[,2][edge.status[1:nrow(edge)] == FALSE][1]; #new current
			t <- birth.time[1:nrow(edge)][edge[,2] == current];
			
#			print(max(unique(lambda)))
			index <- 0;
			repeat{
				index <- index+1; # indicator...
				
				currvec <- c(currvec, current);
				
				lambda.current <- lambda[1:nrow(edge)][edge[,2] == current];
				mu.current <- mu[1:nrow(edge)][edge[,2] == current];
				#cat('lam-mu', lambda.current, mu.current, 'index', index, '\n');
				
				dt <- rexp(1, rate=(lambda.current+mu.current));
				dt.check<-c(dt.check,dt)
				#cat('index', index,t, dt, '\n');
				t <- t + dt;
				# exceed time?
				if (t >= time.stop) {
					end.time[1:nrow(edge)][edge[,2] == current] <- time.stop;
					edge.status[1:nrow(edge)][edge[,2] == current] <- TRUE;
					break;
				}				
				
				#birth or death:
				temp <- runif(1);
				if (temp <= lambda.current /(lambda.current + mu.current)){
					# birth
					edge.status[1:nrow(edge)][edge[,2] == current] <- TRUE;
					edge <- rbind(edge, c(current, next.node), c(current, next.node + 1));
					birth.time[1:nrow(edge)][edge[,1] == current] <- t; 
					end.time[1:nrow(edge)][edge[,2] == current] <- t;
					## evolving lambda (holding eps constant)
					lambda[1:nrow(edge)][edge[,1] == current] <- rlnorm(1, log(lambda[1:nrow(edge)][edge[,2] == current]), sqrt(dt)*stdev);	
					#this evolves the net div rate
					#rtemp <- rlnorm(1, log(lambda[1:nrow(edge)][edge[,2] == current] - mu[1:nrow(edge)][edge[,2] == current]), sqrt(dt)*stdev);
					#lambda[1:nrow(edge)][edge[,1] == current] <- rtemp / (1-eps);
					mu[1:nrow(edge)][edge[,1] == current] <- lambda[1:nrow(edge)][edge[,1] == current]*eps;						
				}else{
					# death
					min.lambda<-c(min.lambda,min(lambda))
					edge.status[1:nrow(edge)][edge[,2] == current] <- TRUE;
					end.time[1:nrow(edge)][edge[,2] == current] <- t;
					break;
					
				}
				#print(end.time[1:nrow(edge)][edge[,2 ] == current]);
				#cat('curr:', current, 'next:', next.node, 'time:', t, '\n', sep='\t');
				#z <- data.frame(v1=edge[,1], v2=edge[,2], es=edge.status[1:nrow(edge)], bt=birth.time[1:nrow(edge)], et=end.time[1:nrow(edge)]);
				#print(z);
				
				current <- next.node;
				next.node <- next.node + 2;
				
				if (nrow(edge) >= maxtax)
				break;
				
			}#repeat
			if (nrow(edge) >= maxtax){
				too.many<-c(too.many,1)
				print('maxtax exceeded');
				print(max(lambda))
				break;	
			}			
		}#while (sum(edge.status))
		birth.time <- birth.time[1:nrow(edge)];
		end.time <- end.time[1:nrow(edge)];
		lambda <- lambda[1:nrow(edge)];
		edge.length <- end.time - birth.time;
		
		if (nrow(edge) >= 2*mintax & nrow(edge) < maxtax)
		break;
	}# while (1)
    
	n <- -1
	for (i in 1:max(edge)) {
		if (any(edge[, 1] == i)) {
			edge[which(edge[, 1] == i), 1] <- n;
			edge[which(edge[, 2] == i), 2] <- n;
			n <- n - 1;
		}
	}
	edge[edge > 0] <- 1:sum(edge > 0);
	tip.label <- 1:sum(edge > 0);
	mode(edge) <- "character";
	mode(tip.label) <- "character";
	obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label, birth.time= birth.time, end.time = end.time, lambda=lambda, dt.check=dt.check, currvec = currvec, too.many = too.many);
    class(obj) <- "phylo";
    obj <- old2new.phylo(obj);
    obj;
    
}

#this works. Number of species blows up very, very quickly with evolving rate
rateShiftTree <- function (b, pb, time.stop = 0, seed = 0, maxtax = 5000, type='decrease', shiftmag=2, mode='no_evolve') 
{
	
	if (time.stop == 0) 
	stop("Must have stopping criterion\n");
    
	breaker <- 0;
	while (1) {
 		edge <- rbind(c(1, 2), c(1, 3));
		
        
		birth.time <- rep(0, maxtax);
		end.time <- rep(-1, maxtax);
  		lambda <- rep(b, maxtax);
		edge.status <- rep(FALSE, maxtax);
		
		current <- 2;
		next.node <- 4;
		
		currvec <- current;
		numberDecreases <- 0;
		numberIncreases <- 0;
		templam <- b;
		tvec <- 0;
		
#while loop makes sure all lineages have run to completion
		while (sum(edge.status[1:nrow(edge)]) < nrow(edge)){
			
			
			current <- edge[,2][edge.status[1:nrow(edge)] == FALSE][1]; #new current
			t <- birth.time[1:nrow(edge)][edge[,2] == current];
			
			repeat{
				
				currvec <- c(currvec, current);
				
				dt <- rexp(1, lambda[1:nrow(edge)][edge[,2] == current]);
				
				t <- t + dt;
				
#templam <- c(templam, lambda[1:nrow(edge)][edge[,2] == current]);
#tvec <- c(tvec, t);
#print(c(lambda[1:nrow(edge)][edge[,2] == current], t));
				
				if (t >= time.stop) {
					end.time[1:nrow(edge)][edge[,2] == current] <- time.stop;
					edge.status[1:nrow(edge)][edge[,2] == current] <- TRUE;
					
					break;
				}
##
				
				edge.status[1:nrow(edge)][edge[,2] == current] <- TRUE;
				
				edge <- rbind(edge, c(current, next.node), c(current, next.node + 1));
				
				birth.time[1:nrow(edge)][edge[,1] == current] <- t; 
				
				end.time[1:nrow(edge)][edge[,2] == current] <- t;
				
				
				
#binomial rate decrease: constant probability of decrease per speciation event
				
				
				event <- rbinom(1, 1, pb);
				
#lambda[1:nrow(edge)][edge[,1] == current] <- rlnorm(1, log(lambda[1:nrow(edge)][edge[,2] == current]), sqrt(dt)*stdev);
				
				if (event){
					if (type == 'decrease'){
						lambda[1:nrow(edge)][edge[,1] == current] <- 0.00001;
						
						numberDecreases <- numberDecreases +1;
					}
					if (type == 'increase'){
						if (mode == 'evolve'){
							lambda[1:nrow(edge)][edge[,1] == current] <- lambda[1:nrow(edge)][edge[,2] == current]*shiftmag;
							numberIncreases <- numberIncreases+1;
						}else if (mode == 'no_evolve' & lambda[1:nrow(edge)][edge[,2] == current] == b){
							lambda[1:nrow(edge)][edge[,1] == current] <- b*shiftmag;
							numberIncreases <- numberIncreases+1;
						} else if (mode == 'no_evolve' & lambda[1:nrow(edge)][edge[,2] == current] == b*shiftmag){
							lambda[1:nrow(edge)][edge[,1] == current] <- b*shiftmag;
						}else{
							stop('error in new rate assignment\n')
						}	  
						
					}
					
				}else{
					lambda[1:nrow(edge)][edge[,1] == current] <- lambda[1:nrow(edge)][edge[,2] == current];
				}
				
				current <- next.node;
				next.node <- next.node + 2;	
				
				
				if (nrow(edge) >= maxtax)
				break;
				
			}#repeat
			
#return(); ### error check ... 
			
			if (nrow(edge) >= maxtax){
				print('maxtax exceeded');
				break;	
			}
			
			
		}#while (sum(edge.status))
		
		birth.time <- birth.time[1:nrow(edge)];
		end.time <- end.time[1:nrow(edge)];
		lambda <- lambda[1:nrow(edge)];
		edge.length <- end.time - birth.time;
		
		if (nrow(edge) > 2 & nrow(edge) < maxtax)
		break;
		
    }# while (1)
    
	
	n <- -1
	for (i in 1:max(edge)) {
		if (any(edge[, 1] == i)) {
			edge[which(edge[, 1] == i), 1] <- n;
			edge[which(edge[, 2] == i), 2] <- n;
			n <- n - 1;
		}
	}
	edge[edge > 0] <- 1:sum(edge > 0);
	tip.label <- 1:sum(edge > 0);
	mode(edge) <- "character";
	mode(tip.label) <- "character";
	obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label, birth.time= birth.time, end.time = end.time, lambda=lambda, currvec = currvec, numberDecreases=numberDecreases, numberIncreases = numberIncreases, templam = templam, tvec = tvec);
    class(obj) <- "phylo";
    obj <- old2new.phylo(obj);
    obj;
    
}




