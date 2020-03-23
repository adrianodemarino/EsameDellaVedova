
### define the grid
nRows <- 6
nCols <- 6
nCells <- nRows*nCols
cellMapping <- matrix(data = NA, nrow = nCells, ncol = 2)
colnames(cellMapping) <- c("x","y")
cellMapping[,1] <- rep(1:(nRows), each = nCols)
cellMapping[,2] <- rep(1:(nCols), times = nRows)
rownames(cellMapping) <- 1:(nCells)

### define the corresponding costs
TSP <- matrix(data = NA, nrow = nCells, ncol = nCells)
rownames(TSP) <- 1:nCells
colnames(TSP) <- 1:nCells
for (rr in 1:nrow(TSP)) {
  for (cc in 1:ncol(TSP)) {
    TSP[rr,cc] <- norm2(cellMapping[rr,],cellMapping[cc,])
  }
}
diag(TSP) <- 9999 #+Inf

minCost <- + Inf
maxIt <- 30
nParents <- 4
nGen <- 4

### initialization of a population
nTot <- 10^4
offspring <- matrix(data = NA, nrow = nTot, ncol = nCells)
for (nPar in 1:nTot) {
  offspring[nPar,] <- sample(x = 1:nCells, nCells, replace = FALSE)
}
nPropSel <- nParents
offspring <- propSelection(population = offspring, n = nPropSel)
# apply(offspring, 1, computeRouteCost)

### GA algorithm
nIt <- 0
EPS <- .01
eps <- +Inf
oldCost <- minCost
while (nIt < maxIt & eps > EPS & oldCost > (nCells - 1)) {
  offspring <- ga(nGen, nParents, offspring, minCost = +Inf, selFrac = .5, mutRate = .1)
  currCost <- computeRouteCost(offspring[1,])
  eps <- abs(oldCost - currCost)
  oldCost <- currCost
  message("\t", currCost)
  nIt <- nIt + 1
}
res <- offspring[1,]
minCost <- computeRouteCost(res)

minCost - (nCells - 1) ## distance from the optimum

## plot results
library(igraph)
TSP_res <- TSP
TSP_res[,] <- NA

for(z in 1:nCells){ 
  if(z < nCells){
    row1 <- res[z]
    col1 <- res[z+1]
    TSP_res[row1,col1] <- TSP[row1,col1] 
  } 
}

g <- graph_from_adjacency_matrix(adjmatrix = TSP_res, mode = "directed")
g <- simplify(g, remove.loops = TRUE)
g <- set_edge_attr(g, "curved", value = .2)
plot(g, layout = layout_on_grid, cex = 0, curved = TRUE,
     edge.width = 7, vertex.size = 25, edge.arrow.width = 3, edge.arrow.size = 1, edge.color = "black", 
     vertex.label.font = 2, vertex.label.cex = 2, vertex.label.color = "black")


#################
### functions ###
#################


##########
### ga ###
##########
##
## genetic algorithm function
##
ga <- function(nGen, nParents, offspring, minCost, selFrac, mutRate){

  gen <- 1

  while (gen <= nGen) { ## loop over generations
    
    message("gen: ", gen, "\tminCost = ", minCost)
    pa <- offspring
    offspring <- cyclicCrossOver(pa) ## compute new offspring
    
    dupOff <- offspring[duplicated(apply(offspring, 1, computeRouteCost)),]
    if(nrow(matrix(dupOff, ncol = nCells)) > 0){
      offspring[duplicated(apply(offspring, 1, computeRouteCost)),] <- t(apply(matrix(dupOff, ncol = nCells),
                                                                               1, FUN = function(x){mutation(x, mutationRate = mutRate)})) ## apply mutation
    }

    paAndOff <- rbind(pa, offspring)
    costs <- apply(paAndOff, 1, computeRouteCost)
    
    currMinCost <- min(costs)
    if(currMinCost <= minCost){
      res <- paAndOff[which(costs == currMinCost)[1],]
      minCost <-  currMinCost
    }
    
    nPropSel <- round(nrow(offspring)*selFrac)
    offspring <- propSelection(population = offspring, n = nPropSel) ## apply selection
    
    gen <- gen + 1
    
  }
  
  return(propSelection(population = offspring, n = nParents)) ## return last offspring (only best nParents elements)
}


#############
### norm2 ###
#############
##
## euclidean norm
##
norm2 <- function(x,y){
  return(sqrt(sum((x - y)^2)))
}

#####################
### propSelection ###
#####################
##
## proportional selection function
##
propSelection <- function(population, n){
  fitness <- apply(population, 1, computeRouteCost)
  names(fitness) <- 1:length(fitness)
  subPop <- population[head(as.numeric(names(sort(fitness, decreasing = FALSE))), n), ]
  
  return(subPop)
}


################
### mutation ###
################
##
## mutation function
##
mutation <- function(x, mutationRate){
  nCells <- length(x)
  
  for (i in 1:nCells) {
    if(runif(1,0,1) < mutationRate){
      j <- sample(x = setdiff(1:nCells, i), size = 1, replace = FALSE)
      
      temp <- x
      temp[i] <- temp[j]
      temp[j] <- x[i]
      x <- temp
    }
  }
  return(x)
}

#######################
### cyclicCrossOver ###
#######################
##
## cyclic crossover function
##
cyclicCrossOver <- function(pa){
  nParents <- nrow(pa)
  nCells <- ncol(pa)
  offspring <- matrix(data = NA, nrow = nParents*(nParents - 1), ncol = nCells)
  
  nOff <- 1
  for (i in 1:nParents) {
    for (j in i:nParents) {
      if(j != i){
        pi <- pa[i,]
        pj <- pa[j,]
        offspring[nOff,] <- crossOver(A = pi, B = pj)
        offspring[nOff + nParents*(nParents - 1)/2,] <- crossOver(A = pj, B = pi)
        
        nOff <- nOff + 1
      }
    }
  }
  
  return(offspring)
}

##################
### cross over ###
##################
##
## crossover function
##
crossOver <- function(A,B){
  os <- rep(NA, length(A))
  terminated <- FALSE
  
  currCol <- 1
  os[1] <- curr <- A[currCol]
  
  terminated <- A[which(B == curr)] == os[1]
  while (!terminated) {
    curr <- os[which(B == curr)] <- A[which(B == curr)]
    terminated <- A[which(B == curr)] == os[1]
  }
  if(sum(is.na(os)) > 0){
    os[is.na(os)] <- B[is.na(os)]
  }
  
  return(os)
}

########################
### computeRouteCost ###
########################
##
## compute the cost of route
##
computeRouteCost <- function(route){
  nCells <- length(route)
  
  cost <- 0
  for(z in 1:nCells){ 
    if(z < nCells){
      row1 <- route[z]
      col1 <- route[z+1]
      cost <- cost + TSP[row1,col1]
    } 
  }
  
  return(cost)
}

