###############################################################################
### Multidimensional Chebyshev Polynomial Approximator                      ###
### Joonhwi Joo                                                             ###
###############################################################################
##### This code runs a Chebyshev polynomial approximation only using the 
##### generic functions embedded in R base package.

##### The function chebyshevApproximate returns the approximated value of the approximator

##### Clear the memory
rm(list=ls())

##### Set the working folder
setwd("C:/Dropbox/2014 Spring/Advanced Quantitave Marketing/Pset 2")

##### Load required libraries

###############################################################################
### 1. Function returning single-dimensional chevyshev polynomial           ###
###############################################################################

##### We write a function which returns the value of Chebyshev polynomial 
##### at the specified value of domain [-1,1], given the degree N
Tn = function(N, x)
{
  return(cos(N*acos(x)))
}

###############################################################################
### 2. Functions approximating multidimensional chevyshev polynomials       ###
###############################################################################
initializeApproximator = function(M,N,a,b,D){
##### We first compute the interpolation nodes
### k is the index vector for nodes
k = c(1:M)

### Now we compute the z vector, which is the vector of nodes in [-1,1]
z = - cos(pi * ((2 * k - 1)/(2 * M)))

##### Next, we create the grid matrix/list
### This is the original z-grid matrix
zmatrix = matrix(z, nrow = length(z), ncol = N)
zlist = split(zmatrix, rep(1:ncol(zmatrix), each = nrow(zmatrix)))
zgridmatrix = t(as.matrix(expand.grid(zlist[1:D])))
zgridlist = split(zgridmatrix, rep(1:ncol(zgridmatrix), each = nrow(zgridmatrix)))
zgridmatrix = t(zgridmatrix)

### The node matrix x has dimension M by D
amatrix = matrix(a, nrow=M, ncol=D, byrow=TRUE)
x = (t(b-a) / 2) %x% (z + 1) + amatrix
xlist = split(x, rep(1:ncol(x), each = nrow(x)))
xgridlist = t(as.matrix(expand.grid(xlist[1:D])))
xgridlist = split(xgridlist, rep(1:ncol(xgridlist), each = nrow(xgridlist)))

### We evaluate f at each point of the grid of x
fevaluated = unlist(lapply(xgridlist, f))

### we take i as the index matrix for dimensions
### and figure out the indexlist
i = matrix(seq(from = 0, to = N), nrow = (N+1), ncol = D)
ilist = split(i, rep(1:ncol(i), each = nrow(i)))
igridmatrix = t(as.matrix(expand.grid(ilist[1:D])))
igridlist = split(igridmatrix, rep(1:ncol(igridmatrix), each = nrow(igridmatrix)))
igridmatrix = t(igridmatrix)

##### We compute the values of coefficients

##### The Numerator part
##### We define a function which returns the numerators
##### which takes the combination of indices i as input
numerator = function(i, x)
{
  i = unlist(i)
  result = fevaluated
  for(j in 1:length(i))
  {
    result = result * Tn(i[j], zgridmatrix[,j])
  }
  return(sum(result))
}

##### The Denominator part
denominator = function(i, x)
{
  sqsum = function(N) sum((cos(N*acos(x)))^2)
  return(prod(sapply(unlist(i), sqsum)))
}


##### Finally save the value of coefficients
coefficients = (sapply(igridlist, numerator, zgridmatrix)) / (sapply(igridlist, denominator, z)) 
result = cbind(igridmatrix, coefficients)

return(result)
}

###############################################################################
### 3. Approximator                                                         ###
###############################################################################
##### The function chebyshevApproximate returns the approximated value of the approximator

chebyshevApproximate = function(evalPoint){
centeredevalPoint = 2 * ((evalPoint - a) / (b - a)) - 1
tempresult = initializeApproximator(M,N,a,b,D)
return(sum(Tn(tempresult[,1], centeredevalPoint[1]) *
             Tn(tempresult[,2], centeredevalPoint[2]) * 
             Tn(tempresult[,3], centeredevalPoint[3]) * tempresult[,4]))
}


###############################################################################
### 4. Parameters                                                           ###
###############################################################################
##### We take the degree N, the (per dimension) number of nodes M, 
##### and the rectangle (a,b) as input. D is automatically determined by (a,b)

### The per-dimension number of nodes M
M = 7
### The degree of approximator N
N = 3
### The box
a = c(-5,-2,-3)
b = c(2,4,3)
### The dimension D
D = length(a)

##### We take the function input here.
f = function(x) {
  return(x[1] * (x[3])^3 + x[2] * x[3] + (x[1])^2 * x[2] * (x[3])^2)
}

###############################################################################
### 5. Test and Evaluation                                                  ###
###############################################################################
### Evaluation point
evalPoint = t(cbind(runif(10000, min = a[1], max = b[1]), runif(10000, min = a[2], max = b[2]), 
                  runif(10000, min = a[3], max = b[3])))
evalPoint = split(evalPoint, rep(1:ncol(evalPoint), each = nrow(evalPoint)))

absDifference = function(evalPoint) {abs(chebyshevApproximate(evalPoint) - f(evalPoint))}

Mlist = c(7,15,30)
Nlist = c(3,5,7)

for(k in 1:3)
{
  N = Nlist[k]
  for(l in 1:3)
  {
    M = Mlist[l]
    difference = sapply(evalPoint, absDifference)
    temp = rbind(system.time(sapply(evalPoint, absDifference), gcFirst=TRUE), max(difference), mean(difference))
    write.csv(temp,paste(c("f-","M-",M,"N-",N,".csv"), collapse=""))
  }
}


##### for g
f = function(x) {
  return(x[1] * log(5 + x[2] * x[3]))
}

for(k in 1:3)
{
  N = Nlist[k]
  for(l in 1:3)
  {
    M = Mlist[l]
    difference = sapply(evalPoint, absDifference)
    temp = rbind(system.time(sapply(evalPoint, absDifference), gcFirst=TRUE), max(difference), mean(difference))
    write.csv(temp,paste(c("g-","M-",M,"N-",N,".csv"), collapse=""))
  }
}

##### for h
f = function(x) {
  return((x[1]) ^ 2 * cos(x[2]) * exp(x[3]))
}

for(k in 1:3)
{
  N = Nlist[k]
  for(l in 1:3)
  {
    M = Mlist[l]
    difference = sapply(evalPoint, absDifference)
    temp = rbind(system.time(sapply(evalPoint, absDifference), gcFirst=TRUE), max(difference), mean(difference))
    write.csv(temp,paste(c("h-","M-",M,"N-",N,".csv"), collapse=""))
  }
}