# Part 1 ------------------------------------------------------------------

library(clrdag)

# Part 2 ------------------------------------------------------------------

set.seed(2018)
p <- 25; n <- 2000;sparsity <- 2/p
## generate a random lower triangular adjacnecy matrix
A <- matrix(rbinom(p*p, 1, sparsity)*sign(runif(p*p, -1, 1))*runif(p*p, 0.7, 1), p, p)
A[upper.tri(A, diag=TRUE)] <- 0

## num of edges in A
sum(A != 0) # 26
## data matrix
X <- matrix(rnorm(n*p), n, p) %*% t(solve(diag(p) - A))

'''
#Loglikelihood
loglikelihood=function(A,sigma_square,X){
  n=ncol(X)
  p=nrow(X)
  tot=0
  for (j in 1:n){
    somma2=0
    for (i in 1:p){
      somma1=0
      for (k in 1:n){
        if (k!=j){
          somma1=somma1+A[j,k]*X[j,k]
        }
      }
      partial=X[i,j]-somma1
    }
    somma2=somma2+partial^2
  }
  tot=tot+(1/(2*sigma_square)*somma2+n*log(sigma_square)/2)
  return(-tot)
}
'''


# Part 3 ------------------------------------------------------------------

compute.middle = function(X,A,j, p){
  # Finding the indexes K = {k s.t. k is in [1,p]\j}
  K = which((1:p)!=j)
  # in k we have the indexes of the row without j
  
  inners = apply(X[, K] * A[j, K], 1, sum) # n length vector containing the inner sums for each i-th sample
  middle = sum((X[, j] - inners)^2)  # as a matter of fact in the middle we are working on the j-th column
  return(middle)
}

sigma.estimator = function(X, A){
  n = dim(X)[1]
  p = dim(X)[2]
  return(sum(sapply(1:p, function(j)compute.middle(X,A,j, p)))/(n*p))
}

log_likelihood = function(X, A, sigma.square){
  n = dim(X)[1]
  p = dim(X)[2]
  return(-sum(sapply(1:p, function(j) (compute.middle(X,A,j,p)/(2*sigma.square)+(n/2*log(sigma.square))))))
}


inner_LRT_function.paths = function(X, X.tr, X.te, D, mu, alpha ){
  
  # Computing all the A.h0 over the sparsity parameters 'k'
  A.mles = lapply(1:sum(D), function(mu)MLEdag(X, D=D, tau=alpha, mu=mu, rho=1.2, trace_obj = F))
  A.h0s  = lapply(A.mles, function(a)a$A.H0)
  A.h1s  = lapply(A.mles, function(a)a$A.H1)
  
  # For each A.h0 compute the associated sigma
  sigmas.0 = t(sapply(A.h0s, function(a)sigma.estimator(X.tr, a)))
  sigmas.1 = t(sapply(A.h1s, function(a)sigma.estimator(X.te, a)))
  
  # Computing the likelihoods of each pair (A.h0, sigma.h0)
  likelihoods.0 = t(sapply(1:length(sigmas.0), function(idx)log_likelihood(X.tr, A.h0s[[idx]], sigmas.0[idx])))
  likelihoods.1 = t(sapply(1:length(sigmas.1), function(idx)log_likelihood(X.tr, A.h1s[[idx]], sigmas.1[idx])))
  
  # Finding the idx that give us the maximum likelihood
  mle_idx.0 = which(likelihoods.0==max(likelihoods.0))[1] # maximum likelihood over h0
  mle_idx.1 = which(likelihoods.1==max(likelihoods.1))[1] # maximum likelihood over h0
  
  return (likelihoods.1[mle_idx.1]-likelihoods.0[mle_idx.0])
  #return(list(U.n=exp(log_likelihood(X.tr, A, sigma.unconstrained)-likelihoods[mle_idx]), A.h0=A.h0s[[mle_idx]]))
  
}

inner_LRT_function.links = function(X, X.tr, X.te, D, mu, alpha){
  tmp = MLEdag(X = X, D = D, tau = alpha, mu = mu, 
               rho = 1.2, trace_obj = FALSE)
  A.h0 = tmp$A.H0
  A = tmp$A.H1
  
  sigma.h0 = sigma.estimator(X.tr, A.h0)
  sigma.unconstrained = sigma.estimator(X.te, A)
  
  return(log_likelihood(X.tr, A, sigma.unconstrained)-log_likelihood(X.tr, A.h0, sigma.h0))
}

log.LRT = function(X,D, links=T, mu,alpha){
  n=dim(X)[1]
  X.tr = X[1:(n%/%2),]
  X.te = X[(n%/%2+1):n,]
  
  if(links){
    U_n = inner_LRT_function.links(X, X.tr, X.te, D, mu, alpha)
    U_n.swap = inner_LRT_function.links(X, X.te, X.tr, D, mu, alpha)
  }
  else{
    U_n = inner_LRT_function.paths(X, X.tr, X.te, D, mu, alpha)
    U_n.swap = inner_LRT_function.paths(X, X.te, X.tr, D,mu, alpha)
  }
  return(list(links = links, U_n = U_n, W_n = (log((exp(U_n)+exp(U_n.swap))/2))))
}


# SIZE --------------------------------------------------------------------

#proportion of rejections we get by applying our test on M dataset generated from
#a model compatible with H0=A[j,k]=0 for all (j,k) belong to F

### Parameters
set.seed(2018)
p <- 10
n <- 200
M <- 1000

### H0: F = { (p,2) }, and A[F] = 0
D <- matrix(0, p, p)
D[p,2] = 1

### Adjacency Matrix >> Hub
# All connected to 1, NO EDGE between p-2 >> **COMPATIBLE** with H0
A      <- matrix(0, p, p)     
A[, 1] <- sign( runif( p, min = -1, max = 1 ) )
A[1,1] <- 0

# Simulation --------------------------------------------------------------

M=1000
cont=0
alpha=0.05
contU=0
contW=0
U_n <- rep(NA, M)
W_n <- rep(NA, M)
for(i in 1:M) {
  X   <- matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A) )
  temp=log.LRT(X,D, links=T,mu=1.2,alpha=alpha)
  U_n[i]=temp$U_n
  W_n[i]=temp$W_n
  
  if (U_n[i]>log(1/alpha)) contU=contU+1
  if (W_n[i]>log(1/alpha)) contW=contW+1
}

cat('alpha is equal to for U ', contU/M,'alpha is equal to for W',contW/M)

# 0.049 and 0.053

# POWER -------------------------------------------------------------------

#Beta= the probability of failing to reject the null hypothesis when the null 
# hypothesis is false
p=10
n=100
sparsity=2/p
### H0: F = { (p,1) }, and A[F] = 0
D <- matrix(0, p, p)
D[2:p,1]=1
D[3:p,2]=1
D[2,8=1]
D[5,6]=1
D[8,9]=1

### Adjacency Matrix >> Hub
# All connected to 1, NO EDGE between p-2 >> **NOT COMPATIBLE** with H0
#A      <- matrix(0, p, p)     
#A[2:p, 1] <- sign( runif( p-1, min = -1, max = 1 ) )
set.seed(2018)
A=matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[2,8]=1
A[5,6]=1
A[8,9]=1
A[2:p,1]=1
A[3:p,2]=1


# SIMULATION
contU=0
contW=0
alpha=0.05
M=1000
U_n <- rep(NA, M)
W_n <- rep(NA, M)
for(i in 1:M) {
  X   <- matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A) )
  temp=log.LRT(X,D, links=T,mu=1.2, alpha=alpha)
  U_n[i]=temp$U_n
  W_n[i]=temp$W_n
  
  if (U_n[i]>log(1/alpha)) contU=contU+1
  if (W_n[i]>log(1/alpha)) contW=contW+1
}

cat('1-beta is equal to for U: ', contU/M,'   ,1-beta is equal to for W:',contW/M)
#0.741 and 0.859


# SIZE pathway ------------------------------------------------------------

p = 10
n = 100
sparsity = 2/p

# Building the linkages set
D.pathway = matrix(0, p,p)
D.pathway[3, 4]=1
D.pathway[4, 5]=1
D.pathway[5, 6]=1
D.pathway[6, 7]=1
# Building A accordingly to h0
A.pathway.hidden=matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A.pathway.hidden[3,4]=1
A.pathway.hidden[4,5]=0
A.pathway.hidden[5,6]=0
A.pathway.hidden[6,7]=0

out_U = 0
out_W = 0
m=1000

for(i in 1:m){
  set.seed(2018)
  X = matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A.pathway.hidden) )
  LRT = log.LRT(X, D.pathway, links=F,mu=1.2,alpha=0.05)
  if(LRT$links)
    cat('ERRORE, chiamata LRT con links=', links,'\n')
  out_U = out_U+(LRT$U_n>-log(alpha))
  out_W = out_W+(LRT$W_n>-log(alpha))
}

cat("U_n =", out_U/m, "   W_n = ",out_W/m)
#0 and 0

# POWER pathway -----------------------------------------------------------

p = 10
n= 100
sparsity = 2/p

D.pathway = matrix(0, p,p)
D.pathway[3, 4]=1
D.pathway[4, 5]=1
D.pathway[5, 6]=1

# Building A  NOT accordingly to h0
A.pathway.hidden=matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A.pathway.hidden[3,4]=1
A.pathway.hidden[4,5]=1
A.pathway.hidden[5,6]=1


out_U = 0
out_W = 0
m=1000

for(i in 1:m){
  set.seed(2018)
  X = matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A.pathway.hidden) )
  LRT = log.LRT(X, D.pathway, links=F,mu=1.2,alpha=0.05)
  if(LRT$links)
    cat('ERRORE, chiamata LRT con links=', links,'\n')
  out_U = out_U+(LRT$U_n>-log(alpha))
  out_W = out_W+(LRT$W_n>-log(alpha))
}

cat("U_n =", out_U/m, "   W_n = ",out_W/m)
#1 e 1

# PART 4 ------------------------------------------------------------------

#import concatenated
#data=read.csv("data/cell_signaling/concatenated.csv",header=T,sep=",")
#X=data.matrix(data, rownames.force = NA)
#data will be our X matrix
#p=11

# Let's check the missed connections with the linkage type hypothesis
# So, we are working with a matrix NOT COMPATIBLE with H0:
# we already know that we have some problems in this case 

#FIRST CONNECTION: PIP2->PKC (columns 4  and column 9)
### H0: F = { (4,9) }, and A[F] = 0
'''
D <- matrix(0, p, p)
D[4,9] = 1
alpha=0.05

temp=log.LRT(X,D, links=T,mu=1.2)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
## OKKKK, we accept the null hypothesis
'''

#Let's continue
#SECOND CONNECTION: Plc-gamma -> PKC (column 3 and 9)
### H0: F = { (3,9) }, and A[F] = 0
D <- matrix(0, p, p)
D[3,9] = 1
alpha=0.05

'''
temp=log.LRT(X,D, links=T)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
## OKKKK
'''

#Let's continue
#THIRD CONNECTION: PIP3-> AKT (column 5 and 7)
### H0: F = { (4,7) }, and A[F] = 0
D <- matrix(0, p, p)
D[4,7] = 1
alpha=0.05

'''
temp=log.LRT(X,D, links=T)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
## OKKKK
'''


#PATHWAY LINKAGE
#check all the connections from PKC
#PKC->(RAF;MEK,ERK) which means columns (9->1,2,6)
D <- matrix(0, p, p)
D[9,1] = 1
D[1,2] = 1
D[2,6] = 1
alpha=0.05
'''
temp=log.LRT(X,D, links=F)
U_n=temp$U_n
W_n=temp$W_n
cat("U_n:", U_n, "   W_n:", W_n)

##GOOD
'''

# PART 5 FINAL ------------------------------------------------------------
##
#PROVE BOXCOX
library(MASS)
data=read.csv("data/cell_signaling/normalized_and_scaled/pma_normalized.csv",header=T,sep=",")
X=data.matrix(data, rownames.force = NA)
#data will be our X matrix
p=11

box_cox=function(column){
  bc=boxcox(column~1,lambda=seq(-5,5));
  best_lambda=bc$x[which(bc$y==max(bc$y))]
  if (best_lambda==0) {
    column=sapply(column,function(x) log(x))
  }
  else {
    column=sapply(column,function(x) (x^best_lambda-1)/best_lambda )
  }
  return(column)
}

normalize=function(column){
   column=scale(box_cox(column))
   return(column)
 }


#DATA TRANSFORMATION
for ( i in 1:dim(X)[2]){
  X[,i]=normalize(X[,i])
}

hist(X[,6],breaks=50)


#FIRST CONNECTION: PIP2->PKC (columns 4  and column 9) ### MISSED: ci dovrebbe essere
### H0: F = { (4,9) }, and A[F] = 0
D <- matrix(0, p, p)
D[4,9] = 1
alpha=0.05
MU=100

temp=log.LRT(X,D, links=T,mu=MU,alpha=alpha)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
## OKKKK, we accept the null hypothesis

out=MLEdag(X=X,D=D,tau=alpha, mu=MU, rho=1.2, trace_obj = F)
out$pval
# they reject the null hypothesis
## DIFFERENCE
#first connection
tab <- matrix(c(1,-16.19,-15.48,
                1,-16.19,-15.48,
                0.21,-18.48,-14.64,
                0.21,-18.48,-14.64), ncol=3, byrow=TRUE)
colnames(tab) <- c('pvalue_MLE_dag','Un','Wn')
rownames(tab) <- c('k=0.01','k=0.1','k=1','k=10')
tab <- as.table(tab)



#SECOND CONNECTION: PKC->Jnk  # C'é
### H0: F = { (9,11) }, and A[F] = 0
D <- matrix(0, p, p)
D[9,11] = 1
alpha=0.05
MU=0.01

U_n=rep(NA,41)
W_n=rep(NA,41)
pvalues=rep(NA,41)
a=seq(0,2,0.05)
j=1
for (i in a){
  temp=log.LRT(X,D, links=T,mu=i, alpha=alpha)
  U_n[j]=temp$U_n
  W_n[j]=temp$W_n
  out=MLEdag(X,D=D,tau=alpha, mu=i, rho=1.2, trace_obj = F)
  pvalues[j]=out$pval
  j=j+1
}

par(mfrow=c(1,1),mar=c(5,2,4.1,2.1))
# plot the first curve by calling plot() function
# First curve is plotted
plot(a, pvalues, type="o", col="blue", pch="o", ylim=c(-20,20),lty=1,main="Statistics" )

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(a, U_n, col="red", pch="*",ylim=c(-20,20))
lines(a, U_n, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(a, W_n, col="green",ylim=c(-20,20),pch="+")
lines(a, W_n, col="green", lty=3)

# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend("bottomright",legend=c("pvalues","U_n","W_n"), col=c("blue","red","green"),
       pch=c("o","*","+"),lty=c(1,2,3),cex=0.5, ncol=1)

'''
temp=log.LRT(X,D, links=T,mu=MU, alpha=alpha)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
#dovremmo rigettare
#male
out=MLEdag(X,D=D,tau=alpha, mu=MU, rho=1.2, trace_obj = F)
out$pval
'''
tab <- matrix(c(1,-16.19,-15.48,
                ,1,-16.19,-15.48
                ,6.0e-11,17.05,16.36
                ,6.0e-11,17.05,16.36,exp(-3.189)), ncol=3, byrow=TRUE)
colnames(tab) <- c('pvalue_MLE_dag','Un','Wn')
rownames(tab) <- c('k=0.01','k=0.1','k=1','k=10')
tab <- as.table(tab)


#Third CONNECTION: PKC->PIP3 (columns 4  and column 9) # NON C'è
### H0: F = { (9,5) }, and A[F] = 0
D <- matrix(0, p, p)
D[9,5] = 1
alpha=0.05
#U_n=rep(NA,100)
#W_n=rep(NA,100)
#pvalues=rep(NA,100)
MU=100

#for ( i in 1:100){
#MU=i
temp=log.LRT(X,D, links=T,mu=MU,alpha=alpha)
U_n=temp$U_n
W_n=temp$W_n

#if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
#dovremmo accettare H0
#ok
out=MLEdag(X,D=D,tau=alpha, mu=MU, rho=1.2, trace_obj = F)
out$pval
  
#}


tab <- matrix(c(1,-16.19,-15.48
                ,0.60,-15.80,-14.80
                ,0.60,-15.80,-14.80
                ,0.60,-15.80,-14.80), ncol=3, byrow=TRUE)
colnames(tab) <- c('pvalue_MLE_dag','Un','Wn')
rownames(tab) <- c('k=0.01','k=0.1','k=1','k=10')
tab <- as.table(tab)

#PATHWAY
#PATHWAY LINKAGE
#check all the connections from PKC
#PKC->(RAF;MEK,ERK) which means columns (9->1,2,6)
D <- matrix(0, p, p)
D[9,1] = 1
D[1,2] = 1
D[2,6] = 1
alpha=0.05
MU=0.1
#first connection
'''
temp=log.LRT(X,D, links=F,mu=MU,alpha=alpha)
U_n=temp$U_n
W_n=temp$W_n
cat("U_n:", U_n, "   W_n:", W_n)

out=MLEdag(X,D=D,tau=alpha, mu=MU, rho=1.2, trace_obj = F)
out$pval
'''
U_n=rep(NA,41)
W_n=rep(NA,41)
pvalues=rep(NA,51)
a=seq(0,5,0.1)
j=1
for (i in a){
  #temp=log.LRT(X,D, links=F,mu=i, alpha=alpha)
  #U_n[j]=temp$U_n
  #W_n[j]=temp$W_n
  out=MLEdag(X,D=D,tau=alpha, mu=i, rho=1.2, trace_obj = F)
  pvalues[j]=out$pval
  j=j+1
}

par(mfrow=c(1,1),mar=c(5,2,4.1,2.1))
# plot the first curve by calling plot() function
# First curve is plotted
plot(a, pvalues, type="o", col="blue", pch="o", ylim=c(-0,1),lty=1,main="pvalues")

'''
# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(a, U_n, col="red", pch="*",ylim=c(-100,20))
lines(a, U_n, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(a, W_n, col="green",ylim=c(-100,20),pch="+")
lines(a, W_n, col="green", lty=3)
'''
# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend("bottomleft",legend=c("pvalues"), col=c("blue"),
       pch=c("o"),lty=c(1),cex=0.4, ncol=1)


tab <- matrix(c(0.82,-71.56,-56.09
                ,1,-71.56,-56.09
                ,1.62e-106,-71.56,-56.09
                ,1.62e-106,-71.56,-56.09), ncol=3, byrow=TRUE)
colnames(tab) <- c('pvalue_MLE_dag','Un','Wn')
rownames(tab) <- c('k=0.01','k=0.1','k=1','k=10')
tab <- as.table(tab)


# PART6 -------------------------------------------------------------------
#import concatenated
data=read.csv("data/cell_signaling/normalized_and_scaled/concatenated_normalized.csv",header=T,sep=",")
X=data.matrix(data, rownames.force = NA)
#data will be our X matrixSu
p=11
#transformation

#DATA TRANSFORMATION
for ( i in 1:dim(X)[2]){
  X[,i]=normalize(X[,i])
}

#FIRST CONNECTION: PIP2->PKC (columns 4  and column 9)
### H0: F = { (4,9) }, and A[F] = 0
D <- matrix(0, p, p)
D[4,9] = 1
alpha=0.05
MU=100

temp=log.LRT(X,D, links=T,mu=MU,alpha=alpha)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")

out=MLEdag(X,D=D,tau=alpha, mu=MU, rho=1.2, trace_obj = F)
out$pval

## DIFFERENCE
tab <- matrix(c(1,-625.26,-534.54
                ,1,-625.26,-534.54
                ,1,-466.80,-407.60
                ,1,-466.80,-407.60, ncol=3, byrow=TRUE)
colnames(tab) <- c('pvalue_MLE_dag','Un','Wn')
rownames(tab) <- c('k=0.01','k=0.1','k=1','k=10')
tab <- as.table(tab)



#SECOND CONNECTION: PKC->Jnk (columns 4  and column 9)
### H0: F = { (9,11) }, and A[F] = 0
D <- matrix(0, p, p)
D[9,11] = 1
alpha=0.05
MU=0.2
temp=log.LRT(X,D, links=T,mu=MU,alpha=alpha)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
#dovremmo rigettare
#male
out=MLEdag(X,D=D,tau=alpha, mu=MU, rho=1.2, trace_obj = F)
out$pval


U_n=rep(NA,51)
W_n=rep(NA,51)
pvalues=rep(NA,51)
a=seq(0,10,0.2)
j=1
for (i in a){
  temp=log.LRT(X,D, links=T,mu=i, alpha=alpha)
  U_n[j]=temp$U_n
  W_n[j]=temp$W_n
  out=MLEdag(X,D=D,tau=alpha, mu=i, rho=1.2, trace_obj = F)
  pvalues[j]=out$pval
  j=j+1
}

par(mfrow = c(3,1),mar=c(2,2,1,1))
# plot the first curve by calling plot() function
# First curve is plotted
plot(a, pvalues, type="o", col="blue", pch="o", ylim=c(-0.5,1.5),main="pvalues" )

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
plot(a, U_n, col="red",main="U_n",ylim=c(-500,500))
#points(a, U_n, col="red", pch="*",ylim=c(-300,100),main="U_n")
lines(a, U_n, col="red")

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
plot(a, W_n, col="green",main="W_n",ylim=c(-500,500))
#points(a, W_n, col="green",ylim=c(-300,100),pch="+")
lines(a, W_n, col="green")





tab <- matrix(c(0,-297.12,-297.81
                ,0,-297.12,-297.81
                ,0,36.50,35.81
                ,0,-297.12,-297.81, ncol=3, byrow=TRUE)
colnames(tab) <- c('pvalue_MLE_dag','Un','Wn')
rownames(tab) <- c('k=0.01','k=0.1','k=1','k=10')
tab <- as.table(tab)


#Third CONNECTION: PKC->PIP3 (columns 4  and column 9)
### H0: F = { (9,5) }, and A[F] = 0
D <- matrix(0, p, p)
D[9,5] = 1
alpha=0.05
MU=100

temp=log.LRT(X,D, links=T,mu=MU,alpha=alpha)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
#dovremmo accettare H0
#ok
out=MLEdag(X,D=D,tau=alpha, mu=MU, rho=1.2, trace_obj = F)
out$pval


tab <- matrix(c(1,-625.26,-534.54
                ,1,-625.26,-534.54
                ,1,-466.80,-407.60
                ,0.84,-585.54,-504.55), ncol=3, byrow=TRUE)
colnames(tab) <- c('pvalue_MLE_dag','Un','Wn')
rownames(tab) <- c('k=0.01','k=0.1','k=1','k=10')
tab <- as.table(tab)

#PATHWAY
#PATHWAY LINKAGE
#check all the connections from PKC
#PKC->(RAF;MEK,ERK) which means columns (9->1,2,6)
D <- matrix(0, p, p)
D[9,1] = 1
D[1,2] = 1
D[2,6] = 1
alpha=0.05
MU=100


temp=log.LRT(X,D, links=F,mu=MU,alpha=alpha)
U_n=temp$U_n
W_n=temp$W_n
cat("U_n:", U_n, "   W_n:", W_n)
# non rigettiamo l'hp nulla quindi ok
#out=MLEdag(X,D=D,tau=alpha, mu=MU, rho=1.2, trace_obj = F)
#out$pval



U_n=rep(NA,41)
W_n=rep(NA,41)
pvalues=rep(NA,41)
a=seq(0,100,1)
j=1
for (i in a){
  #temp=log.LRT(X,D, links=F,mu=i, alpha=alpha)
  #U_n[j]=temp$U_n
  #W_n[j]=temp$W_n
  out=MLEdag(X,D=D,tau=alpha, mu=i, rho=1.2, trace_obj = F)
  pvalues[j]=out$pval
  j=j+1
}

par(mfrow=c(1,1),mar=c(5,2,4.1,2.1))
# plot the first curve by calling plot() function
# First curve is plotted
plot(a, pvalues, type="l", col="blue", pch="o", ylim=c(0,1),lty=1 ,main="pvalues")

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(a, U_n, col="red", pch="*",ylim=c(-100,20))
lines(a, U_n, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(a, W_n, col="green",ylim=c(-100,20),pch="+")
lines(a, W_n, col="green", lty=3)

# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend("bottomleft",legend=c("pvalues"), col=c("blue"),
       pch=c("o"),lty=c(1),cex=0.4, ncol=1)


tab <- matrix(c(1,-924.26,-Inf
                ,1,-924.26,-Inf
                ,0,-924.26,-Inf
                ,0,-924.26,-Inf), ncol=3, byrow=TRUE)
colnames(tab) <- c('pvalue_MLE_dag','Un','Wn')
rownames(tab) <- c('k=0.01','k=0.1','k=1','k=10')
tab <- as.table(tab)


### rispondere alle domande 
