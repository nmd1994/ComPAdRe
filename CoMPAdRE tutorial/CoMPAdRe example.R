#CoMPAdRe simple example 
#Inputs are X, Y, rho.factor = 0.75, niter, rho = 0.9 

#Dependencies 
library(splines)
library('MASS')
library('gglasso')
library('glmnet')
#library('mSSL')
library('glasso')
library('CVglasso')
library('rlist')
library('doSnow')
library('doParallel')
library('mgcv')
library('gamsel')
library('matrixStats')
library('nlme')
#X is nxp matrix of predictors, Y is nxq matrix of responses, 
#rho.factor is penalty on non-linear seelction, niter is number of iterations


#----------------------------------------------
#STEP 1: SIMULATE DATASETS 
#----------------------------------------------
#Simulate Data, n = 250, p = 10, q = 10, second order dependence rho = 0.9
n=250
q=10
p=10

X=matrix(runif(n*p,-1,1),nrow=n,ncol=p)

#Generate noise in a Toeplitz structure 
rho=0.9
bd <- rho^(abs(outer(1:q, 1:q, "-")))
E=mvrnorm(n,rep(0,q),bd)


Y=matrix(0,nrow=n,ncol=q)

#SELECTION OF NON-LINEAR RELATIONSHIPS - Define functions 
f1=function(x){1*(1-exp(-2*x))}
f2=function(x){1*x^(2)}
f3=function(x){1*x^(3)}
f4=function(x,sigma){1*((1/(sqrt(2*pi)*sigma))*exp(-(x^(2)/(2*sigma^(2)))))}
f5=function(x){1*x}


#SELECTION OF NON-LINEAR FUNCTION POSTIIONS (pxq) - levels of sparsity and positions 
response=seq(1,5,1)
covs=seq(1,p,1)

#Two levels 
#First layer... select which responses not sparse
spar=4
response=sample(response,spar)
response=sort(response)
#Second layer of sparsity, select # of covariates for each chosen response
predictors=sample(c(1:5),spar,replace=TRUE)
cov4response=list(0)
funccov4response=list(0)

#Given selected predictors, sample functions 
for(i in 1:length(response)){
  inte=sample(1:5,predictors[i],prob=c(0.125,0.125,0.125,0.125,0.5),replace=TRUE)
  #SELECT WHICH COVARIATES TO APPLY IT TO 
  ct=sample(1:p,length(inte))
  cov4response[[i]]=ct
  funccov4response[[i]]=inte
  ct=cov4response[[i]]
  inte=funccov4response[[i]]
  
  for(j in 1:length(inte)){
    if(inte[j]==1){
      Y[,response[i]]=Y[,response[i]]+f1(X[,ct[j]])
    }
    
    if(inte[j]==2){
      Y[,response[i]]=Y[,response[i]]+f2(X[,ct[j]])
    }
    
    if(inte[j]==3){
      Y[,response[i]]=Y[,response[i]]+f3(X[,ct[j]])
    }
    
    if(inte[j]==4){
      Y[,response[i]]=Y[,response[i]]+f4(X[,ct[j]],sigma=0.1)
    }
    
    if(inte[j]==5){
      Y[,response[i]]=Y[,response[i]]+f5(X[,ct[j]])
    }
  }
}

Y_truth=Y
Y=Y + E

#Run Compadre for 5 iterations, rho.factor = 0.25
Results=Compadre(X=X,Y=Y,lambda.fac=0.25,niter=5)
