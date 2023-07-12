#NEW SIMULATIONS WITH GAMSEL 3/23/23
library(splines)
library('MASS')
library('gglasso')
library('glmnet')
#library('mSSL')
library('glasso')
library('CVglasso')
library('rlist')
#install.packages('doSnow')
#install.packages('doParallel')
library('mgcv')
library('gamsel')
library('matrixStats')
library('nlme')
#install.packages('madstat')
#library('madstat')
#-----------------------------------------
#Pre-specified functions 
#----------------------------------------
pen_mat <- function(inner_knots) {
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots))
  d <- diff(inner_knots)  # The vector of knot differences; b - a 
  g_ab <- splineDesign(knots, inner_knots, derivs = 2) 
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 
      4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 
}
#--------------------------------------------

#--------------------------------------
#STEP 1: Simulate Data 
#--------------------------------------
library(splines)
library('MASS')
library('gglasso')
library('glmnet')
library('glasso')
library('CVglasso')  
library('gtools')
Bov=list(0)
Btrue=list(0)
Ytrue=list(0)
#Yov=list(0)
#Yovtrue=list(0)
B_TRUE_res=list(0)
B_TRUE_gamsel=list(0)
B_FALSE_res=list(0)
B_TRUE_init=list(0)
B_FALSE_gamsel=list(0)
B_FALSE_init=list(0)
#SEL_TRUE_COMAD_ov=list(0)
#SEL_FALSE_COMAD_ov=list(0)
#SEL_TRUE_gamsel_ov=list(0)
#SEL_FALSE_gamsel_ov=list(0)
SEL_TRUE_init_ov=list(0)
SEL_FALSE_init_ov=list(0)
SEL_TRUE_gamsel_ov=list(0)
SEL_FALSE_gamsel_ov=list(0)
SEL_TRUE_res_ov=list(0)
SEL_FALSE_res_ov=list(0)
SEL_TRUE_init_ov_both=list(0)
SEL_FALSE_init_ov_both=list(0)

NLIN_TP_res=list(0)
NLIN_FP_res=list(0)
NLIN_TP_init=list(0)
NLIN_FP_init=list(0)
NLIN_TP_gamsel=list(0)
NLIN_FP_gamsel=list(0)

COV_MAD=list(0)
MAD_OV=list(0)

for(m in 1:50){
  n=250
  q=20
  p=50
  #X=cbind(runif(n,-1,1),runif(n,-1,1),runif(n,-1,1))
  X=matrix(runif(n*p,-1,1),nrow=n,ncol=p)
  #E=mvrnorm(n,rep(0,q),diag(q))
  #for(oe in 1:ncol(X)){
  #  X[,oe]=runif(1,0.5,10)*X[,oe]
  #}
  
  rho=0.9
  bd <- rho^(abs(outer(1:q, 1:q, "-")))
  
  E=mvrnorm(n,rep(0,q),bd)
  Y=matrix(0,nrow=n,ncol=q)
  
  
  #Y[,1]=5*X[,1]^(3) + 3*X[,2]^(2)
  #Y[,1]=5*sin(5*X[,1])
  #Y[,1]=1*(X[,1]^(2))
  #Y[,2]=X[,2]*1
  
  #SELECTION OF NON-LINEAR RELATIONSHIPS - Define functions 
  
  f1=function(x){0.5*(1-exp(-2*x))}
  f2=function(x){0.5*x^(2)}
  f3=function(x){0.5*x^(3)}
  f4=function(x,sigma){0.5*((1/(sqrt(2*pi)*sigma))*exp(-(x^(2)/(2*sigma^(2)))))}
  #Will Use 0.1, 0.01 as sigma peaks
  #f5=function(x){10*(x^(4)-2*x^(2))}
  f5=function(x){0.5*x}
  
  
  
  #SELECTION OF NON-LINEAR FUNCTION POSTIIONS (pxq) - levels of sparsity and positions 
  #response=seq(1,q/2,1)
  response=seq(1,5,1)
  covs=seq(1,p,1)
  
  #Two levels - 
  #First layer... select which responses not sparse
  spar=4
  response=sample(response,spar)
  response=sort(response)
  #Second layer of sparsity, select # of covariates for each chosen response
  predictors=sample(c(1:5),spar,replace=TRUE)
  #predictors=rep(10,spar)
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
      
      #if(inte[j]==6){
      #Y[,response[i]]=Y[,response[i]]+f6(X[,ct[j]])
      #}
      
    }
  }
  
  Y_truth=Y
  Y=Y + E
  
  #--------------------------------------------
  #STEP 1: COMPUTE BIG KEY 
  key=matrix(0,nrow=p,ncol=q)
  for(i in 1:length(response)){
    key[cov4response[[i]],response[i]]=1
  }
  
  for(i in 1:length(funccov4response)){
    indi=which(funccov4response[[i]]==5)
    econd=cov4response[[i]][indi]
    key[econd,response[i]]=2
  }
  #------------------------------------------
  lambda=matrix(0,nrow=p,ncol=q)
  for(i in 1:q){
    for(j in 1:p){
      lambda[j,i]=smooth.spline(X[,j],Y[,i])$lambda
    }
  }
  #ptm=proc.time()
  
  #X=X[,1:2]
  #GET BASES FOR EACH SIMULATION USING THIS CODE 
  rot=rep(list(rep(list(0))),q)
  U_tilde=list(0)
  U_tildee=list(0)
  U_lme=list(0)
  Omega_tilde_d=list(0)
  vo=list(0)
  for(j in 1:q){
    for(i in 1:p){
      valz=quantile(X[,i],probs=seq(0,1,0.1))
      vo[[i]]=round(valz,1)
      #Creation of Spline Basis 
      knots <- c(valz[1], valz[1], valz[1], valz, valz[11], valz[11], valz[11])
      M=length(valz)-2
      inner_knots <- valz
      Phi <- splineDesign(c(rep(range(inner_knots), 3), inner_knots), X[,i])
      Omega <- pen_mat(inner_knots)
      Phi_svd <- svd(Phi)
      Omega_tilde <- t(crossprod(Phi_svd$v, Omega %*% Phi_svd$v)) / Phi_svd$d
      Omega_tilde <- t(Omega_tilde) / Phi_svd$d
      Omega_tilde_svd <- svd(Omega_tilde) 
      Omega_tilde_d[[i]]=diag(Omega_tilde_svd$d)
      root=diag(Omega_tilde_svd$d)
      root=root[1:(ncol(root)-2),1:(ncol(root)-2)]
      U_tildee[[i]] <- Phi_svd$u %*% Omega_tilde_svd$u
      U_lme[[i]]= U_tildee[[i]][,1:11]%*%sqrt(solve(root))
      Mj=(1/n)*(t(U_tildee[[i]])%*%U_tildee[[i]]) + lambda[i,j]*Omega_tilde_d[[i]]
      hed=chol(Mj)
      rot[[j]][[i]]=U_tildee[[i]]%*%solve(hed)
      rot[[j]][[i]]=rot[[j]][[i]][,1:11]
      print(i)
    }
    print(j)
  }
  
  #-------------------------------------------
  #STEP 2: ON OTHER PAGE - Set up BASES 
  #-------------------------------------------
  B=matrix(0,nrow=p,ncol=q)
  library(doSNOW)
  NumberOfCluster <- 10
  cl <- makeCluster(NumberOfCluster)
  registerDoSNOW(cl)
  x1=foreach(i=1:q)%dopar%{
    library('glmnet')
    hey=cv.glmnet(X,Y[,i],intercept=FALSE)
    return(coef(hey,s=hey$lambda.1se)[-1,])
  }
  stopCluster(cl)
  
  B=matrix(unlist(x1),nrow=p,ncol=q)
  por=which(B!=0,arr.ind=TRUE)
  if(length(por)!=0){
    for(i in 1:nrow(por)){
      or=lm(Y[,por[i,2]]~X[,por[i,1]]-1)
      B[por[i,1],por[i,2]]=coef(or)
    }}
  

  
  resd=Y-X%*%B
  re=rot
  #re=cbind(rot[[1]][,1:11],rot[[2]][,1:11],rot[[3]][,1:11])
  comat=matrix(0,nrow=11*p,ncol=q)
  
  grp=rep(1:p,each=11) 
  library(doSNOW)
  NumberOfCluster <- 10
  cl <- makeCluster(NumberOfCluster)
  registerDoSNOW(cl)
  x2=foreach(i=1:q)%dopar%{
    library('gglasso')
    library('rlist')
    eds=cv.gglasso(x=list.cbind(re[[i]]),y=resd[,i],group=rep(1:p,each=11),intercept=FALSE,lambda.factor=0.75)
    return(coef(eds,s=eds$lambda.1se)[-1,])
  }
  stopCluster(cl)
  
  
  comat=matrix(unlist(x2),nrow=11*p,ncol=q)
  indicesz=matrix(0,nrow=p,ncol=q)
  for(i in 1:q){
    derd=comat[,i]
    dr=which(derd!=0,arr.ind=TRUE)
    indix=unique(grp[dr])
    indicesz[indix,i]=1
    print(i)
  }
  
  pod=which(indicesz==1,arr.ind=TRUE)
  pad=sort(unique(pod[,2]))
  nlind=matrix(0,nrow=n,ncol=q)
  
  #RIGHT HERE.... YOU ARE CURRENTLY SUBTRACTING OFF RESIDUAL AND 
  #FITTING A GAM TO THE RESULTS.... WHAT YOU WANT TO DO
  #IS TAKE OUT RESIDUAL FROM Y, AND RE-FIT FOR NEW SELECTED COMPONENTS
  linind=matrix(0,nrow=p,ncol=q)
  linind[which(B!=0,arr.ind=TRUE)]=1
  
  #Create linear terms
  #por is the linear one.... need to 
  rosy=sort(unique(c(por[,2],pad)))
  repido=c(rep(1,nrow(U_lme[[1]])))
  fixpred=matrix(0,nrow=n,ncol=q)
  
  if(length(rosy)!=0){
  for(i in 1:length(rosy)){
    cv=which(linind[,rosy[i]]!=0)
    idic=which(pod[,2]==rosy[i])
    nlinpos=pod[idic,1]
    cv=unique(c(cv,nlinpos))
    if(length(cv)!=0){
      e=""
      for(j in 1:length(cv)){
        e=paste(e,"+","X[,",cv[j],"]",sep="")
      }
      e=substring(e,2)}
    if(length(idic)==1){
      UL=as.matrix(U_lme[[pod[idic,1]]])
      mod_gam2=lme(as.formula(paste('Y[,',rosy[i],']','~',e,"+1",sep='')),random=list(repido=pdIdent(~UL-1)))
      fixpred[,rosy[i]]=fitted(mod_gam2,level=0)
    }
    if(length(idic)>1){
      UL=list(0)
      for(k in 1:length(idic)){
        UL[[k]]=as.matrix(U_lme[[pod[idic[k],1]]])
      }
      rusy=list(0)
      for(l in 1:length(UL)) {# assign function within loop
        assign(paste0('UNC_',l),as.matrix(UL[[l]]))
        #assign(paste0("variable_", i), pdIdent(as.formula(paste('~UNC_',i,'-1',sep=''))))
        rusy[[l]]=pdIdent(as.formula(paste('~UNC_',l,'-1',sep='')))
      }
      
      mod_gam2=lme(as.formula(paste('Y[,',rosy[i],']','~',e,'+1',sep='')),
                   random=list(repido=pdBlocked(rusy)))
      
      fixpred[,rosy[i]]=fitted(mod_gam2,level=0)
    }
    if(length(idic)!=0){
      nlind[,rosy[i]]=fitted(mod_gam2)
    }
    if(length(idic)==0){
      nlind[,rosy[i]]=X%*%B[,rosy[i]]
      fixpred[,rosy[i]]=X%*%B[,rosy[i]]
    }
  }
  }
  
  resd=Y-nlind
  
  est=cov(resd)
  d=CVglasso(X=resd,diagonal=TRUE,path=TRUE)
  prec=d$Omega
  resd=Y-nlind
  
  B_init=B
  fixpred_init=fixpred
  pod_init=pod
  nonlin_init=nlind
  prec_init=prec
  
  #--------------------------------------------
  #Run Gamsel for Comparisons at the End 
  #--------------------------------------------
  Gam_B=matrix(0,nrow=p,ncol=q)
  Gam_nlin=matrix(0,nrow=p,ncol=q)
  for(i in 1:ncol(Y)){
    ed=cv.gamsel(X,Y[,i],degrees=11)
    Gam_B[,i]=ed$gamsel.fit$alphas[,ed$index.1se]
    dr=which(ed$gamsel.fit$betas[,ed$index.1se]!=0,arr.ind=TRUE)
    grp=rep(1:p,each=11) 
    indix=unique(grp[dr])
    Gam_nlin[indix,i]=1
    print(i)
  }
  por=which(Gam_B!=0,arr.ind=TRUE)
  Gam_B[por]=1
  #-------------------------------------------
  #STEP 3: Run Selection Procedure 
  #-------------------------------------------
  mad=c()
  mad[1]=mean(abs(resd))
  
  Bt=list(0)
  #nt=list(0)
  pc=list(0)
  pods=list(0)
  nlinds=list(0)
  
  for(l in 2:5){
    
    Ynew=matrix(0,nrow=n,ncol=q)
    
    for(i in 1:ncol(Y)){
      h=prec[i,i]
      rest=prec[i,-i]
      alpha=-1*(rest/h) 
      Ystar=Y[,-i]
      
      ##Set the Y
      Ynew[,i]=Y[,i]-(nlind[,i]-fixpred[,i])-Ystar%*%alpha+(nlind[,-i]%*%alpha)
      #print(i)
    }
    
    #Run Lasso (Using 1se)
    #B=matrix(0,nrow=p,ncol=q)
    library(doSNOW)
    NumberOfCluster <- 10
    cl <- makeCluster(NumberOfCluster)
    registerDoSNOW(cl)
    x3=foreach(i=1:q)%dopar%{
      library('glmnet')
      hey=cv.glmnet(X,Ynew[,i],intercept=FALSE)
      return(coef(hey,s=hey$lambda.1se)[-1,])
    }
    stopCluster(cl)
    B=matrix(unlist(x3),nrow=p,ncol=q)
    
    por=which(B!=0,arr.ind=TRUE)
    if(length(por)!=0){
      for(i in 1:nrow(por)){
        or=lm(Ynew[,por[i,2]]~X[,por[i,1]]-1)
        B[por[i,1],por[i,2]]=coef(or)
      }
    }
    
    Bt[[l]]=B
    
    #rename Ynew to subtract out the linear term 
    for(i in 1:ncol(Y)){
      h=prec[i,i]
      rest=prec[i,-i]
      alpha=-1*(rest/h) 
      Ystar=Y[,-i]
      Ynew[,i]=Y[,i]-X%*%B[,i]-Ystar%*%alpha+(nlind[,-i]%*%alpha)
      print(i)
    }
    
    re=rot
    comat=matrix(0,nrow=11*p,ncol=q)
    grp=rep(1:p,each=11) 
    
    library(doSNOW)
    NumberOfCluster <- 10
    cl <- makeCluster(NumberOfCluster)
    registerDoSNOW(cl)
    x4=foreach(i=1:q)%dopar%{
      library('gglasso')
      library('rlist')
      eds=cv.gglasso(x=list.cbind(re[[i]]),y=Ynew[,i],group=rep(1:p,each=11),intercept=FALSE)
      return(coef(eds,s=eds$lambda.1se)[-1,])
    }
    stopCluster(cl)
    
    comat=matrix(unlist(x4),nrow=11*p,ncol=q)
    
    indicesz=matrix(0,nrow=p,ncol=q)
    for(i in 1:q){
      derd=comat[,i]
      dr=which(derd!=0,arr.ind=TRUE)
      indix=unique(grp[dr])
      indicesz[indix,i]=1
      print(i)
    }
    
    pod=which(indicesz==1,arr.ind=TRUE)
    pods[[l]]=pod
    pad=sort(unique(pod[,2]))
    #nlind=matrix(0,nrow=n,ncol=q)
    
    
    #RIGHT HERE IS WHERE THE CORRECTION IS NEEDED... 
    #REWRITE YNEW TO NOT HAVE LINEAR OR NON-LINEAR TERM SUBTRACTED OUT 
    #Need to rewrite this then 
    for(i in 1:ncol(Y)){
      h=prec[i,i]
      rest=prec[i,-i]
      alpha=-1*(rest/h) 
      Ystar=Y[,-i]
      Ynew[,i]=Y[,i]-Ystar%*%alpha+(nlind[,-i]%*%alpha)
      print(i)
    }
    
    linind=matrix(0,nrow=p,ncol=q)
    linind[which(B!=0,arr.ind=TRUE)]=1
    
    #Create linear terms
    #por is the linear one.... need to 
    rosy=unique(c(por[,2],pad))
    repido=c(rep(1,nrow(U_lme[[1]])))
    fixpred=matrix(0,nrow=n,ncol=q)
    nlind=matrix(0,nrow=n,ncol=q)
    
    if(length(rosy)!=0){
    for(i in 1:length(rosy)){
      cv=which(linind[,rosy[i]]!=0)
      idic=which(pod[,2]==rosy[i])
      nlinpos=pod[idic,1]
      cv=unique(c(cv,nlinpos))
      if(length(cv)!=0){
        e=""
        for(j in 1:length(cv)){
          e=paste(e,"+","X[,",cv[j],"]",sep="")
        }
        e=substring(e,2)}
      if(length(idic)==1){
        UL=as.matrix(U_lme[[pod[idic,1]]])
        mod_gam2=lme(as.formula(paste('Ynew[,',rosy[i],']','~',e,"+1",sep='')),random=list(repido=pdIdent(~UL-1)))
        fixpred[,rosy[i]]=fitted(mod_gam2,level=0)
      }
      if(length(idic)>1){
        UL=list(0)
        for(k in 1:length(idic)){
          UL[[k]]=as.matrix(U_lme[[pod[idic[k],1]]])
        }
        rusy=list(0)
        for(o in 1:length(UL)) {# assign function within loop
          assign(paste0('UNC_',o),as.matrix(UL[[o]]))
          #assign(paste0("variable_", i), pdIdent(as.formula(paste('~UNC_',i,'-1',sep=''))))
          rusy[[o]]=pdIdent(as.formula(paste('~UNC_',o,'-1',sep='')))
        }
        #dq=pdBlocked(rusy,nam=namu
        mod_gam2=lme(as.formula(paste('Ynew[,',rosy[i],']','~',e,'+1',sep='')),
                     random=list(repido=pdBlocked(rusy)))
        fixpred[,rosy[i]]=fitted(mod_gam2,level=0)
      }
      if(length(idic)!=0){
        nlind[,rosy[i]]=fitted(mod_gam2)
      }
      if(length(idic)==0){
        nlind[,rosy[i]]=X%*%B[,rosy[i]]
        fixpred[,rosy[i]]=X%*%B[,rosy[i]]
      }
    }}
    #----END HERE 
    #for(i in 1:length(pad)){
    #tr=which(indicesz[,pad[i]]!=0)
    #if(length(tr)!=0){
    #  e=""
    #  for(j in 1:length(tr)){
    #    e=paste(e,"+","s(X[,",tr[j],"])",sep="")
    #  }
    #  e=substring(e,2)
    #  f=paste("Ynew[,",pad[i],"]","~",sep="")
    #  edline=paste(f,e)
    #  mod_gam2=gam(as.formula(edline))
    #}
    #if(length(pad!=0)){
    #nlind[,pad[i]]=mod_gam2$fitted.values
    #}
    #}
    nlinds[[l]]=nlind
    
    #resd=resd-nlind
    #nt[l]=sqrt(sum((X%*%B[,1]+nlind[,1]-5*(X[,1]^2))^(2))/n)
    
    resd=Y-nlind
    
    est=cov(resd)
    d=CVglasso(X=resd,diagonal=TRUE)
    prec=d$Omega
    
    pc[l]=mean(abs((bd-solve(prec))[upper.tri(bd-solve(prec))]))
    mad[l]=mean(abs(resd))
    print(l)
  }
  Bt[[1]]=B_init
  
  #Create True Parameter Matrices 
  #We have responses,chosen covariates, and chosen functions
  #Of markers, can use Bt and pod to demarcate these 
  
  
  #Non-Zero Positions
  
  #keyscatterplots
  
  #lab=which(key!=0,arr.ind=TRUE)
  #par(mfrow=c(4,2))
  #for(i in 1:nrow(lab)){
  #  plot(X[,lab[i,1]],Y[,lab[i,2]],main=paste0("Cov"," ",lab[i,1]," ","vs",
  #                                             " ","Y-",lab[i,2]))
  #}
  
  
  #LINEAR TRUTH
  marks=which(key==2,arr.ind=TRUE)
  B_TRUTH=matrix(0,nrow=p,ncol=q)
  B_TRUTH[marks]=0.5
  B_first=mean(abs(B_init[marks]-B_TRUTH[marks]))
  B_res=mean(abs(Bt[[5]][marks]-B_TRUTH[marks]))
  
  #B_TRUTH RATIO
  Btrue[[m]]=B_res/B_first
  
  #B_OVERALL RATIO 
  B_ov=abs(Bt[[5]]-B_TRUTH)
  B_ov[which(key==1,arr.ind=TRUE)]=0
  
  B_ov_init=abs(B_init-B_TRUTH)
  B_ov_init[which(key==1,arr.ind=TRUE)]=0
  Bov[[m]]=mean(B_ov)/mean(B_ov_init)
  
  #Y_TRUTH mean absolute deviation 
  
  one=mean(abs(Y_truth - (nlinds[[5]])))
  two=mean(abs(Y_truth - (nonlin_init)))
  Ytrue[[m]]=one/two

  #TRUE POSITIVE (B)
  yar=Bt[[5]][which(key==2,arr.ind=TRUE)]
  sumr=rep(0,length(yar))
  for(i in 1:length(yar)){
    if(yar[i]!=0){
      sumr[i]=1
    }
  }
  B_TRUE_res[[m]]=mean(sumr)
  
  
  yar=B_init[which(key==2,arr.ind=TRUE)]
  sum0=rep(0,length(yar))
  for(i in 1:length(yar)){
    if(length(yar)!=0){
      if(yar[i]!=0){
        sum0[i]=1
      }
    }
  }
  
  B_TRUE_init[[m]]=mean(sum0)
  
  yar=Gam_B[which(key==2,arr.ind=TRUE)]
  sumr=rep(0,length(yar))
  for(i in 1:length(yar)){
    if(yar[i]!=0){
      sumr[i]=1
    }
  }
  
  B_TRUE_gamsel[[m]]=mean(sumr)
  
  
  yar=Gam_B[which(key==0,arr.ind=TRUE)]
  sumr=rep(0,length(yar))
  for(i in 1:length(yar)){
    if(yar[i]==0){
      sumr[i]=1
    }
  }
  B_FALSE_gamsel[[m]]=1-mean(sumr)
  
  
  yar=Bt[[5]][which(key==0,arr.ind=TRUE)]
  sumr=rep(0,length(yar))
  for(i in 1:length(yar)){
    if(yar[i]==0){
      sumr[i]=1
    }
  }
  B_FALSE_res[[m]]=1-mean(sumr)
  
  yar=B_init[which(key==0,arr.ind=TRUE)]
  sum0=rep(0,length(yar))
  for(i in 1:length(yar)){
    if(yar[i]==0){
      sum0[i]=1
    }
  }
  B_FALSE_init[[m]]=1-mean(sum0)
  
  
  yar=B_init[which(key!=0,arr.ind=TRUE)]
  sumr=rep(0,length(yar))
  #yanxi=nonlin_init[which(key!=0),arr.ind=TRUE]
  #sumo=rep(0,length(yanxi))
  for(i in 1:length(yar)){
    if(yar[i]!=0){
      sumr[i]=1
    }
  }
  SEL_TRUE_init_ov[[m]]=mean(sumr)
  
  yar=B_init[which(key==0,arr.ind=TRUE)]
  sumr=rep(0,length(yar))
  for(i in 1:length(yar)){
    if(yar[i]==0){
      sumr[i]=1
    }
  }
  SEL_FALSE_init_ov[[m]]=1-mean(sumr)
  
  
  #NON-LINEAR TP
  key_np=matrix(0,nrow=p,ncol=q)
  key_np[which(key==1,arr.ind=TRUE)]=1
  
  nlin_resy=matrix(0,nrow=p,ncol=q)
  nlin_resy[pods[[5]]]=1
  
  rav=key_np-nlin_resy
  NLIN_TP_res[[m]]= 1-mean(rav[which(key==1,arr.ind=TRUE)])
  
  
  nlin_resy=matrix(0,nrow=p,ncol=q)
  nlin_resy[pod_init]=1
  
  rav=key_np-nlin_resy
  NLIN_TP_init[[m]]=1-mean(rav[which(key==1,arr.ind=TRUE)])
  
  rav=key_np-Gam_nlin
  NLIN_TP_gamsel[[m]]=1-mean(rav[which(key==1,arr.ind=TRUE)])
  
  #mean((key_np-nlin_resy))

  #FALSE POSITIVE 
  nlin_resy=matrix(0,nrow=p,ncol=q)
  nlin_resy[pods[[5]]]=1
  
  rav=key_np-nlin_resy
  NLIN_FP_res[[m]]= mean(abs(rav[which(key!=1,arr.ind=TRUE)]))
  
  nlin_resy=matrix(0,nrow=p,ncol=q)
  nlin_resy[pod_init]=1
  
  rav=key_np-nlin_resy
  NLIN_FP_init[[m]]=mean(abs(rav[which(key!=1,arr.ind=TRUE)]))
  
  
  rav=key_np-Gam_nlin
  NLIN_FP_gamsel[[m]]=mean(abs(rav[which(key!=1,arr.ind=TRUE)]))
  
  
  COV_MAD[[m]] = pc[[5]]/mean(abs((bd-solve(prec_init))[upper.tri(bd-solve(prec_init))]))
  MAD_OV[[m]] = mad[[5]]/mad[[1]]
  
  
  #OVERALL TP - Co-MAdRE
  yar=Bt[[5]][which(key!=0,arr.ind=TRUE)]
  yar[which(yar!=0,arr.ind=TRUE)]=1
  nlin_resy=matrix(0,nrow=p,ncol=q)
  nlin_resy[pods[[5]]]=1
  summy=yar+nlin_resy[which(key!=0,arr.ind=TRUE)]
  summy[which(summy>1)]=1

  SEL_TRUE_res_ov[[m]]=mean(summy)
  
  #Overall FP - Co-MAdRE
  yar=Bt[[5]][which(key==0,arr.ind=TRUE)]
  yar[which(yar!=0,arr.ind=TRUE)]=1
  nlin_resy=matrix(0,nrow=p,ncol=q)
  nlin_resy[pods[[5]]]=1
  summy=yar+nlin_resy[which(key==0,arr.ind=TRUE)]
  summy[which(summy>1)]=1
  SEL_FALSE_res_ov[[m]]=mean(summy)
  
  #OVERALL TP - Gamsel 
  yar=Gam_B[which(key!=0,arr.ind=TRUE)]
  nlin_resy=Gam_nlin
  summy=yar+nlin_resy[which(key!=0,arr.ind=TRUE)]
  summy[which(summy>1)]=1
  SEL_TRUE_gamsel_ov[[m]]=mean(summy)

  #OVERALL FP - Gamsel 
  yar=Gam_B[which(key==0,arr.ind=TRUE)]
  nlin_resy=Gam_nlin
  summy=yar+nlin_resy[which(key==0,arr.ind=TRUE)]
  summy[which(summy>1)]=1
  SEL_FALSE_gamsel_ov[[m]]=mean(summy)
  
  
  #OVERALL TP - Independent
  yar=B_init[which(key!=0,arr.ind=TRUE)]
  yar[which(yar!=0,arr.ind=TRUE)]=1
  nlin_resy=matrix(0,nrow=p,ncol=q)
  nlin_resy[pod_init]=1
  summy=yar+nlin_resy[which(key!=0,arr.ind=TRUE)]
  summy[which(summy>1)]=1
  SEL_TRUE_init_ov_both[[m]]=mean(summy)

  #OVERALL FP - Independent  
  yar=B_init[which(key==0,arr.ind=TRUE)]
  yar[which(yar!=0,arr.ind=TRUE)]=1
  nlin_resy=matrix(0,nrow=p,ncol=q)
  nlin_resy[pod_init]=1
  summy=yar+nlin_resy[which(key==0,arr.ind=TRUE)]
  summy[which(summy>1)]=1
  SEL_FALSE_init_ov_both[[m]]=mean(summy)
  
  print(m)}








