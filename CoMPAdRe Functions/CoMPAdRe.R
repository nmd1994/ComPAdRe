#CoMPAdRe function 
#Uses helper function pen_mat and relies on the dependencies shown in the worked tutorial
# X is n x p matrix of p predictors, Y is n x q matrix of q responses, 
#lambda.fac is penalization term for cross validation tuning parameter selection for non-linear predictor-response associations 
#niter is number of iterations desired for algorithm 



Compadre=function(X, Y, lambda.fac, niter){
  n=nrow(Y)
  q=ncol(Y)
  p=ncol(X)
  
  lambda=matrix(0,nrow=p,ncol=q)
  for(i in 1:q){
    for(j in 1:p){
      lambda[j,i]=smooth.spline(X[,j],Y[,i])$lambda
    }
  }
  
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
      #print(i)
    }
    #print(j)
  }
  
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
  comat=matrix(0,nrow=11*p,ncol=q)
  
  grp=rep(1:p,each=11) 
  library(doSNOW)
  NumberOfCluster <- 10
  cl <- makeCluster(NumberOfCluster)
  registerDoSNOW(cl)
  x2=foreach(i=1:q)%dopar%{
    library('gglasso')
    library('rlist')
    eds=cv.gglasso(x=list.cbind(re[[i]]),y=resd[,i],group=rep(1:p,each=11),intercept=FALSE,lambda.factor=lambda.fac)
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
  
  Bt=list(0)
  pods=list(0)
  nlinds=list(0)
  
  for(l in 2:niter){
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
      #print(i)
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
      eds=cv.gglasso(x=list.cbind(re[[i]]),y=Ynew[,i],group=rep(1:p,each=11),intercept=FALSE,
                     lambda.factor=lambda.fac)
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
      #print(i)
    }
    
    pod=which(indicesz==1,arr.ind=TRUE)
    pods[[l]]=pod
    pad=sort(unique(pod[,2]))
    
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
            rusy[[o]]=pdIdent(as.formula(paste('~UNC_',o,'-1',sep='')))
          }
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
    nlinds[[l]]=nlind
    
    resd=Y-nlind
    
    est=cov(resd)
    d=CVglasso(X=resd,diagonal=TRUE)
    prec=d$Omega
    #paste0("iteration", " ", l)
  }
  values=list(Bt[[niter]],nlinds[[niter]],pods[[niter]],solve(prec))
  names(values)=c("linear coefficients","estimated functions f-hat","non-linear selection","estimated covariance")
  return(values)  
}
