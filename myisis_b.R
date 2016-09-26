# Iterative Sure Independence Screening with Elastic Net penalty.
# This code is modified based on R package "SIS".
# 
# y: response variable ~ binomial(m,p)
# x: design matrix (n x p matrix)
# 
# Input
###### Need a format of two columns for input for y: cbind(m-y,y) ###### 

# a simple example
##### Generating data:
# n <- 300
# p <- 1000
# m <- 5
# trueb <- rep(0,p)
# trueb[1:4] <- c(3,3,3,3) # coefficients
# x <- matrix(rnorm(n*p), nrow=n, ncol=p) # design matrix
# pr <- 1/(1+exp(-x%*%trueb)) # probability for binomial distribution
# y <- rbinom(n, size=m, prob=pr)
#
##### Iterative Sure Independent Screening with elastic net:
# res <- myisis_b(x,cbind(m-y,y), alpha=0.5)
# ix <- res$ix
# print(ix)
#
##### Applying to bootstrap samples:
# library(caret)
# library(plyr)
# train <- createResample(y,times=10,list=TRUE)
# coef <- list("vector")
# set.seed(1)
# for ( i in 1:10 ){
# coef[[i]] <- myisis_b(x[train[[i]],], cbind(m-y[train[[i]]], y[train[[i]]]), alpha=0.5, nsis=10)$ix
# }
# cc <- c(unlist(coef))
# idx <- sort(count(cc)$freq,decreasing=T,index.return=T)$ix
#
##### Printing final results:
# print(count(cc)[idx[1:10],])





library(glmnet)
# alpha=1, lasso
# alpha=0, ridge

standardize <- function(X){
  center = colMeans(X)
  X.c = sweep(X, 2, center)
  unit.var = sqrt(apply(X.c, 2, crossprod))
  val = sweep(X.c, 2, unit.var, "/")
  return(val)
}

mg <- function(candind, x=x, y=y, ix1){
  ones <- rep(1,dim(x)[1])
  margfit <- coef( glm.fit(cbind(ones, x[,candind], x[,ix1]), y, family=binomial() ) )[2]
  return(margfit)
}

myisis_b <- function(x, y, nfolds = 10, alpha = 0.5, nsis=NULL, standardize=TRUE){

  n <- dim(x)[1];	p <- dim(x)[2]
  
  models = vector("list")
  if(is.null(nsis)==TRUE){ 
    nsis = floor(n/(4*log(n)))
    if(p < n){ nsis = floor(p/3) }
  }
  
  if(standardize == TRUE){ 
    old.x = x
    x = standardize(x)
  }
  
  iterind = 0
  repeat{
    if (iterind==0){
      margc <- abs(cor(x,y[,2]))
      rankc <- sort(margc,decreasing=T,index.return=T)
      d <- floor(2/3*nsis)
      ix0 <- rankc$ix[1:d]
      ix0 <- sort(ix0)
    }
    cat("Iter ",iterind,": ", ix0, "\n")
    
    iterind <- iterind+1
    cv.fit <- cv.glmnet(x[,ix0],y,nfolds=nfolds,alpha=alpha,family="binomial")
      beta <- coef(cv.fit,s="lambda.min")[-1]
      ix1 <- ix0[which(abs(beta)>1e-10)]

    cat("Iter ",iterind,": ix1: ", ix1, "\n")
    
        if(length(ix1) >= nsis){
          cat("Maximum number of variables selected \n")        
          break
        }
    
    models[[iterind]] = ix1; flag.models = 0
    if(iterind > 1){
      for(j in 1:(iterind-1)){
        if(identical(models[[j]],ix1) == TRUE) flag.models = 1
      }
    }
    if(flag.models==1){
      ix0 = ix1
      cat("Model already selected \n")
      break   
    }
    
    candind = setdiff(1:p, ix1)
    
    margc <- abs(sapply(candind,mg,x,y,ix1))
    rankc <- sort(margc,decreasing=T,index.return=T)
    newix <- rankc$ix[1:(nsis-length(ix1))]
    newix <- sort(candind[newix])
    cat("Iter ",iterind,": newix: ", newix, "\n")
    
    ix1 <- sort(c(ix1,newix))
    
    if(setequal(ix1,ix0)){
      ix0 = sort(setdiff(ix1,newix))
      cat("Model already selected \n")
      break
    }

    ix0 <- ix1 
  }
  
  return(list(ix=ix0))
  
}

