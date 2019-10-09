varsel1 <- function(X,Y,tt,order){
  
  nobs = nrow(X)
  ndim = ncol(X)
  dfo = rep(order-1,ndim)
  index = rep(1:ndim,times=dfo)
  aa2 <- c()
  for(k in 1:ndim){
    aa1 <- c()
    for(j in 1:(order-1)){
      aa <- c()
      for(i in 1:length(tt)){
        aa <- c(aa,Legendre.modelve(tt[i],np.order = j,tmin = min(tt), tmax = max(tt)))
      }
      aa1 <- cbind(aa1,aa*X[,k])
    }
    aa2 <- cbind(aa2,aa1)
  }
 

  Xc = scale(aa2,center=T,scale=T)
  n = nrow(Xc)
  
  connect = matrix(0,nrow=ndim,ncol=ndim)
  coefest = matrix(0,nrow=sum(dfo),ncol=ndim)
  regfun = vector("list",length=ndim)
  for(i in 1:ndim)
  {
    yc <- Y[,i]-mean(Y[,i])
    
    out1 <- GrpLasso(X=Xc,y=yc,index=index,lambda=30,crit="BIC")
    var.grp <- out1$var.select  # genes selected
    coef.grp <- out1$coef
    
    ### Adaptive Group Lasso
    index.adp <- index[is.element(index,var.grp)]
    W.adp = sapply(1:length(var.grp),function(j) sqrt(sum(coef.grp[index.adp==var.grp[j]]^2)))
    Xc.adp = Xc[,is.element(index,var.grp)]
    Xcs.adp = scale(Xc.adp,center=F,scale=rep(1/W.adp,times=dfo[var.grp]))
    init.adp = coef.grp/rep(W.adp,times=dfo[var.grp])
    lambda = lambdamax(Xcs.adp,yc,index=index.adp,coef.ini=init.adp,
                       penscale=sqrt,center=F,standardize=F,model=LinReg())*0.7^(1:18)
    out2 = GrpLasso(X=Xc.adp,y=yc,W=W.adp,index=index.adp,ini=coef.grp,
                    lambda=lambda,crit="BIC")
    var.adp = out2$var.select
    coef.adp = out2$coef
    connect[i,var.adp] <-  1
    coefest[is.element(index,var.adp),i] <-  coef.adp
    regfun[[i]] <-  sapply(var.adp,function(ii) Xc[,is.element(index,ii)]%*%coefest[is.element(index,ii),i])
    connect[i,var.adp] = 1
    coefest[is.element(index,var.adp),i] = coef.adp
    regfun[[i]] = sapply(var.adp,function(ii) Xc[,is.element(index,ii)]%*%coefest[is.element(index,ii),i])
    cat("var=",i,var.adp,"\n")
  }
  return(list(connect=connect,regfun=regfun,coefest=coefest))
}