varsel <- function(X,Y,tt){
  
  nobs = nrow(X)
  ndim = ncol(X)
  dfo = rep(ceiling(nobs^(1/3)+4),ndim)
  index = rep(1:ndim,times=dfo)
  BS = vector("list",length=ndim)
  for(i in 1:ndim){
    knot = quantile(unique(X[,i]),probs=seq(0,1,length=dfo[i]-2))
    allknot = c(rep(knot[1],4),knot[2:(dfo[i]-3)],rep(knot[dfo[i]-2],4))
    BS[[i]] = SplineBasis(allknot)
    if(i == 1){
      XX = evaluate(BS[[i]],X[,i])
    }else{
      XX = cbind(XX,evaluate(BS[[i]],X[,i]))
    }
  }
  Xc = scale(XX,center=T,scale=F)
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


GrpLasso <- function(X,y,W=NULL,index,ini=rep(0,ncol(X)),lambda=NULL,
                     crit=c("BIC","EBIC"),center=F)
{
  if(center==T){
    y = y-mean(y)
    X = scale(X,center=T,scale=F)
  }
  n = nrow(X)
  ind = unique(index)
  p = length(ind)
  dfo = sapply(1:p,function(j) sum(index==ind[j]))
  
  # fit model for a sequence of penalty parameters
  if(!is.null(W)){
    X = scale(X,center=F,scale=rep(1/W,times=dfo))
    ini = ini/rep(W,times=dfo)
  }
  
  # set up the candidates for penalty parameter
  if(is.null(lambda)){
    lambda = lambdamax(X,y,index=index,coef.ini=ini,penscale=sqrt,center=F,
                       standardize=F,model=LinReg())*0.9^(1:20)
  }else if(length(lambda)==1){
    lambda = lambdamax(X,y,index=index,coef.ini=ini,penscale=sqrt,center=F,
                       standardize=F,model=LinReg())*0.9^(1:lambda)
  }
  
  fit = grplasso(X,y,index=index,lambda=lambda,model=LinReg(),center=F,
                 standardize=F,coef.ini=ini,penscale=sqrt,
                 control=grpl.control(update.hess="lambda",tol=10^-8,trace=0))
  
  # calculate BIC/EBIC
  nlambda = length(lambda)
  rss = sapply(1:nlambda,function(j) sum((y-fit$fitted[,j])^2))
  var.select = sapply(1:nlambda,function(j) unique(index[abs(fit$coef[,j])>0]))
  dfo.lambda = sapply(1:nlambda,function(j) sum(dfo[is.element(ind,var.select[[j]])]))
  if(crit!="BIC" & crit!="EBIC"){
    cat("Error: Criterion not implemented. Reset to BIC!\n")
    crit = "BIC"
  }
  if(crit=="BIC"){
    bic = log(rss)+dfo.lambda*log(n)/n
  }else if(crit == "EBIC"){
    bic = log(rss)+dfo.lambda*log(n)/n+0.5*dfo.lambda*log(p)/n
  }
  
  # select model with smallest value of selection criterion
  id.ss = which.min(bic)
  var.ss = var.select[[id.ss]]
  fit.ss = fit$fitted[,id.ss]
  coef.ss = fit$coef[,id.ss]
  if(!is.null(W)){
    coef.ss = coef.ss*rep(W,times=dfo)
  }
  coef.ss = coef.ss[is.element(index,var.ss)]
  
  return(list(var.select=var.ss,coefficients=coef.ss,fitted=fit.ss,BIC=bic,
              lambda=lambda,fit.org=fit))
}

Legendre <-function( t, np.order, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- np.order-1
  L <- 1 
  if (np.order>=1)
    L <- ti
  if (np.order>=2)
    L <- 0.5*(3*ti*ti-1)
  if (np.order>=3)
    L <- 0.5*(5*ti^3-3*ti)
  if (np.order>=4)
    L <- 0.125*(35*ti^4-30*ti^2+3)
  if (np.order>=5)
    L <- 0.125*(63*ti^5-70*ti^3+15*ti)
  if (np.order>=6)
    L <- (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)
  if (np.order>=7)
    L <- (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)
  if (np.order>=8)
    L <- (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)
  if (np.order>=9)
    L <- (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)
 
  return(L);
}




LegendreN <-function( t, mu,np.order, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- np.order-1
  L <- 1*mu[1]
  if (np.order>=1)
    L <- ti*mu[2]
  if (np.order>=2)
    L <- 0.5*(3*ti*ti-1)*mu[3]
  if (np.order>=3)
    L <- 0.5*(5*ti^3-3*ti)*mu[4]
  if (np.order>=4)
    L <- 0.125*(35*ti^4-30*ti^2+3)*mu[5]
  if (np.order>=5)
    L <- 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6]
  if (np.order>=6)
    L <- (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)*mu[7]
  if (np.order>=7)
    L <- (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)*mu[8]
  if (np.order>=8)
    L <- (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)*mu[9]
  if (np.order>=9)
    L <- (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)*mu[10]
  
  return(L);
}

dLegendreN <-function( t, mu,np.order, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- np.order-1
  L <- 1*mu[1]
  if (np.order>=1)
    L <- mu[2]
  if (np.order>=2)
    L <- 0.5*(6*ti)*mu[3]
  if (np.order>=3)
    L <- 0.5*(15*ti^2-3)*mu[4]
  if (np.order>=4)
    L <- 0.125*(35*4*ti^3-60*ti)*mu[5]
  if (np.order>=5)
    L <- 0.125*(63*5*ti^4-210*ti^2+15)*mu[6]
  if (np.order>=6)
    L <- (1/16)*(231*6*ti^5-4*315*ti^3+210*ti)*mu[7]
  
  return(L);
}

Legendre.modelve <- function(t, np.order, tmin = NULL, tmax = NULL)
{
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- NA
  L <- 1;
  if (np.order >= 2)
    L <- 0.5 * (6 * ti) 
  if (np.order >= 3)
    L <- 0.5 * (15 * ti ^ 2 - 3)
  if (np.order >= 4)
    L <- 0.125 * (35 * 4 * ti ^ 3 - 60 * ti)
  if (np.order >= 5)
    L <- 0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)
  if (np.order >= 6)
    L <- (1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *
                       ti)
  if (np.order >= 7)
    L <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *
                       ti ^ 2 - 35) 
  return(L);
}