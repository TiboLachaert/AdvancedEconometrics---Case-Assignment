####################################################################
##################               OLS              ##################
####################################################################
OLS_own = function (y, x) 
{
  n  <- length(y)
  k  <- ncol(x)
  df <- n-k
  
  ## Run OLS
  xy     <- t(x)%*%y
  xxi    <- t(x)%*%x
  coefs  <- as.vector(solve(xxi)%*%xy)

  yhat   <- as.vector(x%*%coefs)
  res    <- y-yhat
  sigma2 <- as.vector(t(res)%*%res/df)
  
  stdvs  <- sqrt(sigma2)*sqrt(diag(solve(xxi)))
  tstats <- coefs/stdvs
  pvals  <- 2*(1-pt(abs(tstats),df))

  ## Save output
  names(coefs) <- colnames(x)
  
  coefs  <- round(coefs,3)
  stdvs  <- round(stdvs,3)
  tstats <- round(tstats,3)
  pvals  <- round(pvals,3)
  
  out = rbind(coefs, stdvs, tstats, pvals)
  out = t(out)
  return(out)
}

####################################################################
##################       Fixed Effects (FE)       ##################
####################################################################
FE_own = function (y, x) 
{
  T  <- dim(x)[1]
  N  <- dim(x)[2]
  K  <- dim(x)[3]
  df <- N*T-K
  
  ## Run FE
  Md  <- diag(T) - 1/T
  XDX <- matrix(data = 0, nrow = K, ncol = K)
  XDy <- matrix(data = 0, nrow = K, ncol = 1)
  
  for (i in 1:N){
    xi   <- x[, i, ]
    yi   <- y[, i]
    
    xdxi <- t(xi)%*%Md%*%xi
    XDX  <- XDX + xdxi
    
    xdyi <- t(xi)%*%Md%*%yi
    XDy  <- XDy + xdyi
  }

  coefs  <- as.vector(Ginv(XDX)%*%XDy)
    
  y_long <- matrix(y, ncol = 1)
  x_long <- matrix(x, ncol = K)
  
  yhat   <- as.vector(x_long%*%coefs)
  res    <- y_long - yhat
  sigma2 <- as.vector(t(res)%*%res/df)

  stdvs  <- sqrt(sigma2)*sqrt(diag(Ginv(XDX)))
  tstats <- coefs/stdvs
  pvals  <- 2*(1-pt(abs(tstats),df))
  
  ## Save output
  
  names(coefs) <- colnames(x[1, ,])
  
  coefs  <- round(coefs,3)
  stdvs  <- round(stdvs,3)
  tstats <- round(tstats,3)
  pvals  <- round(pvals,3)
  
  out = rbind(coefs, stdvs, tstats, pvals)
  out = t(out)
  return(out)
}


