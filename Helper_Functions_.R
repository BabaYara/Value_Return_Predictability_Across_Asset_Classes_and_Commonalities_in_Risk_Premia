myrename <- function(x,y){
  names(x) <- y
  return(x)
}

formatN2 <- function(x, y = 3){
  sprintf(paste("%.",y,"f", sep=""), round(x, y))
}


ComRet <- function(R, RF, H, type = 1){  
  
  data <- merge(R, RF, join='left')
  #data <- data[complete.cases]
  data$Add <- data[,1] + data[, 2]
  if(type == 1){
    R_c01 = rollapply(data$Add, width = H, FUN = function(x){prod(1 +  x)}, align = "left")
    R_c02 = rollapply(data$RF, width = H, FUN = function(x){prod(1 +  x)}, align = "left")
    R_c   = R_c01 - R_c02
  }
  
  if(type == 2){
    R_c = rollapply(data$R, H, function(x){
      sum(x[, 1])})}
  return(R_c)
}

olsNW02 <- function(reg = reg,lag){
  # %Syntax: olsnwsdbeta(Y,X,Const,lag)
  # %
  # %Calculates the Newey-West standard deviation of the
  # %estimated OLS betas an OLS regression of Y on X, with
  # %lag lengt lag. 
  # %Y may contain more than one variable; Beta assumes that
  # %each column in Y refers to a different variable. If you
  # %wantto include a constant in your regression, then the
  # %variable Const should be equal to 1; If it is zero then
  # %Beta assumes that you do not want to include a constant.
  
  
  Y <- matrix(as.matrix(reg$model)[, 1])
  if(ncol(reg$model) > 2){
    X <- as.matrix(reg$model)[, 2:ncol(reg$model)]
  } else {
    X <- matrix(as.matrix(reg$model)[, 2:ncol(reg$model)])
  }
  
  Const <- 0
  if(length(names(reg$coefficients)) > I(ncol(reg$model)-1)){
    Const = 1
  } 
  
  Y.dim = dim(Y)
  X.dim = dim(X)
  
  if(Const == 1){
    c = array(1, c(Y.dim[1],1))
    X = cbind(c, X)
    X.dim[2] = X.dim[2] + 1
  }
  
  Omega = array( 0, c(Y.dim[2]*X.dim[2],Y.dim[2]*X.dim[2]))
  b = solve(t(X)%*%(X))%*%(t(X)%*%(Y))
  e = Y - X%*%b
  d = 0.5
  
  for(k in 0:lag){
    if(k > 0){
      d = 1
    }
    
    for(i in 1:Y.dim[2]){
      Zi = array(0, c(Y.dim[1], X.dim[2]))
      
      for(ii in 1:X.dim[2]){
        Zi[ ,ii] = e[ ,i]*X[ ,ii]
      }
      
      for(j in 1:Y.dim[2]){
        Z1j = array(0, c(Y.dim[1], X.dim[2]))
        Z2j = array(0, c(Y.dim[1], X.dim[2]))
        
        for(jj in 1:X.dim[2]){
          Z1j[I(1+k):Y.dim[1], jj] = e[1:I(Y.dim[1]-k), j]*X[1:I(Y.dim[1]-k), jj]
          Z2j[1:I(Y.dim[1]-k), jj] = e[I(1+k):Y.dim[1], j]*X[I(1+k):Y.dim[1], jj] 
        }
        
        Omega[I(I(i-1)*X.dim[2] + 1):I(i*X.dim[2]), I(I(j-1)*X.dim[2]+1):I(j*X.dim[2])]  = Omega[
          I(I(i-1)*X.dim[2]+1):I(i*X.dim[2]), I(I(j-1)*X.dim[2]+1):I(j*X.dim[2])] + (d*(1-k/(1+lag)))*(
            solve(t(X)%*%X)%*%(t(Zi)%*%Z1j)%*%solve(t(X)%*%X) + solve(t(X)%*%X)%*%t(t(Z2j)%*%Zi)%*%solve(t(X)%*%X))
        
      }
    }
  }
  
  V = Omega*(Y.dim[1]/(Y.dim[1]-X.dim[2]))
  
  var = diag(V)
  
  se = sqrt(var) 
  return(V)
}


hodrickse <- function(reg,lags){
  # single equation Hodrick (1992) standard errors
  # construct time series of (uncorrelated) et+1 from regression of ln rt+1 on constant
  
  Y <- matrix(as.matrix(reg$model)[, 1])
  X <- X <- as.matrix(reg$model)[, 2:ncol(reg$model)]
  #  if(ncol(reg$model) > 2){
  #    X <- as.matrix(reg$model)[, 2:ncol(reg$model)]
  #  } else {
  #    X <- matrix(as.matrix(reg$model)[, 2:ncol(reg$model)])
  #  }
  
  Tn = nrow(as.matrix(Y))
  model = lm(Y ~ 1);
  resid = as.matrix(model$residuals)
  
  T2 = nrow(as.matrix(X))
  X = cbind(1, X)
  N = ncol(X)
  A = solve(t(X)%*%X)
  SbT = matrix(0, nrow = N, ncol = N)
  for(t in (1+lags - 1):(Tn - 1)){
    if(is.null(dim(X[(t- (lags - 1)):t, ]))){
      wkt = resid[t+1 , 1]%*%X[(t- (lags - 1)):t, ];
      SbT = SbT + t(wkt)%*%wkt
    } else {
      wkt = resid[t+1 , 1]%*%apply(X[(t- (lags - 1)):t, ], 2, sum);
      SbT = SbT + t(wkt)%*%wkt
    }
    
  }
  omega = A%*%SbT%*%A;
  return(unname(omega))
}
 
