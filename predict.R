#Use the trained Bayesian bridge model to generate predictions.
bayesian_bridge_predict=function(x,y,G,newx,newG,Gauss_bayes){
  beta=Gauss_bayes$beta
  omega=Gauss_bayes$omega
  tau_square=Gauss_bayes$tau_square
  eta=Gauss_bayes$eta
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.matrix(newx)) newx <- as.matrix(newx)
  n <- nrow(x)
  p <- ncol(x)
  m <- nrow(newx)
  K_inv=Gauss_bayes$K_inv
  ### Compute the pairwised distance between newx and x
  mzmax <- m*n
  D <- matrix(0,nrow=mzmax,ncol=p)
  for (k in 1:m) {
    lines_ix <- ((k-1)*n+1):(k*n)
    D[lines_ix,] <- (newx[rep(k,n),]-x)^2
  }
  newK <- matrix(0,nrow=m,ncol=n)
  for (i in 1:m) {
    line_ix <- ((i-1)*n+1):(i*n)
    newK[i, ]= exp(-D[line_ix, ]%*%omega^2)
  }
  predy=newG%*%beta+newK%*%K_inv%*%(y-G%*%beta)
  return(drop(predy))
}