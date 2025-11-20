# Metropolis algorithm for MCMC sampling
metropolis=function(q_cur,u_cur,U,eps,L=5){
  accepted=0
  q=q_cur
  len_q=length(q)
  u=u_cur
  
  randL = ceiling(runif(1)*L)
  for (i in 1:randL){
    stepp=rmvnorm(1,mean=rep(0,len_q),sigma=diag(len_q))
    stepp=drop(stepp)
    q_step=q+eps*abs(q)*stepp
    u_step=U(q_step,dd=F)
    log_alpha=u-u_step
    if (is.nan(log_alpha)) 
      log_alpha=-Inf
    if (log_alpha >=0.0) 
      alpha=1.0
    else 
      alpha=exp(log_alpha)
    if (runif(1)<=alpha) {
      q=q_step
      u=u_step
      accepted=accepted+1
    }
  }
  accepted=accepted/randL
  return(list(q=q,u=u,Ind=accepted))
}