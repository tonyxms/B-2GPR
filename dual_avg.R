#Implement the dual-averaging scheme used to adaptively tune the step size during the burn-in phase.
dual_avg=function(eps,eps_loghn,eps_An,eps_mu,iter,an,gamma=0.05,n0=10,kappa=0.75,a0=0.65){
  eps_An=(1-1/(iter+n0))*eps_An+(a0-an)/(iter+n0)
  logh=eps_mu-sqrt(iter)/gamma*eps_An
  eps_loghn=iter^(-kappa)*logh+(1-iter^(-kappa))*eps_loghn
  eps=exp(logh)
  return(list(eps=eps,eps_loghn=eps_loghn,eps_An=eps_An))
}