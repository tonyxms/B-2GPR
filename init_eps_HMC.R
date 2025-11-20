#Choose the initial step size for Algorithm 1 (Hamiltonian Monte Carlo).
init_eps_HMC=function(q,u,du,U){
  q0=q
  u0=u
  du0=du
  len_q0=length(q0)
  v0=rnorm(len_q0)
  v0 = v0 - q*drop(t(q)%*%v0)
  
  E_cur=u0 + .5*sum(v0^2)
  eps=1
  
  leap=leapfrog_HMC(q0,v0,u0,du0,U,eps)
  q=leap$q
  v=leap$v
  u=leap$u
  du=leap$du
  
  E_prp=u + .5*sum(v^2)
  logAP = -E_prp + E_cur
  a=2*(exp(logAP)>0.5)-1
  while (a*logAP>-a*log(2)){
    eps=eps*(2^a)
    leap=leapfrog_HMC(q0,v0,u0,du0,U,eps)
    q=leap$q
    v=leap$v
    u=leap$u
    du=leap$du
    E_prp=u + .5*sum(v^2)
    logAP = -E_prp + E_cur
  }
  return (eps)
  
}