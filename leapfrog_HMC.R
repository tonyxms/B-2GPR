#Implement the leapfrog integrator (for-loop part) used in Algorithm 1 (Hamiltonian Monte Carlo).
leapfrog_HMC=function(q,v,u,du,U,eps){
  len_q=length(q)
  # Make a half step for velocity
  v = v - eps/2 * du
  
  # Make a full step for the position
  q = q + eps * v
  
  u = U(q)
  du = U(q,T)
  
  
  # Make a half step for velocity
  v = v - eps/2 * du
  
  #		if (abs(t(q)%*%v)>1e-6){ # calibrate direction possibly deviated by error accumulation 
  #			v <- v - q*(t(q)%*%v)
  #			cat('Direction calibrated!')
  #		}
  return(list(q=q,v=v,u=u,du=du))
}