#Implement the leapfrog integrator (for-loop part) used in Algorithm 2 (Spherical Hamiltonian Monte Carlo).
leapfrog=function(q,v,u,du,U,eps){
  len_q=length(q)
  # Make a half step for velocity
  g = c(du,0) - q*drop(t(q[-len_q])%*%du)
  v = v - eps/2 * g
  
  # Make a full step for the position
  q0 = q; v_nom = sqrt(sum(v^2))
  cosvt = cos(v_nom*eps); sinvt=sin(v_nom*eps)
  q = q0*cosvt + v/v_nom*sinvt
  v <- -q0*v_nom*sinvt + v*cosvt
  
  u = U(q[-len_q])
  du = U(q[-len_q],dd=T); 
  g = c(du,0) - q*drop((t(q[-len_q])%*%du))
  
  # Make a half step for velocity
  v = v - eps/2 * g
  
  #		if (abs(t(q)%*%v)>1e-6){ # calibrate direction possibly deviated by error accumulation 
  #			v <- v - q*(t(q)%*%v)
  #			cat('Direction calibrated!')
  #		}
  return(list(q=q,v=v,u=u,du=du))
}