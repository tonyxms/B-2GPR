#Implement Algorithm 2 (Spherical Hamiltonian Monte Carlo) with adaptive step size tuning.
=function(q_cur, u_cur, du_cur, U,eps,eps_loghn,eps_An,eps_mu, L=5,iter,Nadpt){
  # initialization
  q = q_cur; len_q = length(q)
  u = u_cur; du = du_cur
  
  # sample velocity
  v = rnorm(len_q) # standard multi-normal
  v <- v - q*drop(t(q)%*%v) # force v to lie in tangent space of sphere
  
  # Evaluate potential and kinetic energies at start of trajectory
  E_cur = u + .5*sum(v^2)
  
  randL = ceiling(L*runif(1))
  
  q_record=q
  v_record=v
  
  
  
  for (l in 1:randL){
    leap=leapfrog(q,v,u,du,U,eps)
    q=leap$q
    v=leap$v
    u=leap$u
    du=leap$du
    
    q_record=cbind(q_record,q,deparse.level = 0)
    v_record=cbind(v_record,v,deparse.level = 0)
    
    
    #		if (abs(t(q)%*%v)>1e-6){ # calibrate direction possibly deviated by error accumulation 
    #			v <- v - q*(t(q)%*%v)
    #			cat('Direction calibrated!')
    #		}
    
    #if (drop(t(q)%*%q_cur)<0) break
    val <- drop(t(q) %*% q_cur)
    tryCatch({
      if (val < 0) break
    }, error = function(e) {
      cat("q: ", q, "\n")
      cat("q_cur: ", q_cur, "\n")
      cat("val", val, "\n")
      stop(e)
    })
  }
  E_prp = u + .5*sum(v^2)
  logAP = -E_prp + E_cur
  
  if( is.finite(logAP)&&(log(runif(1))<min(0,logAP)) ){
    q = q
    u = u
    du = du
    Ind = 1
  }
  else {
    q = q_cur
    u = u_cur
    du = du_cur
    Ind = 0
  }
  
  if (iter<=Nadpt){
    eps_set=dual_avg(eps,eps_loghn,eps_An,eps_mu,iter,an=exp(min(c(0,logAP))))
    eps=eps_set$eps
    eps_loghn=eps_set$eps_loghn
    eps_An=eps_set$eps_An
  }
  
  if (iter==Nadpt){
    eps=exp(eps_loghn)
  }
  
  return (list(q=q,u=u,du=du,Ind=Ind,eps=eps,eps_loghn=eps_loghn,eps_An=eps_An,q_record=q_record,v_record=v_record))
  
}