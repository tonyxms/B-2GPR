#Implement Algorithm 3 (sampling procedure for sph-B^2GPR) for models with a constant fixed-effect term (scalar).
bayesian_bridge_sphmc_adp_convergence2=function(x,y,prior_var,G,n_in,n_adpt,n_out_min,n_out_max,q_norm,r){
  
  n=nrow(x)
  d=ncol(x)
  d_beta=ncol(G)
  df_tau_square=prior_var[5]
  a_eta=prior_var[6]
  b_eta=prior_var[7]
  q_norm=q_norm
  r_beta=r[1]
  r_omega=r[2]
  r_beta_c=r[1]
  r_omega_c=r[2]
  
  
  beta=rmvnorm(1,mean=rep(0,d_beta),sigma=sqrt(10)*diag(d_beta))
  beta=drop(beta)
  omega=rnorm(d,mean=0,sd=sqrt(10))
  omega=abs(omega)
  #omega=as.matrix(omega)
  tau_square=rinvchisq(1,df_tau_square)
  eta=rgamma(1,shape=a_eta,rate=b_eta)
  
  beta_c=rmvnorm(1,mean=rep(0,d_beta),sigma=sqrt(10)*diag(d_beta))
  beta_c=drop(beta_c)
  omega_c=rnorm(d,mean=0,sd=sqrt(10))
  omega_c=abs(omega_c)
  #omega=as.matrix(omega)
  tau_square_c=rinvchisq(1,df_tau_square)
  eta_c=rgamma(1,shape=a_eta,rate=b_eta)
  
  I=diag(n)
  one=drop(rep(1,n))
  ### Compute the Distance matrix 
  mzmax <- n*(n-1)/2
  D <- matrix(0,nrow=mzmax,ncol=d)
  ll <- 0
  for (k in 1:(n-1)){
    ll <- tail(ll,1)+(1:(n-k))
    D[ll,] <- (x[rep(k,n-k),]-x[(k+1):n,])^2
  }
  
  Dis2=D%*%(omega^2)
  Dis <- exp(-Dis2)
  K <- matrix(0,nrow=n,ncol=n)
  K[lower.tri(K)] <- Dis
  K <- t(K)+K
  diag(K) <- one
  
  K_inv=solve(K+eta*I)
  
  Dis2_c=D%*%(omega_c^2)
  Dis_c <- exp(-Dis2_c)
  K_c <- matrix(0,nrow=n,ncol=n)
  K_c[lower.tri(K)] <- Dis_c
  K_c <- t(K_c)+K_c
  diag(K_c) <- one
  
  K_inv_c=solve(K_c+eta_c*I)
  
  #energy function for beta
  U_beta=function(beta,sigma,sigma_inv,beta_hat,y=y,dd=F){
    if (dd==F){
      return (drop((1/2)*t(beta-beta_hat)%*%sigma_inv%*%(beta-beta_hat)))
    }
    else if (dd==T){
      return (drop(sigma_inv%*%(beta-beta_hat)))
    }
    else stop("Wrong choice!")
  }
  
  Ut_beta=function(theta_beta,sigma,sigma_inv,beta_hat,y=y,q_norm,r_beta,dd=F){
    beta=sign(theta_beta)*abs(theta_beta)^(2/q_norm)*r_beta
    if(dd==F){
      ut=U_beta(beta,sigma,sigma_inv,beta_hat,y)
      return (drop(ut))
    }
    else if(dd==T){
      du=U_beta(beta,sigma,sigma_inv,beta_hat,y,T)
      dut=du*r_beta*2/q_norm*abs(theta_beta)^(2/q_norm-1)
      return (drop(dut))
    }
    else stop("Wrong choice!")
  }
  
  #energy function for omega
  U_omega=function(omega,D,eta,tau_square,beta,y=y,dd=F){
    Dis2_omega=D%*%(omega^2)
    Dis_omega <- exp(-Dis2_omega)
    K_omega <- matrix(0,nrow=n,ncol=n)
    K_omega[lower.tri(K)] <- Dis_omega
    K_omega <- t(K_omega)+K_omega
    diag(K_omega) <- one
    K_inv_omega=solve(K_omega+eta*I)
    #cat("1","\n")
    if (dd==F){
      return (drop(0.5*as.numeric(determinant(K_omega+eta*I,logarithm=TRUE)$modulus)+(1/(2*tau_square))*t((y-G%*%beta))%*%K_inv_omega%*%(y-G%*%beta)))
    }
    else if (dd==T){
      omega_grad=matrix(0,nrow=d,ncol=1)
      for (i in 1:d){
        Dis_m_omega=-D[ ,i]*2*omega[i]
        K_m_omega=matrix(0,nrow=n,ncol=n)
        K_m_omega[lower.tri(K_m_omega)] = Dis_m_omega
        K_m_omega = t(K_m_omega)+K_m_omega
        diag(K_m_omega)=0
        K_p_omega=K_omega*K_m_omega
        omega_grad[i]=0.5*tr(K_inv_omega%*%K_p_omega)-(1/(2*tau_square))*t(y-G%*%beta)%*%K_inv_omega%*%K_p_omega%*%K_inv_omega%*%(y-G%*%beta)
      }
      return (drop(omega_grad))
    }
    else stop("Wrong choice!")
  }
  
  Ut_omega=function(theta_omega,D,eta,tau_square,beta,y=y,q_norm,r_omega,dd=F){
    omega=sign(theta_omega)*abs(theta_omega)^(2/q_norm)*r_omega
    if(dd==F){
      ut=U_omega(omega,D,eta,tau_square,beta,y)
      return (drop(ut))
    }
    else if(dd==T){
      du=U_omega(omega,D,eta,tau_square,beta,y,T)
      dut=du*r_omega*2/q_norm*abs(theta_omega)^(2/q_norm-1)
      return (drop(dut))
    }
    else stop("Wrong choice!")
  }
  
  
  U_eta=function(eta,K,tau_square,beta,y=y,dd=F){
    K_inv_eta=solve(K+eta*I)
    if (dd==F) {
      return (drop((1-a_eta)*log(eta)+b_eta*eta+0.5*as.numeric(determinant(K+eta*I,logarithm=TRUE)$modulus)+(1/(2*tau_square))*drop(t((y-G%*%beta))%*%K_inv_eta%*%(y-G%*%beta))))
    }
    else if (dd==T) {
      return (drop((1-a_eta)/eta+b_eta+0.5*tr(K_inv_eta)-(1/(2*tau_square))*drop(t(y-G%*%beta)%*%K_inv_eta%*%K_inv_eta%*%(y-G%*%beta))))
    }
    else stop("Wrong choice!")
  }
  
  U_r_beta=function(theta_beta,r_beta,K_inv,tau_square,q_norm,y=y,dd=F){
    if (dd==F){
      return (drop(1/(2*tau_square)*t(y-G%*%(sign(theta_beta)*abs(theta_beta)^(2/q_norm)*r_beta))%*%K_inv%*%(y-G%*%(sign(theta_beta)*abs(theta_beta)^(2/q_norm)*r_beta))))
    }
  }
  
  U_r_omega=function(theta_omega,r_omega,D,eta,beta,tau_square,q_norm,y=y,dd=F){
    omega=sign(theta_omega)*abs(theta_omega)^(2/q_norm)*r_omega
    Dis2_omega=D%*%(omega^2)
    Dis_omega <- exp(-Dis2_omega)
    K_omega <- matrix(0,nrow=n,ncol=n)
    K_omega[lower.tri(K)] <- Dis_omega
    K_omega <- t(K_omega)+K_omega
    diag(K_omega) <- one
    K_inv_omega=solve(K_omega+eta*I)
    if (dd==F){
      return (drop(0.5*as.numeric(determinant(K_omega+eta*I,logarithm=TRUE)$modulus)+(1/(2*tau_square))*t((y-G%*%beta))%*%K_inv_omega%*%(y-G%*%beta)))
    }
    
  }
  
  
  accp_beta=0
  accp_omega=0
  accp_eta=0
  accp_r_beta=0
  accp_r_omega=0
  
  theta_beta=beta/(1.05*norm(beta,"2"))
  theta_beta_til=c(theta_beta,sqrt(1-sum(theta_beta^2)))
  theta_omega=omega/(1.05*norm(omega,"2"))
  theta_omega_til=c(theta_omega,sqrt(1-sum(theta_omega^2)))
  
  sigma_inv=(1/tau_square)*t(G)%*%K_inv%*%G
  sigma=solve(sigma_inv)
  beta_hat=(1/tau_square)*sigma%*%(t(G)%*%K_inv%*%y)
  
  u_beta=Ut_beta(theta_beta,sigma,sigma_inv,beta_hat,y=y,q_norm,r_beta)
  du_beta=Ut_beta(theta_beta,sigma,sigma_inv,beta_hat,y=y,q_norm,r_beta,T)
  
  eps_beta=init_eps(theta_beta_til,u_beta,du_beta,function(theta_beta,dd=F) Ut_beta(theta_beta,sigma,sigma_inv,beta_hat,y=y,q_norm,r_beta,dd=dd))
  
  eps_beta_mu=log(10*eps_beta)
  eps_beta_loghn=0
  eps_beta_An=0
  
  u_omega=Ut_omega(theta_omega,D,eta,tau_square,beta,y=y,q_norm,r_omega)
  du_omega=Ut_omega(theta_omega,D,eta,tau_square,beta,y=y,q_norm,r_omega,T)
  
  eps_omega=init_eps(theta_omega_til,u_omega,du_omega,function(theta_omega,dd=F) Ut_omega(theta_omega,D,eta,tau_square,beta,y=y,q_norm,r_omega,dd=dd))
  
  eps_omega_mu=log(10*eps_omega)
  eps_omega_loghn=0
  eps_omega_An=0
  
  beta_record=beta
  omega_record=omega
  tau_square_record=tau_square
  eta_record=eta
  eps_beta_record=eps_beta
  eps_omega_record=eps_omega
  r_beta_record=r_beta
  r_omega_record=r_omega
  
  accp_beta_c=0
  accp_omega_c=0
  accp_eta_c=0
  accp_r_beta_c=0
  accp_r_omega_c=0
  
  theta_beta_c=beta_c/(1.05*norm(beta_c,"2"))
  theta_beta_til_c=c(theta_beta_c,sqrt(1-sum(theta_beta_c^2)))
  theta_omega_c=omega_c/(1.05*norm(omega_c,"2"))
  theta_omega_til_c=c(theta_omega_c,sqrt(1-sum(theta_omega_c^2)))
  
  sigma_inv_c=(1/tau_square_c)*t(G)%*%K_inv_c%*%G
  sigma_c=solve(sigma_inv_c)
  beta_hat_c=(1/tau_square_c)*sigma_c%*%(t(G)%*%K_inv_c%*%y)
  
  u_beta_c=Ut_beta(theta_beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,q_norm,r_beta_c)
  du_beta_c=Ut_beta(theta_beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,q_norm,r_beta_c,T)
  
  eps_beta_c=init_eps(theta_beta_til_c,u_beta_c,du_beta_c,function(theta_beta_c,dd=F) Ut_beta(theta_beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,q_norm,r_beta_c,dd=dd))
  
  eps_beta_mu_c=log(10*eps_beta_c)
  eps_beta_loghn_c=0
  eps_beta_An_c=0
  
  u_omega_c=Ut_omega(theta_omega_c,D,eta_c,tau_square_c,beta_c,y=y,q_norm,r_omega_c)
  du_omega_c=Ut_omega(theta_omega_c,D,eta_c,tau_square_c,beta_c,y=y,q_norm,r_omega_c,T)
  
  eps_omega_c=init_eps(theta_omega_til_c,u_omega_c,du_omega_c,function(theta_omega_c,dd=F) Ut_omega(theta_omega_c,D,eta_c,tau_square_c,beta_c,y=y,q_norm,r_omega_c,dd=dd))
  
  eps_omega_mu_c=log(10*eps_omega_c)
  eps_omega_loghn_c=0
  eps_omega_An_c=0
  
  beta_record_c=beta_c
  omega_record_c=omega_c
  tau_square_record_c=tau_square_c
  eta_record_c=eta_c
  eps_beta_record_c=eps_beta_c
  eps_omega_record_c=eps_omega_c
  r_beta_record_c=r_beta_c
  r_omega_record_c=r_omega_c
  
  Iter=1
  convergence=FALSE
  
  q_record=NULL
  v_record=NULL
  
  q_record_c=NULL
  v_record_c=NULL
  
  while (1){
    
    # Use Spherical HMC to get sample beta
    
    sigma_inv=(1/tau_square)*t(G)%*%K_inv%*%G
    sigma=solve(sigma_inv)
    beta_hat=(1/tau_square)*sigma%*%(t(G)%*%K_inv%*%y)
    
    
    u_beta=Ut_beta(theta_beta,sigma,sigma_inv,beta_hat,y=y,q_norm,r_beta)
    du_beta=Ut_beta(theta_beta,sigma,sigma_inv,beta_hat,y=y,q_norm,r_beta,T)
    
    samp_beta=adp_SphHMC(theta_beta_til,u_beta,du_beta,function(theta_beta,dd=F) Ut_beta(theta_beta,sigma,sigma_inv,beta_hat,y=y,q_norm,r_beta,dd=dd),eps=eps_beta,eps_loghn=eps_beta_loghn,eps_An=eps_beta_An,eps_mu=eps_beta_mu,L=n_in,iter=Iter,Nadpt=n_adpt)
    theta_beta_til=samp_beta$q
    accp_beta=accp_beta+samp_beta$Ind
    theta_beta=theta_beta_til[-(d_beta+1)]
    #beta=sign(theta_beta)*abs(theta_beta)^(2/q_norm)*r_beta
    eps_beta=samp_beta$eps
    eps_beta_loghn=samp_beta$eps_loghn
    eps_beta_An=samp_beta$eps_An
    
    #beta_record=cbind(beta_record,beta,deparse.level = 0)
    eps_beta_record=c(eps_beta_record,eps_beta)
    
    u_r_beta=U_r_beta(theta_beta,r_beta,K_inv,tau_square,q_norm,y=y)
    
    samp_r_beta=metropolis(r_beta,u_r_beta,function(r_beta,dd=F) U_r_beta(theta_beta,r_beta,K_inv,tau_square,q_norm,y=y,dd=dd),eps=.02,L=n_in)
    
    r_beta=samp_r_beta$q
    accp_r_beta=accp_r_beta+samp_r_beta$Ind
    
    r_beta_record=c(r_beta_record,r_beta)
    
    beta=sign(theta_beta)*abs(theta_beta)^(2/q_norm)*r_beta
    
    beta_record=cbind(beta_record,beta,deparse.level = 0)
    
    
    # Use Spherical HMC to get sample omega
    
    u_omega=Ut_omega(theta_omega,D,eta,tau_square,beta,y=y,q_norm,r_omega)
    du_omega=Ut_omega(theta_omega,D,eta,tau_square,beta,y=y,q_norm,r_omega,T)
    
    samp_omega=adp_SphHMC(theta_omega_til,u_omega,du_omega,function(theta_omega,dd=F) Ut_omega(theta_omega,D,eta,tau_square,beta,y=y,q_norm,r_omega,dd=dd),eps=eps_omega,eps_loghn=eps_omega_loghn,eps_An=eps_omega_An,eps_mu=eps_omega_mu,L=n_in,iter=Iter,Nadpt=n_adpt)
    theta_omega_til=samp_omega$q
    accp_omega=accp_omega+samp_omega$Ind
    theta_omega=theta_omega_til[-(d+1)]
    #omega=sign(theta_omega)*abs(theta_omega)^(2/q_norm)*r_omega
    eps_omega=samp_omega$eps
    eps_omega_loghn=samp_omega$eps_loghn
    eps_omega_An=samp_omega$eps_An
    q_record=cbind(q_record,samp_omega$q_record,deparse.level = 0)
    v_record=cbind(v_record,samp_omega$v_record,deparse.level = 0)
    
    
    #omega_record=cbind(omega_record,omega,deparse.level = 0)
    eps_omega_record=c(eps_omega_record,eps_omega)
    
    
    
    
    u_r_omega=U_r_omega(theta_omega,r_omega,D,eta,beta,tau_square,q_norm,y=y)
    
    samp_omega=metropolis(r_omega,u_r_omega,function(r_omega,dd=F) U_r_omega(theta_omega,r_omega,D,eta,beta,tau_square,q_norm,y=y,dd=dd),eps=.02,L=n_in)
    
    r_omega=samp_omega$q
    accp_r_omega=accp_r_omega+samp_omega$Ind
    
    r_omega_record=c(r_omega_record,r_omega)
    
    omega=sign(theta_omega)*abs(theta_omega)^(2/q_norm)*r_omega
    
    omega_record=cbind(omega_record,omega,deparse.level = 0)
    
    
    Dis2=D%*%omega^2
    Dis <- exp(-Dis2)
    K <- matrix(0,nrow=n,ncol=n)
    K[lower.tri(K)] <- Dis
    K <- t(K)+K
    diag(K) <- one
    
    K_inv=solve(K+eta*I)
    
    # sample tau_square
    
    sn_square=t(y-G%*%beta)%*%K_inv%*%(y-G%*%beta)
    tau_square_hat=(1+sn_square)/(df_tau_square+n)
    tau_square=rinvchisq(1,df_tau_square+n,scale=tau_square_hat)
    
    tau_square_record=c(tau_square_record,tau_square)
    # sample eta
    
    u_eta=U_eta(eta,K,tau_square,beta,y=y)
    
    samp_eta=metropolis(eta,u_eta,function(eta,dd=F) U_eta(eta,K,tau_square,beta,y=y,dd=dd),eps=.02,L=n_in)
    eta=samp_eta$q
    accp_eta=accp_eta+samp_eta$Ind
    if (eta<1e-08){
      eta=1e-08
    }
    
    
    eta_record=c(eta_record,eta)
    
    K_inv=solve(K+eta*I)
    
    
    
    
    
    
    
    
    
    # Use Spherical HMC to get sample beta
    
    sigma_inv_c=(1/tau_square_c)*t(G)%*%K_inv_c%*%G
    sigma_c=solve(sigma_inv_c)
    beta_hat_c=(1/tau_square_c)*sigma_c%*%(t(G)%*%K_inv_c%*%y)
    
    u_beta_c=Ut_beta(theta_beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,q_norm,r_beta_c)
    du_beta_c=Ut_beta(theta_beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,q_norm,r_beta_c,T)
    
    samp_beta_c=adp_SphHMC(theta_beta_til_c,u_beta_c,du_beta_c,function(theta_beta_c,dd=F) Ut_beta(theta_beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,q_norm,r_beta_c,dd=dd),eps=eps_beta_c,eps_loghn=eps_beta_loghn_c,eps_An=eps_beta_An_c,eps_mu=eps_beta_mu_c,L=n_in,iter=Iter,Nadpt=n_adpt)
    theta_beta_til_c=samp_beta_c$q
    accp_beta_c=accp_beta_c+samp_beta_c$Ind
    theta_beta_c=theta_beta_til_c[-(d_beta+1)]
    #beta_c=sign(theta_beta_c)*abs(theta_beta_c)^(2/q_norm)*r_beta
    eps_beta_c=samp_beta_c$eps
    eps_beta_loghn_c=samp_beta_c$eps_loghn
    eps_beta_An_c=samp_beta_c$eps_An
    
    #beta_record_c=cbind(beta_record_c,beta_c,deparse.level = 0)
    eps_beta_record_c=c(eps_beta_record_c,eps_beta_c)
    
    u_r_beta_c=U_r_beta(theta_beta_c,r_beta_c,K_inv_c,tau_square_c,q_norm,y=y)
    
    samp_r_beta_c=metropolis(r_beta_c,u_r_beta_c,function(r_beta_c,dd=F) U_r_beta(theta_beta_c,r_beta_c,K_inv_c,tau_square_c,q_norm,y=y,dd=dd),eps=.02,L=n_in)
    
    r_beta_c=samp_r_beta_c$q
    accp_r_beta_c=accp_r_beta_c+samp_r_beta_c$Ind
    
    r_beta_record_c=c(r_beta_record_c,r_beta_c)
    
    beta_c=sign(theta_beta_c)*abs(theta_beta_c)^(2/q_norm)*r_beta_c
    
    beta_record_c=cbind(beta_record_c,beta_c,deparse.level = 0)
    
    
    
    # Use Spherical HMC to get sample omega
    
    u_omega_c=Ut_omega(theta_omega_c,D,eta_c,tau_square_c,beta_c,y=y,q_norm,r_omega_c)
    du_omega_c=Ut_omega(theta_omega_c,D,eta_c,tau_square_c,beta_c,y=y,q_norm,r_omega_c,T)
    
    samp_omega_c=adp_SphHMC(theta_omega_til_c,u_omega_c,du_omega_c,function(theta_omega_c,dd=F) Ut_omega(theta_omega_c,D,eta_c,tau_square_c,beta_c,y=y,q_norm,r_omega_c,dd=dd),eps=eps_omega_c,eps_loghn=eps_omega_loghn_c,eps_An=eps_omega_An_c,eps_mu=eps_omega_mu_c,L=n_in,iter=Iter,Nadpt=n_adpt)
    theta_omega_til_c=samp_omega_c$q
    accp_omega_c=accp_omega_c+samp_omega_c$Ind
    theta_omega_c=theta_omega_til_c[-(d+1)]
    #omega_c=sign(theta_omega_c)*abs(theta_omega_c)^(2/q_norm)*r_omega
    eps_omega_c=samp_omega_c$eps
    eps_omega_loghn_c=samp_omega_c$eps_loghn
    eps_omega_An_c=samp_omega_c$eps_An
    
    eps_omega_record_c=c(eps_omega_record_c,eps_omega_c)
    
    q_record_c=cbind(q_record_c,samp_omega_c$q_record,deparse.level = 0)
    v_record_c=cbind(v_record_c,samp_omega_c$v_record,deparse.level = 0)
    
    #omega_record_c=cbind(omega_record_c,omega_c,deparse.level = 0)
    
    u_r_omega_c=U_r_omega(theta_omega_c,r_omega_c,D,eta_c,beta_c,tau_square_c,q_norm,y=y)
    
    samp_omega_c=metropolis(r_omega_c,u_r_omega_c,function(r_omega_c,dd=F) U_r_omega(theta_omega_c,r_omega_c,D,eta_c,beta_c,tau_square_c,q_norm,y=y,dd=dd),eps=.02,L=n_in)
    
    r_omega_c=samp_omega_c$q
    accp_r_omega_c=accp_r_omega_c+samp_omega_c$Ind
    
    r_omega_record_c=c(r_omega_record_c,r_omega_c)
    
    omega_c=sign(theta_omega_c)*abs(theta_omega_c)^(2/q_norm)*r_omega_c
    
    omega_record_c=cbind(omega_record_c,omega_c,deparse.level = 0)
    
    Dis2_c=D%*%omega_c^2
    Dis_c <- exp(-Dis2_c)
    K_c <- matrix(0,nrow=n,ncol=n)
    K_c[lower.tri(K_c)] <- Dis_c
    K_c <- t(K_c)+K_c
    diag(K_c) <- one
    
    K_inv_c=solve(K_c+eta_c*I)
    
    # sample tau_square
    
    sn_square_c=t(y-G%*%beta_c)%*%K_inv_c%*%(y-G%*%beta_c)
    tau_square_hat_c=(1+sn_square_c)/(df_tau_square+n)
    tau_square_c=rinvchisq(1,df_tau_square+n,scale=tau_square_hat_c)
    
    tau_square_record_c=c(tau_square_record_c,tau_square_c)
    
    # sample eta
    
    u_eta_c=U_eta(eta_c,K_c,tau_square_c,beta_c,y=y)
    
    samp_eta_c=metropolis(eta_c,u_eta_c,function(eta_c,dd=F) U_eta(eta_c,K_c,tau_square_c,beta_c,y=y,dd=dd),eps=.02,L=n_in)
    
    eta_c=samp_eta_c$q
    accp_eta_c=accp_eta_c+samp_eta_c$Ind
    if (eta_c<1e-08){
      eta_c=1e-08
    }
    
    eta_record_c=c(eta_record_c,eta_c)
    
    K_inv_c=solve(K_c+eta_c*I)
    
    if(Iter%%100==0){
      cat('Iteration ',Iter,' completed!\n')
      cat('Acceptance Rate for beta', accp_beta/100,' \n')
      cat('Acceptance Rate for omega', accp_omega/100,' \n')
      cat('Acceptance Rate for eta', accp_eta/100,' \n')
      cat('Acceptance Rate for r_beta', accp_r_beta/100,' \n')
      cat('Acceptance Rate for r_omega', accp_r_omega/100,' \n')
      accp_beta=0
      accp_omega=0
      accp_eta=0
      accp_r_beta=0
      accp_r_omega=0
    }
    
    if(Iter%%100==0){
      cat('Iteration ',Iter,' completed!\n')
      cat('Acceptance Rate for beta_c', accp_beta_c/100,' \n')
      cat('Acceptance Rate for omega_c', accp_omega_c/100,' \n')
      cat('Acceptance Rate for eta_c', accp_eta_c/100,' \n')
      cat('Acceptance Rate for r_beta_c', accp_r_beta_c/100,' \n')
      cat('Acceptance Rate for r_omega_c', accp_r_omega_c/100,' \n')
      accp_beta_c=0
      accp_omega_c=0
      accp_eta_c=0
      accp_r_beta_c=0
      accp_r_omega_c=0
    }
    
    
    if (Iter>=n_out_min){
      if (Iter<n_out_max){
        if (Iter%%100==0){
          chain_beta=mcmc(t(beta_record))
          chain_beta_c=mcmc(t(beta_record_c))
          chain_omega=mcmc(t(abs(omega_record)))
          chain_omega_c=mcmc(t(abs(omega_record_c)))
          chain_tau_square=mcmc(tau_square_record)
          chain_tau_square_c=mcmc(tau_square_record_c)
          chain_eta=mcmc(eta_record)
          chain_eta_c=mcmc(eta_record_c)
          combinedchains_beta = mcmc.list(chain_beta, chain_beta_c)
          combinedchains_omega = mcmc.list(chain_omega, chain_omega_c)
          combinedchains_tau_square = mcmc.list(chain_tau_square, chain_tau_square_c)
          combinedchains_eta = mcmc.list(chain_eta, chain_eta_c)
          if (gelman.diag(combinedchains_beta)$psrf[1,1]<1.1 & gelman.diag(combinedchains_omega)$mpsrf<1.1 & gelman.diag(combinedchains_tau_square)$psrf[1,1]<1.1 & gelman.diag(combinedchains_eta)$psrf[1,1]<1.1){
            convergence=TRUE
            break
          }
        }
      }
    }
    
    if (Iter>=n_out_max){
      break
    }
    
    
    
    Iter=Iter+1
    
  }
  
  beta=median(beta_record[(n_adpt+1):length(tau_square_record)])
  omega=apply(omega_record[ ,(n_adpt+1):dim(omega_record)[2]],1,median)
  tau_square=median(tau_square_record[(n_adpt+1):length(tau_square_record)])
  eta=median(eta_record[(n_adpt+1):length(eta_record)])
  
  Dis2=D%*%omega^2
  Dis <- exp(-Dis2)
  K <- matrix(0,nrow=n,ncol=n)
  K[lower.tri(K)] <- Dis
  K <- t(K)+K
  diag(K) <- one
  
  K_inv=solve(K+eta*I)
  
  beta_c=median(beta_record_c[(n_adpt+1):length(tau_square_record)])
  omega_c=apply(omega_record_c[ ,(n_adpt+1):dim(omega_record_c)[2]],1,median)
  tau_square_c=median(tau_square_record_c[(n_adpt+1):length(tau_square_record_c)])
  eta_c=median(eta_record_c[(n_adpt+1):length(eta_record_c)])
  
  Dis2_c=D%*%omega_c^2
  Dis_c <- exp(-Dis2_c)
  K_c <- matrix(0,nrow=n,ncol=n)
  K_c[lower.tri(K_c)] <- Dis_c
  K_c <- t(K_c)+K_c
  diag(K_c) <- one
  
  K_inv_c=solve(K_c+eta_c*I)
  
  Gauss_bayes=list(beta=beta,omega=omega,tau_square=tau_square,eta=eta,K_inv=K_inv,K=K,eps_beta=eps_beta,eps_omega=eps_omega,beta_record=beta_record,omega_record=omega_record,tau_square_record=tau_square_record,eta_record=eta_record, eps_beta_record=eps_beta_record,eps_omega_record=eps_omega_record,beta_c=beta_c,omega_c=omega_c,tau_square_c=tau_square_c,eta_c=eta_c,K_inv_c=K_inv_c,K_c=K_c,eps_beta_c=eps_beta_c,eps_omega_c=eps_omega_c,beta_record_c=beta_record_c,omega_record_c=omega_record_c,tau_square_record_c=tau_square_record_c,eta_record_c=eta_record_c, eps_beta_record_c=eps_beta_record_c,eps_omega_record_c=eps_omega_record_c,convergence=convergence,q_record=q_record,v_record=v_record,q_record_c=q_record_c,v_record_c=v_record_c,r_beta=r_beta,r_beta_c=r_beta_c,r_omega=r_omega,r_omega_c=r_omega_c,r_beta_record=r_beta_record,r_beta_record_c=r_beta_record_c,r_omega_record=r_omega_record,r_omega_record_c=r_omega_record_c)
  class(Gauss_bayes)="Gaussian_bayes"
  return(Gauss_bayes)
  
}