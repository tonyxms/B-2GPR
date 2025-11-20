#Implement Algorithm 4 (sampling procedure for hmc-B^2GPR) for models with a constant fixed-effect term (scalar).
bayesian_bridge_sphmc_adp_hmc_convergence2=function(x,y,prior_var,G,n_in,n_adpt,n_out_min,n_out_max,l,rou){
  n=nrow(x)
  d=ncol(x)
  d_beta=ncol(G)
  a_beta=prior_var[1]
  b_beta=prior_var[2]
  a_omega=prior_var[3]
  b_omega=prior_var[4]
  df_tau_square=prior_var[5]
  a_eta=prior_var[6]
  b_eta=prior_var[7]
  #Set the initial R matrix
  if (l==0){
    R=diag(1)
  }
  else if (l==1){
    R=diag(c(1,rep(rou,d)))
  }
  else if (l==2){
    R=diag(c(1,rep(rou,d),rep(rou^2,d*(d+1)/2)))
  }
  else stop("Wrong choice for l!")
  p=dim(R)[1]
  
  
  v_beta_square=rinvgamma(1,shape=a_beta,scale=b_beta)
  v_omega_square=rinvgamma(1,shape=a_omega,scale=b_omega)
  beta=rmvnorm(1,mean=rep(0,d_beta),sigma=v_beta_square*R)
  beta=drop(beta)
  omega=rnorm(d,mean=0,sd=sqrt(v_omega_square))
  omega=abs(omega)
  #omega=as.matrix(omega)
  tau_square=rinvchisq(1,df_tau_square)
  eta=rgamma(1,shape=a_eta,rate=b_eta)
  
  v_beta_square_c=rinvgamma(1,shape=a_beta,scale=b_beta)
  v_omega_square_c=rinvgamma(1,shape=a_omega,scale=b_omega)
  beta_c=rmvnorm(1,mean=rep(0,d_beta),sigma=v_beta_square_c*R)
  beta_c=drop(beta_c)
  omega_c=rnorm(d,mean=0,sd=sqrt(v_omega_square_c))
  omega_c=abs(omega_c)
  #omega_c=as.matrix(omega_c)
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
  K_c[lower.tri(K_c)] <- Dis_c
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
  
  #energy function for omega
  U_omega=function(omega,D,eta,tau_square,beta,v_omega_square,y=y,dd=F){
    Dis2_omega=D%*%(omega^2)
    Dis_omega <- exp(-Dis2_omega)
    K_omega <- matrix(0,nrow=n,ncol=n)
    K_omega[lower.tri(K)] <- Dis_omega
    K_omega <- t(K_omega)+K_omega
    diag(K_omega) <- one
    K_inv_omega=solve(K_omega+eta*I)
    #cat("1","\n")
    if (dd==F){
      return (drop(sum(-omega^2/(2*v_omega_square))+0.5*as.numeric(determinant(K_omega+eta*I,logarithm=TRUE)$modulus)+(1/(2*tau_square))*t((y-G%*%beta))%*%K_inv_omega%*%(y-G%*%beta)))
    }
    else if (dd==T){
      omega_grad=matrix(0,nrow=d,ncol=1)
      for (i in 1:d){
        Dis_m_omega=-D[ ,i]*2*omega[i]
        K_m_omega=matrix(0,nrow=n,ncol=n)
        K_m_omega[lower.tri(K_m_omega)] <- Dis_m_omega
        K_m_omega <- t(K_m_omega)+K_m_omega
        diag(K_m_omega)=0
        K_p_omega=K_omega*K_m_omega
        omega_grad[i]=-omega[i]/v_omega_square+0.5*tr(K_inv_omega%*%K_p_omega)-(1/(2*tau_square))*t(y-G%*%beta)%*%K_inv_omega%*%K_p_omega%*%K_inv_omega%*%(y-G%*%beta)
      }
      return (drop(omega_grad))
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
  
  accp_beta=0
  accp_omega=0
  accp_eta=0
  
  sigma_inv=(1/tau_square)*t(G)%*%K_inv%*%G+(1/v_beta_square)*solve(R)
  sigma=solve(sigma_inv)
  beta_hat=(1/tau_square)*sigma%*%(t(G)%*%K_inv%*%y)
  
  u_beta=U_beta(beta,sigma,sigma_inv,beta_hat,y=y)
  du_beta=U_beta(beta,sigma,sigma_inv,beta_hat,y=y,T)
  
  eps_beta=init_eps_HMC(beta,u_beta,du_beta,function(beta,dd=F) U_beta(beta,sigma,sigma_inv,beta_hat,y=y,dd=dd))
  
  eps_beta_mu=log(10*eps_beta)
  eps_beta_loghn=0
  eps_beta_An=0
  
  u_omega=U_omega(omega,D,eta,tau_square,beta,v_omega_square,y=y)
  du_omega=U_omega(omega,D,eta,tau_square,beta,v_omega_square,y=y,T)
  
  eps_omega=init_eps_HMC(omega,u_omega,du_omega,function(omega,dd=F) U_omega(omega,D,eta,tau_square,beta,v_omega_square,y=y,dd=dd))
  
  eps_omega_mu=log(10*eps_omega)
  eps_omega_loghn=0
  eps_omega_An=0
  
  beta_record=beta
  omega_record=omega
  tau_square_record=tau_square
  eta_record=eta
  eps_beta_record=eps_beta
  eps_omega_record=eps_omega
  v_beta_square_record=v_beta_square
  v_omega_square_record=v_omega_square
  
  accp_beta_c=0
  accp_omega_c=0
  accp_eta_c=0
  
  sigma_inv_c=(1/tau_square_c)*t(G)%*%K_inv_c%*%G+(1/v_beta_square_c)*solve(R)
  sigma_c=solve(sigma_inv_c)
  beta_hat_c=(1/tau_square_c)*sigma_c%*%(t(G)%*%K_inv_c%*%y)
  
  u_beta_c=U_beta(beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y)
  du_beta_c=U_beta(beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,T)
  
  eps_beta_c=init_eps_HMC(beta_c,u_beta_c,du_beta_c,function(beta_c,dd=F) U_beta(beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,dd=dd))
  
  eps_beta_mu_c=log(10*eps_beta_c)
  eps_beta_loghn_c=0
  eps_beta_An_c=0
  
  u_omega_c=U_omega(omega_c,D,eta_c,tau_square_c,beta_c,v_omega_square_c,y=y)
  du_omega_c=U_omega(omega_c,D,eta_c,tau_square_c,beta_c,v_omega_square_c,y=y,T)
  
  eps_omega_c=init_eps_HMC(omega_c,u_omega_c,du_omega_c,function(omega_c,dd=F) U_omega(omega_c,D,eta_c,tau_square_c,beta_c,v_omega_square_c,y=y,dd=dd))
  
  eps_omega_mu_c=log(10*eps_omega_c)
  eps_omega_loghn_c=0
  eps_omega_An_c=0
  
  beta_record_c=beta_c
  omega_record_c=omega_c
  tau_square_record_c=tau_square_c
  eta_record_c=eta_c
  eps_beta_record_c=eps_beta_c
  eps_omega_record_c=eps_omega_c
  v_beta_square_record_c=v_beta_square_c
  v_omega_square_record_c=v_omega_square_c
  
  Iter=1
  convergence=FALSE
  
  while (1){
    
    sigma_inv=(1/tau_square)*t(G)%*%K_inv%*%G+(1/v_beta_square)*solve(R)
    sigma=solve(sigma_inv)
    beta_hat=(1/tau_square)*sigma%*%(t(G)%*%K_inv%*%y)
    
    # Use HMC to get sample beta
    
    u_beta=U_beta(beta,sigma,sigma_inv,beta_hat,y=y)
    du_beta=U_beta(beta,sigma,sigma_inv,beta_hat,y=y,T)
    
    samp_beta=adp_HMC(beta,u_beta,du_beta,function(beta,dd=F) U_beta(beta,sigma,sigma_inv,beta_hat,y=y,dd=dd),eps=eps_beta,eps_loghn=eps_beta_loghn,eps_An=eps_beta_An,eps_mu=eps_beta_mu,L=n_in,iter=Iter,Nadpt=n_adpt)
    beta=samp_beta$q
    accp_beta=accp_beta+samp_beta$Ind
    eps_beta=samp_beta$eps
    eps_beta_loghn=samp_beta$eps_loghn
    eps_beta_An=samp_beta$eps_An
    
    beta_record=cbind(beta_record,beta,deparse.level = 0)
    eps_beta_record=c(eps_beta_record,eps_beta)
    
    v_beta_square=rinvgamma(1,shape=a_beta+(p/2),scale=b_beta+0.5*t(beta)%*%solve(R)%*%beta)
    v_beta_square_record=c(v_beta_square_record,v_beta_square)
    
    # Use HMC to get sample omega
    
    u_omega=U_omega(omega,D,eta,tau_square,beta,v_omega_square,y=y)
    du_omega=U_omega(omega,D,eta,tau_square,beta,v_omega_square,y=y,T)
    
    samp_omega=adp_HMC(omega,u_omega,du_omega,function(omega,dd=F) U_omega(omega,D,eta,tau_square,beta,v_omega_square,y=y,dd=dd),eps=eps_omega,eps_loghn=eps_omega_loghn,eps_An=eps_omega_An,eps_mu=eps_omega_mu,L=n_in,iter=Iter,Nadpt=n_adpt)
    omega=abs(samp_omega$q)
    accp_omega=accp_omega+samp_omega$Ind
    eps_omega=samp_omega$eps
    eps_omega_loghn=samp_omega$eps_loghn
    eps_omega_An=samp_omega$eps_An
    
    omega_record=cbind(omega_record,omega,deparse.level = 0)
    eps_omega_record=c(eps_omega_record,eps_omega)
    
    Dis2=D%*%omega^2
    Dis <- exp(-Dis2)
    K <- matrix(0,nrow=n,ncol=n)
    K[lower.tri(K)] <- Dis
    K <- t(K)+K
    diag(K) <- one
    
    K_inv=solve(K+eta*I)
    
    v_omega_square=rinvgamma(1,shape=a_omega+(d/2),scale=b_omega+0.5*sum(omega^2))
    
    v_omega_square_record=c(v_omega_square_record,v_omega_square)
    
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
    
    
    
    
    
    
    
    sigma_inv_c=(1/tau_square_c)*t(G)%*%K_inv_c%*%G+(1/v_beta_square_c)*solve(R)
    sigma_c=solve(sigma_inv_c)
    beta_hat_c=(1/tau_square_c)*sigma_c%*%(t(G)%*%K_inv_c%*%y)
    
    # Use HMC to get sample beta
    
    u_beta_c=U_beta(beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y)
    du_beta_c=U_beta(beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,T)
    
    samp_beta_c=adp_HMC(beta_c,u_beta_c,du_beta_c,function(beta_c,dd=F) U_beta(beta_c,sigma_c,sigma_inv_c,beta_hat_c,y=y,dd=dd),eps=eps_beta_c,eps_loghn=eps_beta_loghn_c,eps_An=eps_beta_An_c,eps_mu=eps_beta_mu_c,L=n_in,iter=Iter,Nadpt=n_adpt)
    beta_c=samp_beta_c$q
    accp_beta_c=accp_beta_c+samp_beta_c$Ind
    eps_beta_c=samp_beta_c$eps
    eps_beta_loghn_c=samp_beta_c$eps_loghn
    eps_beta_An_c=samp_beta_c$eps_An
    
    beta_record_c=cbind(beta_record_c,beta_c,deparse.level = 0)
    eps_beta_record_c=c(eps_beta_record_c,eps_beta_c)
    
    v_beta_square_c=rinvgamma(1,shape=a_beta+(p/2),scale=b_beta+0.5*t(beta_c)%*%solve(R)%*%beta_c)
    v_beta_square_record_c=c(v_beta_square_record_c,v_beta_square_c)
    
    # Use HMC to get sample omega
    
    u_omega_c=U_omega(omega_c,D,eta_c,tau_square_c,beta_c,v_omega_square_c,y=y)
    du_omega_c=U_omega(omega_c,D,eta_c,tau_square_c,beta_c,v_omega_square_c,y=y,T)
    
    samp_omega_c=adp_HMC(omega_c,u_omega_c,du_omega_c,function(omega_c,dd=F) U_omega(omega_c,D,eta_c,tau_square_c,beta_c,v_omega_square_c,y=y,dd=dd),eps=eps_omega_c,eps_loghn=eps_omega_loghn_c,eps_An=eps_omega_An_c,eps_mu=eps_omega_mu_c,L=n_in,iter=Iter,Nadpt=n_adpt)
    omega_c=abs(samp_omega_c$q)
    accp_omega_c=accp_omega_c+samp_omega_c$Ind
    eps_omega_c=samp_omega_c$eps
    eps_omega_loghn_c=samp_omega_c$eps_loghn
    eps_omega_An_c=samp_omega_c$eps_An
    
    omega_record_c=cbind(omega_record_c,omega_c,deparse.level = 0)
    eps_omega_record_c=c(eps_omega_record_c,eps_omega_c)
    
    Dis2_c=D%*%omega_c^2
    Dis_c <- exp(-Dis2_c)
    K_c <- matrix(0,nrow=n,ncol=n)
    K_c[lower.tri(K_c)] <- Dis_c
    K_c <- t(K_c)+K_c
    diag(K_c) <- one
    
    K_inv_c=solve(K_c+eta_c*I)
    
    v_omega_square_c=rinvgamma(1,shape=a_omega+(d/2),scale=b_omega+0.5*sum(omega_c^2))
    
    v_omega_square_record_c=c(v_omega_square_record_c,v_omega_square_c)
    
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
      accp_beta=0
      accp_omega=0
      accp_eta=0
    }
    
    if(Iter%%100==0){
      cat('Iteration ',Iter,' completed!\n')
      cat('Acceptance Rate for beta_c', accp_beta_c/100,' \n')
      cat('Acceptance Rate for omega_c', accp_omega_c/100,' \n')
      cat('Acceptance Rate for eta_c', accp_eta_c/100,' \n')
      accp_beta_c=0
      accp_omega_c=0
      accp_eta_c=0
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
          chain_v_beta_square=mcmc(v_beta_square_record)
          chain_v_beta_square_c=mcmc(v_beta_square_record_c)
          chain_v_omega_square=mcmc(v_omega_square_record)
          chain_v_omega_square_c=mcmc(v_omega_square_record_c)
          combinedchains_beta = mcmc.list(chain_beta, chain_beta_c)
          combinedchains_omega = mcmc.list(chain_omega, chain_omega_c)
          combinedchains_tau_square = mcmc.list(chain_tau_square, chain_tau_square_c)
          combinedchains_eta = mcmc.list(chain_eta, chain_eta_c)
          combinedchains_v_beta_square = mcmc.list(chain_v_beta_square, chain_v_beta_square_c)
          combinedchains_v_omega_square = mcmc.list(chain_v_omega_square, chain_v_omega_square_c)
          if (gelman.diag(combinedchains_beta)$psrf[1,1]<1.1 & gelman.diag(combinedchains_omega)$mpsrf<1.1 & gelman.diag(combinedchains_tau_square)$psrf[1,1]<1.1 & gelman.diag(combinedchains_eta)$psrf[1,1]<1.1 & gelman.diag(combinedchains_v_beta_square)$psrf[1,1]<1.1 & gelman.diag(combinedchains_v_omega_square)$psrf[1,1]<1.1){
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
  
  
  
  
  
  
  Gauss_bayes=list(beta=beta,omega=omega,tau_square=tau_square,eta=eta,v_beta_square=v_beta_square,v_omega_square=v_omega_square,K_inv=K_inv,K=K,eps_beta=eps_beta,eps_omega=eps_omega,beta_record=beta_record,omega_record=omega_record,tau_square_record=tau_square_record,eta_record=eta_record,v_beta_square_record=v_beta_square_record,v_omega_square_record=v_omega_square_record, eps_beta_record=eps_beta_record,eps_omega_record=eps_omega_record,beta_c=beta_c,omega_c=omega_c,tau_square_c=tau_square_c,eta_c=eta_c,v_beta_square_c=v_beta_square_c,v_omega_square_c=v_omega_square_c,K_inv_c=K_inv_c,K_c=K_c,eps_beta_c=eps_beta_c,eps_omega_c=eps_omega_c,beta_record_c=beta_record_c,omega_record_c=omega_record_c,tau_square_record_c=tau_square_record_c,eta_record_c=eta_record_c,v_beta_square_record_c=v_beta_square_record_c,v_omega_square_record_c=v_omega_square_record_c, eps_beta_record_c=eps_beta_record_c,eps_omega_record_c=eps_omega_record_c,convergence=convergence)
  class(Gauss_bayes)="Gaussian_bayes"
  return(Gauss_bayes)
  
  
  
  
  
  
  
  
  
  
  
  
  
}