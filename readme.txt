###
header.R: Include different R packages and R files we need to use.

adp_HMC.R: Implement Algorithm 1 (Hamiltonian Monte Carlo) with adaptive step size tuning.

adp_SphHMC.R: Implement Algorithm 2 (Spherical Hamiltonian Monte Carlo) with adaptive step size tuning.

bayesian_bridge_adp_convergence.R: Implement Algorithm 3 (sampling procedure for sph-B^2GPR) for models with a non-constant fixed-effect term (vector-valued).

bayesian_bridge_adp_convergence2.R: Implement Algorithm 3 (sampling procedure for sph-B^2GPR) for models with a constant fixed-effect term (scalar).

bayesian_bridge_adp_HMC_convergence.R: Implement Algorithm 4 (sampling procedure for hmc-B^2GPR) for models with a non-constant fixed-effect term (vector-valued).

bayesian_bridge_adp_HMC_convergence2.R: Implement Algorithm 4 (sampling procedure for hmc-B^2GPR) for models with a constant fixed-effect term (scalar).


dual_avg.R: Implement the dual-averaging scheme used to adaptively tune the step size during the burn-in phase.

init_eps_HMC.R: Choose the initial step size for Algorithm 1 (Hamiltonian Monte Carlo).

init_eps.R: Choose the initial step size for Algorithm 2 (Spherical Hamiltonian Monte Carlo).

leapfrog_HMC.R: Implement the leapfrog integration (for-loop part) used in Algorithm 1 (Hamiltonian Monte Carlo).

leapfrog.R: Implement the leapfrog integration (for-loop part) used in Algorithm 2 (Spherical Hamiltonian Monte Carlo).

metropolis.R: Implement the Metropolis-Hastings algorithm.

predict.R: Use the trained Bayesian bridge model to generate predictions.

predict_c.R: Use an alternative trained Bayesian bridge model to generate predictions for comparison.

second_order.R: Construct the second-order polynomial basis.

###

