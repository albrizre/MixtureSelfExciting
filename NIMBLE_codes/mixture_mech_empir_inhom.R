mixture_mech_empir_inhom <- nimbleCode({
  
  lambda0 ~ dgamma(10*lambda0_base,10)
  
  gamma1 ~ dnorm(0,0.1)
  gamma2 ~ dnorm(0,0.1)
  
  for (i in 1:N) {
    
    lambda[i] <- exp(lambda0+gamma1*sin((2*pi*time_all[i])/30)+gamma2*cos((2*pi*time_all[i])/30)+delta[time_point[i]])
    log(L[i]) <- I[i]*log(lambda[i])-w[i]*lambda[i] # I_i = 1 if event, 0 otherwise; w[i] = 'interval length' if partition point, 0 otherwise
    z[i] <- -log(L[i])+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  # RW2 prior on the effect of days
  delta[1:N_integral] ~ dcar_normal(adj[1:N_W_adj], weights[1:N_W_adj], num[1:N_integral], tau.delta, zero_mean = 1)
  sigma2.delta ~ dgamma(1,1)
  tau.delta <- 1/sigma2.delta
  
})