mixture_mech_empir_basic <- nimbleCode({
  
  lambda0 ~ dgamma(10*lambda0_base,10)
  theta ~ dunif(0,20)
  phi ~ dunif(0,50)
  
  for (i in 1:N) {
    
    lambda[i] <- lambda0+theta*inprod(matrix_ones[i,],exp(-matrix_distances[i,]/phi))
    log(L[i]) <- I[i]*log(lambda[i])-w[i]*lambda[i] # I_i = 1 if event, 0 otherwise; w[i] = 'interval length' if partition point, 0 otherwise
    z[i] <- -log(L[i])+C
    zeros[i] ~ dpois(z[i])
    
  }
  
})
