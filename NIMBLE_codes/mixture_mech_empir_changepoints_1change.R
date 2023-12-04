mixture_mech_empir_changepoints_1change <- nimbleCode({
  
  lambda0 ~ dgamma(10*lambda0_base,10)
  lambda0_Mech ~ dgamma(10*lambda0_base,10)
  theta ~ dunif(0,20)
  phi ~ dunif(0,50)
  # phi ~ dgamma(7,1)
  
  for (i in 1:N) {
    
    Mech_i[i] <- Mech[step(time_point[i] - change_point) + 1]
    lambda[i] <- lambda0+Mech_i[i]*(theta*inprod(matrix_ones[i,],exp(-matrix_distances[i,]/phi)))
    log(L[i]) <- I[i]*log(lambda[i])-w[i]*lambda[i] # I_i = 1 if event, 0 otherwise; w[i] = 'interval length' if partition point, 0 otherwise
    z[i] <- -log(L[i])+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  change_point ~ dunif(0,Tmax)
  
  partition_changes[1] <- 0 # if we assume K change points, then we have K+2 points in the partition
  for (k in 1:K){
    partition_changes[k+1] <- change_point
  }
  partition_changes[K+2] <- Tmax
  
  Mech[1] ~ dbern(0.5)
  for (k in 2:(K+1)){ # if we assume K change points, then we have K+1 intervals defined
    Mech[k] <- 1 - Mech[k-1]
  }
  
  
})