mixture_mech_empir_poisCP <- nimbleCode({
  
  lambda0 ~ dgamma(10*lambda0_base,10)
  theta ~ dunif(0,20)
  phi ~ dunif(0,50)
  # phi ~ dgamma(7,1)
  
  for (i in 1:N) {
    
    for (k in 1:(K+1)){ 
      interval_aux[i,k] <- step(time_point[i] - partition_changes[k])*step(partition_changes[k+1] - time_point[i])
    }
    aux[i] <- inprod(interval_aux[i,1:(K+1)],vector_indices[1:(K+1)])
    lambda[i] <- lambda0+lambdac[aux[i]]
    log(L[i]) <- I[i]*log(lambda[i])-w[i]*lambda[i] # I_i = 1 if event, 0 otherwise; w[i] = 'interval length' if partition point, 0 otherwise
    z[i] <- -log(L[i])+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in 1:(K+1)){
    lambdac[i] ~ dgamma(0.01,0.01)
  }
  
  change_point[1] ~ dunif(0,Tmax)
  incrc[1] <- 0
  for (k in 2:K){
    incrc[k] ~ dunif(0,Tmax - change_point[k-1])
    change_point[k] <- change_point[k-1] + incrc[k]
  }
  
  partition_changes[1] <- 0 # if we assume K change points, then we have K+2 points in the partition
  for (k in 1:K){
    partition_changes[k+1] <- change_point[k]
  }
  partition_changes[K+2] <- Tmax
  
#  Mech[1] ~ dbern(0.5)
#  for (k in 2:(K+1)){ # if we assume K change points, then we have K+1 intervals defined
#    Mech[k] <- 1 - Mech[k-1]
#  }
  
  
  
  
})