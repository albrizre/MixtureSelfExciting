library(nimble)
source("NIMBLE_codes/mixture_mech_empir_changepoints_Kchanges.R")
  
alphas=c(0.2) # background intensity of the process
alpha=alphas[1]
Tmax=200 # Length of the temporal window
offspring_force = 0.8 # Average number of offspring events
lambda_exp = 1

K = 2 # Choose number of change points

for (sim in 1:100){

  load(paste0("Simulated_datasets/Sim_",sim,"_alpha_",alpha,"_force_",offspring_force,"_lambda_exp_",lambda_exp,"_B.rda"))
  
  # Sort by Time

  all_times = sort(all_times)

  # Times of the events and partition_changes of the time interval

  time_events=all_times
  time_points=seq(0,Tmax,1)

  # Model construction and call

  n_previous_aux=c()
  for (i in 1:length(time_points)){
    if (length(which(time_events<time_points[i]))==0){
      n_previous_aux=c(n_previous_aux,0)
    } else {
      n_previous_aux=c(n_previous_aux,which(time_events<time_points[i])[length(which(time_events<time_points[i]))])
    }
  }

  constants <- list(
  
    N = length(all_times)+length(time_points),
    N_events = length(all_times),
    N_integral = length(time_points),
    w=c(rep(0,length(all_times)),rep(1,length(time_points))), # 1 is the length of the time interval (in days)
    I=c(rep(1,length(all_times)),rep(0,length(time_points))),
    n_previous=c(0:(length(all_times)-1),n_previous_aux),
    time_events = time_events,
    time_points = time_points,
    time_all = c(time_events,time_points),
    C=1000000,
    time_point = floor(c(time_events,time_points))+1,
    lambda0_base=length(all_times)/Tmax,
    
    Tmax = Tmax + 1,
    
    probs_cat = rep(0.01,Tmax),
    vector_indices = 1:(K+1),
    
    K=K
  
  )

  matrix_ones=matrix(0,nrow=constants$N,constants$N_events)
  for (i in 1:constants$N){
    if (constants$n_previous[i]>=1){
      matrix_ones[i,1:constants$n_previous[i]]=1 
    }
  }
  constants$matrix_ones=matrix_ones
  
  matrix_distances=matrix(0,nrow=constants$N,constants$N_events)
  for (i in 1:constants$N){
    if (constants$n_previous[i]>=1){
      matrix_distances[i,1:constants$n_previous[i]]=constants$time_all[i]-constants$time_events[1:constants$n_previous[i]]
    }
  }
  constants$matrix_distances=matrix_distances
  
  # Fit model
  
  data <- list(zeros = rep(0,constants$N))
  pars_monitor=c("lambda0",
               "theta",
               "phi",
               "Mech",
               "Mech_i",
               "lambda",
               "change_point",
               "partition_changes",
               "interval_aux",
               "aux",
               "L")
  model_name="Kchanges"
  # change_point = sort(runif(constants$K))*Tmax
  change_point = sort(runif(constants$K,0,constants$Tmax))
  partition_changes = c(0,change_point,Tmax)
  inits <- function() list(lambda0 = constants$lambda0_base,
                         theta = 0,
                         phi = 50,
                         Mech = c(0,1,0),
                         Mech_const = c(0,1,0),
                         Mech_i = rep(1,constants$N),
                         change_point = change_point,
                         partition_changes = c(0,change_point,constants$Tmax),
                         incrc = rep(0,constants$K),
                         pi = rep(0.5,constants$K+1))

  if (!file.exists(paste0("Simulated_Models/Sim_",sim,"_alpha_",alpha,"_force_",offspring_force,"_lambda_exp_",lambda_exp,"_",model_name,"_",K,"_B.rda"))){
      mcmc.output <- nimbleMCMC(mixture_mech_empir_changepoints_Kchanges, data = data, inits = inits, constants = constants,
                          monitors = pars_monitor, 
                          niter = 200000, nburnin = 100000, nchains = 1, thin = 10,
                          summary = TRUE, WAIC = TRUE)
      save(mcmc.output,file=paste0("Simulated_Models/Sim_",sim,"_alpha_",alpha,"_force_",offspring_force,"_lambda_exp_",lambda_exp,"_",model_name,"_",K,"_B.rda"))
    }


}
