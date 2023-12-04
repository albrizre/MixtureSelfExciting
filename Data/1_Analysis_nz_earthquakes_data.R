library(nimble)
source("NIMBLE_codes/mixture_mech_empir_changepoints_Kchanges.R")
source("NIMBLE_codes/mixture_mech_empir_changepoints_1change.R")
source("NIMBLE_codes/mixture_mech_empir_basic.R")
source("NIMBLE_codes/mixture_mech_empir_poisCP.R")
source("NIMBLE_codes/mixture_mech_empir_poisCP_1change.R")
source("NIMBLE_codes/mixture_mech_empir_inhom.R")
source("Functions/RW2.R")
load("Data/nz_earthquakes.rda")

nz_earthquakes$origintime=gsub("UTC","",nz_earthquakes$origintime)
nz_earthquakes$date=as.numeric(as.Date(substr(nz_earthquakes$origintime,1,10)))
nz_earthquakes$date=nz_earthquakes$date-min(nz_earthquakes$date)+1
for (i in 1:nrow(nz_earthquakes)){
  if (nchar(nz_earthquakes$origintime[i]==10)){
    nz_earthquakes$origintime[i]=paste0(nz_earthquakes$origintime[i]," 00:00:00")
  }
}
nz_earthquakes$Time=nz_earthquakes$date+as.numeric(substr(nz_earthquakes$origintime,12,13))/24+as.numeric(substr(nz_earthquakes$origintime,15,16))/1440
#nz_earthquakes=nz_earthquakes[nz_earthquakes$Time>=200 & nz_earthquakes$Time<=600,]

# Sort by Time
nz_earthquakes=nz_earthquakes[order(nz_earthquakes$Time),]
all_times=nz_earthquakes$Time
all_times=all_times[3000:length(all_times)]
all_times=all_times-590

# Times of the events and partition_changes of the time interval

Tmax=max(all_times)
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

########################
# Mixture K (>1) changes
########################

for (K in 2:5){

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
  
  set.seed(12345)
  data <- list(zeros = rep(0,constants$N))
  change_point = sort(runif(constants$K,0,constants$Tmax))
  partition_changes = c(0,change_point,Tmax)
  inits <- function() list(lambda0 = constants$lambda0_base,
                           theta = 0,
                           phi = 50,
                           Mech = c(0,1,0,1,0,1)[1:(K+1)],
                           Mech_const = c(0,1,0,1,0,1)[1:(K+1)],
                           Mech_i = rep(1,constants$N),
                           change_point = change_point,
                           partition_changes = c(0,change_point,constants$Tmax),
                           incrc = rep(0,constants$K),
                           pi = rep(0.5,constants$K+1))
                           
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
  
  for (model_name in c("Kchanges")){
      if (!file.exists(paste0("Models/nz_earthquakes_",model_name,"_",K,".rda"))){
      mcmc.output <- nimbleMCMC(mixture_mech_empir_changepoints_Kchanges, data = data, inits = inits, constants = constants,
                              monitors = pars_monitor, thin = 10,
                              niter = 200000, nburnin = 100000, nchains = 1, 
                              summary = TRUE, WAIC = TRUE)
      save(mcmc.output,file=paste0("Models/nz_earthquakes_",model_name,"_",K,".rda"))
    }
  }
}

##################
# Mixture 1 change
##################

for (K in c(1)){

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
  
  set.seed(12345)
  data <- list(zeros = rep(0,constants$N))
  change_point = sort(runif(constants$K,0,constants$Tmax))
  partition_changes = c(0,change_point,Tmax)
  inits <- function() list(lambda0 = constants$lambda0_base,
                           theta = 0,
                           phi = 50,
                           Mech = c(0,1),
                           Mech_const = c(0,1),
                           Mech_i = rep(1,constants$N),
                           change_point = change_point,
                           partition_changes = c(0,change_point,constants$Tmax),
                           incrc = rep(0,constants$K),
                           pi = rep(0.5,constants$K+1))
                           
  pars_monitor=c("lambda0",
                 "theta",
                 "phi",
                 "Mech",
                 "Mech_i",
                 "lambda",
                 "change_point",
                 "partition_changes")
  
  for (model_name in c("1change")){
      if (!file.exists(paste0("Models/nz_earthquakes_",model_name,"_",K,".rda"))){
      mcmc.output <- nimbleMCMC(mixture_mech_empir_changepoints_1change, data = data, inits = inits, constants = constants,
                              monitors = pars_monitor, thin = 10,
                              niter = 200000, nburnin = 100000, nchains = 1, 
                              summary = TRUE, WAIC = TRUE)
      save(mcmc.output,file=paste0("Models/nz_earthquakes_",model_name,"_",K,".rda"))
    }
  }
}

#####################
# Self-exciting model
#####################

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
  
  set.seed(12345)
  data <- list(zeros = rep(0,constants$N))
  change_point = sort(runif(constants$K,0,constants$Tmax))
  partition_changes = c(0,change_point,Tmax)
  inits <- function() list(lambda0 = constants$lambda0_base,
                           theta = 0,
                           phi = 50)
                           
  pars_monitor=c("lambda0",
                 "theta",
                 "phi",
                 "lambda")
  
  for (model_name in c("self")){
      if (!file.exists(paste0("Models/nz_earthquakes_",model_name,".rda"))){
      mcmc.output <- nimbleMCMC(mixture_mech_empir_basic, data = data, inits = inits, constants = constants,
                              monitors = pars_monitor, thin = 10,
                              niter = 200000, nburnin = 100000, nchains = 1, 
                              summary = TRUE, WAIC = TRUE)
      save(mcmc.output,file=paste0("Models/nz_earthquakes_",model_name,".rda"))
    }
  }
  
####################################
# Poisson homogeneous 1 change point
####################################

for (K in 1:1){
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
    
  set.seed(12345)
  data <- list(zeros = rep(0,constants$N))
  change_point = sort(runif(constants$K,0,constants$Tmax))
  partition_changes = c(0,change_point,Tmax)
  inits <- function() list(lambda0 = constants$lambda0_base,
                             theta = 0,
                             phi = 50,
                             lambdac = rep(1,constants$K + 1),
                             change_point = change_point,
                             partition_changes = c(0,change_point,constants$Tmax),
                             incrc = rep(0,constants$K),
                             pi = rep(0.5,constants$K+1))
                             
  pars_monitor=c("lambda0",
                   "theta",
                   "phi",
                   "lambda",
                   "lambdac",
                   "change_point",
                   "partition_changes")
    
  for (model_name in c("poisCP")){
    if (!file.exists(paste0("Models/nz_earthquakes_",model_name,"_",K,".rda"))){
    mcmc.output <- nimbleMCMC(mixture_mech_empir_poisCP_1change, data = data, inits = inits, constants = constants,
                                monitors = pars_monitor, thin = 10,
                                niter = 200000, nburnin = 100000, nchains = 1, 
                                summary = TRUE, WAIC = TRUE)
    save(mcmc.output,file=paste0("Models/nz_earthquakes_",model_name,"_",K,".rda"))
    }
  }
}

##########################################
# Poisson homogeneous K (>1) change points
##########################################

for (K in 2:5){
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
    
  set.seed(12345)
  data <- list(zeros = rep(0,constants$N))
  change_point = sort(runif(constants$K,0,constants$Tmax))
  partition_changes = c(0,change_point,Tmax)
  inits <- function() list(lambda0 = constants$lambda0_base,
                             theta = 0,
                             phi = 50,
                             lambdac = rep(1,constants$K + 1),
                             change_point = change_point,
                             partition_changes = c(0,change_point,constants$Tmax),
                             incrc = rep(0,constants$K),
                             pi = rep(0.5,constants$K+1))
                             
  pars_monitor=c("lambda0",
                   "theta",
                   "phi",
                   "lambda",
                   "lambdac",
                   "change_point",
                   "partition_changes")
    
  for (model_name in c("poisCP")){
    if (!file.exists(paste0("Models/nz_earthquakes_",model_name,"_",K,".rda"))){
    mcmc.output <- nimbleMCMC(mixture_mech_empir_poisCP, data = data, inits = inits, constants = constants,
                                monitors = pars_monitor, thin = 10,
                                niter = 200000, nburnin = 100000, nchains = 1, 
                                summary = TRUE, WAIC = TRUE)
    save(mcmc.output,file=paste0("Models/nz_earthquakes_",model_name,"_",K,".rda"))
    }
  }
}

#######################
# Poisson inhomogeneous 
#######################

# Week
W <- length(time_points)
RW2_W <- RW2(W)

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
      
      N_W = length(RW2_W$num),
      N_W_adj = length(RW2_W$adj),
  
      adj = RW2_W$adj,                             
      num = RW2_W$num,
      weights = RW2_W$weights,
      
      pi = pi 
      
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
    
  set.seed(12345)
  data <- list(zeros = rep(0,constants$N))
  inits <- function() list(lambda0 = constants$lambda0_base,
                             gamma1 = 0, gamma2 = 0, delta = rep(0,constants$N_integral), sigma2.delta = 1)
                             
  pars_monitor=c("lambda0",
                  "lambda","gamma1","gamma2","delta","sigma2.delta")
    
  for (model_name in c("inhom")){
    if (!file.exists(paste0("Models/nz_earthquakes_",model_name,".rda"))){
    mcmc.output <- nimbleMCMC(mixture_mech_empir_inhom, data = data, inits = inits, constants = constants,
                                monitors = pars_monitor, thin = 10,
                                niter = 500000, nburnin = 250000, nchains = 1, 
                                summary = TRUE, WAIC = TRUE)
    save(mcmc.output,file=paste0("Models/nz_earthquakes_",model_name,".rda"))
    }
  }