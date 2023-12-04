library(nimble)
library(sp)
library(rgdal)
library(rgeos)
library(sf)
library(spdep)
library(maptools)
library(spatstat)
source("Functions/simulate_self_exciting.R")
source("Functions/simulate_self_exciting_mixture.R")

alphas=c(0.2) # background intensity of the process
alpha=alphas[1]
Tmax=200 # Length of the temporal window
offspring_force = 0.8 # Average number of offspring events
lambda_exp = 1
under_mech="back_times < 40 | back_times > 160"
under_mech="back_times < 40"

if (under_mech=="back_times < 40"){
  scenario="A"
}
if (under_mech=="back_times < 40 | back_times > 160"){
  scenario="B"
}
set.seed(123)
for (alpha in alphas){
  
  for (i in 1:100){
    
    sim_pattern=simulate_self_exciting_mixture(alpha,offspring_force,lambda_exp,Tmax,under_mech)
    back_times=sim_pattern$back_times
    all_times=sim_pattern$all_times
    print(paste0(i,", n = ",length(all_times)))
    one = rep(1, times = length(back_times))
    plot(back_times, one, yaxt = 'n', ann=FALSE, pch=19)
    triggered_times=all_times[!all_times%in%back_times]
    one = rep(1, times = length(triggered_times))
    points(triggered_times, one, pch=19, col="red")
    plot(all_times)
    
    save(all_times,file=paste0("Simulated_datasets/Sim_",i,"_alpha_",alpha,"_force_",offspring_force,"_lambda_exp_",lambda_exp,"_",scenario,".rda"))
    
  }
  
}

