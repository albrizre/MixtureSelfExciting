simulate_self_exciting_mixture <- function(alpha,offspring_force,lambda_exp,Tmax,under_mech){
  # alpha -> background intensity of the process
  # Tmax -> length of the temporal window
  # under_mech -> vector with the underlying latent states (in the form of a condition)
  
  # Simulate events from the background process
  inter = rexp(1000, alpha) # 1000 should be large enough for alpha
  back_times = cumsum(inter)
  back_times = back_times[back_times<Tmax]
  # length(back_times)/Tmax this should be close to alpha
  
  # Simulate offspring
  l = 0
  condition_while=T
  all_times=back_times
  gen_times=back_times[eval(parse(text=under_mech))] # the background times are the times for the starting generation
  while (condition_while){
    l=l+1
    # print(l)
    # offspring_force = rgamma(length(gen_times),1,2)
    offspring_n = rpois(length(gen_times),offspring_force)
    offspring_times=c()
    for (i in 1:length(offspring_n)){
      if (offspring_n[i]>0){
        offspring_times = c(offspring_times, gen_times[i]+rexp(offspring_n[i],lambda_exp))
      }
    }
    offspring_times=offspring_times[offspring_times<Tmax]
    all_times=sort(c(all_times,offspring_times))
    if (length(offspring_times)>0){
      gen_times=offspring_times # set new times for the current generation
    } else {
      condition_while<-FALSE
    }
  }
  output=list()
  output[["back_times"]]=back_times
  output[["all_times"]]=all_times
  return(output)
}
