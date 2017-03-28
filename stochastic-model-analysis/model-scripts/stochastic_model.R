#these models used fixed S
CMVModel_ode = function(t, x, parms, AUCData = NULL){
  with(as.list(c(parms, x)), {
    dS <- lambda - mu * S - beta * S * V
    dI <- beta * S * V  -  delta * I
    dV <- p * I - c * V 
    
    res <- c(dS, dI, dV)
    list(res)
  })
}

CMVModel_latent_ode = function(t, x, parms){
  with(as.list(c(parms, x)), {
    dS <- lambda - mu * S - beta * S * V
    dI0 <- beta * S * V - mu * I0 - alpha * I0
    dI <- alpha * I0 -  delta * I
    dV <- p * I - c * V 
    
    res <- c(dS, dI0, dI, dV)
    list(res)
  })
}

stochastic_model_I = function(max_time, parms, tao = 0.01, initI = 1, seed_set = NULL){
  if(!is.null(seed_set)) set.seed(seed_set) #reproduce results
  
  time_vec = seq(0, max_time, tao)
  out_index = length(time_vec)
  
  
  I_vec = rep(0, out_index)
  V_vec = rep(0, out_index)
  
  #initial
  S = parms$K
  V = 0
  I = initI
  
  I_vec[1] = initI
  for(i in 2:out_index){
    
  events = tao * with(as.list(parms),
                        c(
                          beta * S * V,
                          delta * I, 
                          p * I, 
                          c * V
                        )
    )
    
  if(length(which(events < 0)) > 0) events[which(events < 0)] = 0 #no negative events
  
  infected = rpois(1, events[1])
  infected_death = rpois(1, events[2])
  newV = rpois(1, events[3])
  deathV = rpois(1, events[4]) #death and immune clr    
  
  #if(i > 100) browser()
  
  I = max(0, I + infected - infected_death)
  V = max(0, V + newV - deathV)
  
  #these are updated (or not updated) in the conditionals
  I_vec[i] = I
  V_vec[i] = V
  
  if(is.na(V)) browser()
  
  if((I == 0 & V == 0) | V > 1e9){
    out_index = i
    break
  }
  }
  data.frame(
    time = time_vec[1:out_index],
    I = I_vec[1:out_index],
    V = V_vec[1:out_index]
  )
}

stochastic_model_latent_I = function(max_time, parms, tao = 0.1, initI0 = 1, initI = 0, seed_set = NULL){
  if(!is.null(seed_set)) set.seed(seed_set) #reproduce results
  
  time_vec = seq(0, max_time, tao)
  out_index = length(time_vec)
  
  I0_vec =rep(0, out_index)
  I_vec = rep(0, out_index)
  V_vec = rep(0, out_index)
  
  #initial
  S = parms$K
  V = 0
  I0 = initI0
  I = initI
  
  I_vec[1] = initI
  I0_vec[1] = initI0
  
  for(i in 2:out_index){
    
    events = tao * with(as.list(parms),
                        c(
                          beta * S * V, 
                          alpha * I0,
                          mu * I0,
                          delta * I, 
                          p * I, 
                          c * V
                        )
    )
    
    if(length(which(events < 0)) > 0) events[which(events < 0)] = 0 #no negative events
    
    infected = rpois(1, events[1])
    latent_out = rpois(1, events[2])
    latent_death = rpois(1, events[3])
    infected_death = rpois(1, events[4])
    newV = rpois(1, events[5])
    deathV = rpois(1, events[6]) #death and immune clr    
    
    #if(i > 100) browser()
    
    I0 = max(0, I0 + infected - latent_out - latent_death)
    I = max(0, I + latent_out - infected_death)
    V = max(0, V + newV - deathV)
    
    #these are updated (or not updated) in the conditionals
    I0_vec[i] = I0
    I_vec[i] = I
    V_vec[i] = V
    
    if(is.na(V)) browser()
    
    if((I == 0 & V == 0 & I0 == 0) | V > 1e9){
      out_index = i
      break
    }
  }
  data.frame(
    time = time_vec[1:out_index],
    I0 = I0_vec[1:out_index],
    I = I_vec[1:out_index],
    V = V_vec[1:out_index]
  )
}


stochastic_model = function(max_time, parms, tao = 0.01, initV = 10, seed_set = NULL,
                            exposure_rate = 0, #per day
                            exposure_mean = 0, #per exposure, log10 scale,
                            exposure_variance = 0){ #zero just means constant)
  
  if(!is.null(seed_set)) set.seed(seed_set) #reproduce results
  
  time_vec = seq(0, max_time, tao)
  out_index = length(time_vec)
  
  I_vec = rep(0, out_index)
  V_vec = rep(0, out_index)
  
  #initial
  S = parms$K
  V = initV
  I = 0
  
  V_vec[1] = initV
  
  for(i in 2:out_index){
    
    events = tao * with(as.list(parms),
                        c(
                          beta * V, 
                          delta * I, 
                          p * I, 
                          c * V
                        )
    )
    if(length(which(events < 0)) > 0) events[which(events < 0)] = 0 #no negative events
    
    infected = rbinom(1, S, events[1])
    infected_death = rpois(1, events[2])
    newV = rpois(1, events[3])
    deathV = rpois(1, events[4])     
    
    #exposure of virus
    if(runif(1) < exposure_rate * tao) expV = rlnorm(1, exposure_mean * log(10), exposure_variance)  else expV = 0
    
    I = max(0, I + infected - infected_death)
    V = max(0, V + newV - deathV + expV)
    
    I_vec[i] = I
    V_vec[i] = V
    
    if(is.na(V)) browser()
    
    #just a single exposure analysis then stop early if I and V run out
    #Or R breaks when V > 1e10
    if(((exposure_mean == 0 | exposure_rate == 0) & (I == 0 & V == 0)) | V > 1e9){
      out_index = i
      break
    }
    
  }
  data.frame(
    time = time_vec[1:out_index],
    I = I_vec[1:out_index],
    V = V_vec[1:out_index]
  )
}

#takes stochastic simulation in and samples virus
data_sampling = function(data_in, 
                         sampling_space = 7,#days between sample
                         noise_var = 0, #noise is all variance, noise_mean != 0 means bias
                         noise_mean = 0){
  
  #sampling days
  sample_index = which(data_in$time %% sampling_space == 0)
  
  noise = rnorm(length(sample_index), 10^noise_mean, 10^noise_var)
  sampled_V = data_in$V[sample_index] 
  sampled_V_noise = sampled_V + noise
  
  data.frame(
    time = data_in$time[sample_index],
    count = ifelse(sampled_V_noise < 150, 0, log10(sampled_V_noise)),
    real_count = ifelse(sampled_V < 150, 0, log10(sampled_V))
  )
} 



