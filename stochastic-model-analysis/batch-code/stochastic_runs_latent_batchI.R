args<-(commandArgs(TRUE));
if(length(args)==0){
  print("No arguments supplied.")

}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
  print(args)
}

library(plyr)
library(dplyr, lib = "/home/bmayer/CMVGrowth/libraryFolder")
library(doParallel)
registerDoParallel(cores = 8)
source("stochastic_model.R")


#baseparmsIn = read.csv("model_fit_growth_informed.csv")
starttime <- Sys.time()

#load("~/Dropbox/CMV-Primary-Analysis/Mathematical Model/Initial dynamics/initial_values.RData")
#betas from no latency = c(2.535e-12, 3.315e-12, 5.087e-12 ) #min, median, max

I0s = c(1:10, 15, 10 * 2:9, 1:4 %o% 10^2)

parms = data.frame(
  beta = NA,
  c = 2,
  K = 1e7 * 40,  #Dawes
  alpha = 1,
  mu = 1/4.5, #Dawes  is 1/4.5
  lambda = 4e8 / 4.5, #lambda = mu * S0
  delta = 0.77, #assigned directly in the batch file using prop parameter
  p = 1600 #from PNAS temperature paper
)

if(beta_i == 1) betas = c(1.01, 1.05, 1.1)/with(parms, p * K/ (c * delta *(1+mu/alpha)))
if(beta_i == 2) betas = seq(1.15, 1.3, 0.05)/with(parms, p * K/ (c * delta *(1+mu/alpha)))
if(beta_i == 3) betas = seq(1.35, 1.65, 0.05)/with(parms, p * K/ (c * delta *(1+mu/alpha)))
if(beta_i == 4) betas = seq(1.7, 2.15, 0.05)/with(parms, p * K/ (c * delta *(1+mu/alpha)))

stochastic_sim_batch = ldply(1:length(betas), function(b){
  out = ldply(1:length(I0s), function(I0){

    ldply(1:10000, function(i){
      parms$beta = betas[b]
      model_sim = stochastic_model_latent_I(500, parms, initI0 = 0, initI = I0s[I0])

      #observed blips, just assuming it's consecutive
      observable_index = which(model_sim$V >= 150)
      if(length(observable_index) > 0){
        blip_obs_time = with(model_sim, time[max(observable_index)] - time[min(observable_index)])
        if(length(observable_index) == 1) random_sample =  model_sim$V[observable_index] else random_sample = sample(x = model_sim$V[observable_index], size = 1)
        sample_time_1 = 7*runif(1)
        if(sample_time_1 <= max(blip_obs_time)){
          consecutive_blips = 1 + floor((max(blip_obs_time) - sample_time_1)/7)

        } else  consecutive_blips = 0

      } else{
        blip_obs_time = NA
        random_sample = NA
        consecutive_blips = 0
      }

      data.frame(
        run = i,
        c = 2,
        K = 1e7 * 40,  #Dawes
        mu = 1/4.5, #Dawes  is 1/4.5
        lambda = 4e8 / 4.5, #lambda = mu * S0
        alpha = 1,
        delta = 0.77, #assigned directly in the batch file using prop parameter
        p = 1600, #from PNAS temperature paper
        beta = parms$beta,
        initI = I0s[I0],
        lagIinit = model_sim$I0[1],
        repIinit = model_sim$I[1],
        R0 = with(parms, beta * K * p/(c * delta *(1+mu/alpha))),
        max_time = max(model_sim$time),
        last_V = tail(model_sim$V, 1),
        max_V = max(model_sim$V),
        max_I = max(model_sim$I),
        observed = (length(observable_index) > 0),
        blip_obs_time = blip_obs_time,
        random_blip_sample = random_sample,
        consecutive_blips = consecutive_blips
      )
    }, .parallel = T)
  })
})

file_out = paste("stochastic_sim_batchlatentI", beta_i, ".RData", sep = "")

save(stochastic_sim_batch, file = file_out)

Sys.time() - starttime

