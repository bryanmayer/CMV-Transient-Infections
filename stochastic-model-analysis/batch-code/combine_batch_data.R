library(plyr)

file_name = "stochastic_sim_batchlatentI"

stochastic_sim_latency_I = ldply(1:4, function(beta_i){
  load(paste(file_name, beta_i, ".RData", sep =""))
  stochastic_sim_batch
})

save(stochastic_sim_latency_I, file = paste("stochastic_sim_I_latent.Rdata", sep = ""))
