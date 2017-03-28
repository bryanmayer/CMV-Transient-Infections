#Parameters for stochastic runs
#Beta and initial values (I or V), set in the analysis

#load("~/Dropbox/CMV-Primary-Analysis/Mathematical Model/Initial dynamics/initial_values.RData")
#betas from no latency = c(2.535e-12, 3.315e-12, 5.087e-12 ) #min, median, max
#summary(base_model_fits$beta)
betas = c(2.535e-12, 3.315e-12, 5.087e-12) #min, median, max

#parms fixed from previous CMV primary infection
parms = data.frame(
  beta = NA,
  c = 2,
  K = 1e7 * 40,  #Dawes
  mu = 1/4.5, #Dawes  is 1/4.5
  lambda = 4e8 / 4.5, #lambda = mu * S0
  delta = 0.77, #assigned directly in the batch file using prop parameter
  p = 1600 #from PNAS temperature paper #old 60
)

#load("~/Dropbox/CMV-Primary-Analysis/Mathematical Model/Latency model/Initial dynamics/initial_latent_values.RData")
#summary(base_model_fits$beta)
betas_latent = c(3.192e-12, 4.800e-12, 9.134e-12) #min, median, max

#parms fixed from previous CMV primary infection
parms_latent = data.frame(
  beta = NA,
  c = 2,
  K = 1e7 * 40,  #Dawes
  mu = 1/4.5, #Dawes  is 1/4.5
  lambda = 4e8 / 4.5, #lambda = mu * S0
  alpha = 1,
  delta = 0.77, #assigned directly in the batch file using prop parameter
  p = 1600 #from PNAS temperature paper #old 60
)
