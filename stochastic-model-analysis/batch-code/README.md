These are the scripts I used to run jobs on the server. 

1. stochastic_runs_latent_batchI.R - runs the simulation for a series of input parameters depending on the shells script call.

2. mu_fit_latency.sh - the shell script I used.

3. combine_batch_data.R - a short script that I used to combine the data into a single .RData file.

stochastic_model.R is available in stochastic-model-analysis/

These scripts were run on the scicomp cluster at the Fred Hutch and are being shared for informational purposes. The final data set is large (115mb) and is uploaded to dataverse (stochastic_sim_I_latent.Rdata. TBD: add citation).