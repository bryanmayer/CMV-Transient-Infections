#!/bin/bash

for x in {1,2,3,4}; do
sbatch --cpus-per-task=8 --time=0-24 --wrap="R --no-save --no-restore '--args beta_i=$x' < stochastic_runs_latent_batchI.R"
done