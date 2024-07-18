# Author:         Lucas Jeay-Bizot
# Created:        July 18th 2024
# Last Modified:  July 18th 2024

# Purpose: This code, using the BFDA package, will find the smallest n_max such
# that, using a Bayesian Sequential sampling design, for a true effect size of
# 0.5 (moderate) more than 80% of the studies would return strong evidence and
# more than 95% of the studies would return moderate evidence for the effect. It
# is done assuming data collected using a stopping threshold set for strong
# evidence.

# References: Schönbrodt, F. D. & Stefan, A. M. (2019). BFDA: An R package for
# Bayes factor design analysis (version 0.5.0). Retrieved from
# https://github.com/nicebread/BFDA

# Dependencies
library(BFDA) # Package for Bayes Factor Design Analysis (https://github.com/nicebread/BFDA)

## INPUT ## Effect sizes (rounded to the nearest 2nd decimal place)
es_moderate  <- 0.5

## Parameters ##
n_step    <- 1 # iterate in steps of 1 to have granularity in the sample size
BF_thres  <- 10   # threshold for strong evidence
prior     <- list("Cauchy", list(prior.location=0, prior.scale=sqrt(2)/2)) # uninformed prior
type      <- "t.paired"  # paired samples t-test
alternative <- "greater" # unidirectional test
cores     <- 4 # simulate faster by dividing labor across cores
n_sim     <- 1000 # number of simulations
rdm_seed  <- 2024 # fix the seed (arbitrarily set to current year)


#-#-# Run simulations #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-
sim.H1_moderate <- BFDA.sim(expected.ES = es_moderate, type = type, 
                             prior = prior, alternative = alternative, 
                             boundary = BF_thres, cores = cores, 
                             stepsize = n_step, B=n_sim, seed = rdm_seed)

#-#-# Plot results #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #
BFDA.analyze(sim.H1_moderate, boundary=BF_thres)

#-###

T <- sim.H1_moderate[['sim']]
# Extract unique run IDs
run_id <- unique(T$id)
num_runs <- numeric(length(run_id))

# Calculate number of runs for each ID
for (i in seq_along(run_id)) {
  id <- run_id[i]
  idx_id <- which(T$id == id)
  num_runs[i] <- length(idx_id)
}

# Calculate the 95th percentile for num_runs (- 1 to 1st crossing) (+10 because default n_min is 10)
n_max_bf10_95 <- ceiling(quantile(num_runs - 1, 0.95)) + 10
n_max_bf10_80 <- ceiling(quantile(num_runs - 1, 0.80)) + 10

# Initialize study_pass matrix
study_pass <- matrix(NA, nrow = (n_max_bf10_95 - 10), ncol = 1000) # studies that pass the BF>3 at end point

# Calculate study_pass values
for (i in 1:(n_max_bf10_95 - 10)) {                                             # try different n_max
  for (j in 1:1000) {                                                           # for each simulation check whether BF>=3 with this n_max
    idx_id <- which(T$id == run_id[j])                                          # extract indexes of the j_th simulation
    if (length(idx_id) >= i) {                                                  # if the sim has more rounds than n_max (i) it means BF<10
      temp_dat <- exp(T$logBF[idx_id])                                          # Extract the BF
      if (temp_dat[i] >= 3) {                                                   # Check the BF at n_max (i)
        study_pass[i, j] <- 1                                                   # Here the study reached n_max with 3<=BF<10
      } else {
        study_pass[i, j] <- 0                                                   # Here the study terminated and BF<3
      }
    } else {
      study_pass[i, j] <- 1                                                     # here the BF_10 has been crossed so by default BF>3
    }
  }
}

# Calculate n_max_bf3_95
n_max_bf3_95 <- which(rowMeans(study_pass) > 0.95)[1] - 1 + 10

# Find best nmax to ensure either 80% power at strong evidence or 95% power at moderate evidence
max(n_max_bf3_95, n_max_bf10_80)
