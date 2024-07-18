# Author:         Lucas Jeay-Bizot
# Created:        July 18th 2024
# Last Modified:  July 18th 2024

# Purpose: This code, using the BFDA package, will plot Fig S1 B and Fig S2 A.
# It will take as input the effect sizes found in the pilot analyses and run
# Monte Carlo simulations to determine the distribution of studies that would
# find an effect with such an effect size using a Bayesian Sequential Sampling
# design with a maximal sample size.


# References: Schönbrodt, F. D. & Stefan, A. M. (2019). BFDA: An R package for
# Bayes factor design analysis (version 0.5.0). Retrieved from
# https://github.com/nicebread/BFDA

# Dependencies
library(BFDA) # Package for Bayes Factor Design Analysis (https://github.com/nicebread/BFDA)

## INPUT ## Effect sizes (rounded to the nearest 2nd decimal place)
es_lateRP_M1  <- 0.79
es_Mtime_M1   <- 0.58

## Parameters ##
n_max     <- 61   # maximal sample size as determined in design_SI_nmax.R
n_min     <- 21   # n_min and n_step are implemented to avoid a spurious threshold crossing
n_step    <- 5
BF_thres  <- 10   # threshold for strong evidence
prior     <- list("Cauchy", list(prior.location=0, prior.scale=sqrt(2)/2)) # uninformed prior
type      <- "t.paired"  # paired samples t-test
alternative <- "greater" # unidirectional test
cores     <- 4 # simulate faster by dividing labor across cores
n_sim     <- 1000 # number of simulations
rdm_seed  <- 2024 # fix the seed (arbitrarily set to current year)


#-#-# Run simulations #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-
sim.H1_lateRP_M1 <- BFDA.sim(expected.ES = es_lateRP_M1, type = type, 
                             prior = prior, n.min = n_min, n.max = n_max, 
                             alternative = alternative, boundary = BF_thres,
                             cores = cores, stepsize = n_step, B=n_sim,
                             seed = rdm_seed)

sim.H1_Mtime_M1 <- BFDA.sim(expected.ES = es_Mtime_M1, type = type, 
                             prior = prior, n.min = n_min, n.max = n_max, 
                             alternative = alternative, boundary = BF_thres,
                             cores = cores, stepsize = n_step, B=n_sim,
                             seed = rdm_seed)

#-#-# Plot results #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #-#-# #
BFDA.analyze(sim.H1_lateRP_M1, boundary=BF_thres)

BFDA.analyze(sim.H1_Mtime_M1, boundary=BF_thres)

plot(sim.H1_lateRP_M1)

plot(sim.H1_Mtime_M1)
