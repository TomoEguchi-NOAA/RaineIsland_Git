# Nhat using drone and boat-based painted turtle counts with 
# binomial likelihood
# 

rm(list=ls())
library(jagsUI)
library(tidyverse)
library(bayesplot)

m.vector <- c(62, 19, 7, 149, 131)
n.vector <- c(247, 80, 33, 1351, 792 )
idx.vector <- c(1, 1, 1, 2, 2)

jags.data <- list(M = 587,
                  m = m.vector,
                  n = n.vector,
                  idx = idx.vector,
                  n.survey = length(m.vector),
                  n.methods = 2,
                  N.min = 1000)


#"p.vessel",
MCMC.params <- list(n.samples = 50000,
                    n.thin = 50,
                    n.burnin = 20000,
                    n.chains = 5)

# v1 - one N and both methods try to estimate it
jags.params.1 <- c("N", "deviance", "log.lkhd")
jm.1 <- jags(data = jags.data,
           parameters.to.save = jags.params.1,
           model.file = "models/model_Petersen_UAS_Vessel_v1.txt",
           n.chains = MCMC.params$n.chains,
           n.burnin = MCMC.params$n.burnin,
           n.iter = MCMC.params$n.samples,
           parallel = T,
           DIC = T)  

mcmc_trace(jm.1$samples,
           c( "N"))

mcmc_dens(jm.1$samples,
          c(  "N"))


# v1 - one N and both methods try to estimate it
jags.params.2 <- c("N[1]", "N[2]", "N.B1", "deviance", "log.lkhd")
jm.2 <- jags(data = jags.data,
             parameters.to.save = jags.params.2,
             model.file = "models/model_Petersen_UAS_Vessel_v2.txt",
             n.chains = MCMC.params$n.chains,
             n.burnin = MCMC.params$n.burnin,
             n.iter = MCMC.params$n.samples,
             parallel = T,
             DIC = T)  

mcmc_trace(jm.2$samples,
           c( "N[1]", "N[2]", "N.B1"))

mcmc_dens(jm.2$samples,
           c(  "N[1]", "N[2]", "N.B1"))
