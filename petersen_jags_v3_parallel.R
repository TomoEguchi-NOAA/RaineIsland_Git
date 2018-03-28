#Petersen_jags

rm(list=ls())

tic <- Sys.time()
save.results <- T

source('RaineIsl_functions.R')
library(rjags)
library(R2jags)
library(runjags)

datafile <- 'data/Petersen_data.csv'
modelfile <- "models/Model_Petersen_Binom_v3.txt" 
nest.file = 'data/nesting_success.csv'
tally.file = 'data/Tally.csv'

tally.list <- get.tally.data.v3(nestfile = nest.file,
                                tallyfile = tally.file)

load.module("glm"); load.module("lecuyer")

LP.list <- get.LP.data(datafile)

bugs.data <- list(n = LP.list$n.mat,
                  m = LP.list$m.mat,
                  M = LP.list$M.mat,
                  nt = LP.list$LP.nt.vec,
                  Nmin = apply(LP.list$n.mat, FUN=max, 
                               MARGIN = 1, na.rm=T),
                  Nmax = apply(tally.list$TC, FUN = max,
                               MARGIN = 1, na.rm = T) * 50,
                  T = length(LP.list$uniq.yrs))

parameters <- c("N", "deviance")
#modelName <- paste0(modelfile)
#load.module("dic")
results <- vector(mode =  'list', 
                  length = length(LP.list$uniq.yrs))

jags.fit <- jags.parallel(data = bugs.data,
                          #inits = initsFunction(),
                          parameters.to.save = parameters,
                          model.file = modelfile,
                          n.burnin = 100000,
                          n.chains = 5,
                          n.iter = 200000,
                          n.thin = 50,
                          jags.module = c('dic', 'glm', 'lecuyer'))

tmp <- as.mcmc(jags.fit)

g.diag <- gelman.diag(tmp)

results <- combine.mcmc(tmp)
sum.res <- summary(results)
sum.tmp <- sum.res$quantiles[grep(pattern = 'N',
                                  row.names(sum.res$quantiles)),
                             c('2.5%', '50%', '97.5%')]
results.df <- as.data.frame(sum.tmp)
#results.df$season <- LP.list$uniq.yrs
#colnames(results.df) <- c('lowN', 'modeN', 'highN', 'season')
#saveRDS(results.df, file = 'rdsfiles/petersen_jags_out.rds')

if (save.results)
  save(results.df, results, g.diag, sum.res, bugs.data, tmp,
       file = paste0('RData/Petersen_Jags_v3_parallel_out_',
                     Sys.Date(), '.RData'))

# update(jm, n.iter = n.iter)

toc <- Sys.time()

