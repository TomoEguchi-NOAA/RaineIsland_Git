#Petersen_jags

rm(list=ls())

tic <- Sys.time()
save.results <- F

datafile <- 'data/Petersen_data.csv'
modelfile <- "models/Model_Petersen_Binom_v2.txt" 

source('RaineIsl_functions.R')
library(rjags)
library(runjags)

load.module("glm"); load.module("lecuyer")

LP.list <- get.LP.data(datafile)

bugs.data <- list(n = LP.list$n.mat,
                  m = LP.list$m.mat,
                  M = LP.list$M.mat,
                  nt = LP.list$LP.nt.vec,
                  Nmin = apply(LP.list$n.mat, FUN=max, 
                               MARGIN = 1, na.rm=T),
                  T = length(LP.list$uniq.yrs))

initsFunction <- function(){
  N <- runif(n = length(LP.list$uniq.yrs), 
             apply(LP.list$n.mat, FUN=max, 
                   MARGIN = 1, na.rm=T), 
             100000)
  A <- list(N = N, 
            .RNG.name="lecuyer::RngStream")
  return(A)
  #.RNG.seed=runif(n = 1, min = 1, max = 1e+06)
}  

n.chains <- 5
n.adapt <- 100000
n.iter <- 200000
n.thin <- 50

parameters <- c("N", "deviance")
#modelName <- paste0(modelfile)
load.module("dic")
results <- vector(mode =  'list', 
                  length = length(LP.list$uniq.yrs))

jm <- jags.model(modelfile, 
                 data = bugs.data, 
                 inits = initsFunction(),
                 n.chains = n.chains,
                 n.adapt = n.adapt)

tmp <- coda.samples(jm, 
                    variable.names = parameters,
                    n.iter = n.iter,
                    thin = n.thin)

g.diag <- gelman.diag(tmp)

results <- combine.mcmc(tmp)
sum.res <- summary(results)
sum.tmp <- sum.res$quantiles[grep(pattern = 'N',
                                  row.names(sum.res$quantiles)),
                             c('2.5%', '50%', '97.5%')]
results.df <- as.data.frame(sum.tmp)
results.df$season <- LP.list$uniq.yrs
colnames(results.df) <- c('lowN', 'modeN', 'highN', 'season')
#saveRDS(results.df, file = 'rdsfiles/petersen_jags_out.rds')

if (save.results)
  save(results.df, results, g.diag, sum.res, bugs.data,
       file = paste0('RData/Petersen_Jags_v2_out_',
                     Sys.Date(), '.RData'))

# update(jm, n.iter = n.iter)

toc <- Sys.time()
