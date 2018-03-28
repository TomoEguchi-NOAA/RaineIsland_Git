#Petersen_jags

rm(list=ls())

tic <- Sys.time()
save.results <- F

datafile <- 'data/Petersen_data.csv'
modelfile <- "models/Model_M0_Augment_v1.txt" 

source('RaineIsl_functions.R')
library(rjags)
library(R2jags)
library(runjags)

load.module("glm"); load.module("lecuyer")

LP.list <- get.LP.data(datafile)

# pool all at first:
M.all <- rowSums(LP.list$M.paint.mat, na.rm = TRUE)
n.all <- rowSums(LP.list$n.paint.mat, na.rm = TRUE)
m.all <- rowSums(LP.list$m.paint.mat, na.rm = TRUE)

yaug.arr <- array(data = 0, 
                  dim = c(length(LP.list$uniq.yrs),
                          20000,2))

for (t in 1:length(LP.list$uniq.yrs)){
  yaug.arr[t, 1:M.all[t], 1] <- 1
  yaug.arr[t, 1:m.all[t], 2] <- 1
}

bugs.data <- list(yaug = yaug.arr[1,,],
                  M = 20000)
                  #T = length(LP.list$uniq.yrs))

parameters <- c("N", "p", "Omega", "deviance")
#modelName <- paste0(modelfile)
#load.module("dic")
results <- vector(mode =  'list', 
                  length = length(LP.list$uniq.yrs))

initsFcn <- function(yaug, n.chains = 1){
  tmp <- vector(mode = 'list', length = n.chains)
  for (k in 1:n.chains){
    tmp[[k]] <-  list(z = rep(1, nrow(yaug)),
                      p = runif(1, 0, 1))
  }
  return(tmp)
}

# jags.fit <- jags.parallel(data = bugs.data,
#                           inits = initsFcn(yaug = yaug.arr[1,,], n.chains = 2),
#                           parameters.to.save = parameters,
#                           model.file = modelfile,
#                           n.burnin = 100000,
#                           n.chains = 2,
#                           n.iter = 200000,
#                           n.thin = 50,
#                           jags.module = c('dic', 'glm', 'lecuyer'))

jm <- jags.model(modelfile,
                 data = bugs.data,
                 inits = list(z = rep(1, 20000),
                              p = runif(1, 0, 1)),
                 n.chains = 5,
                 n.adapt = 5000)

load.module("dic")
zm <- coda.samples(jm, variable.names = parameters,
                   n.iter = 10000)

gelman.diag(zm)
bayesplot::mcmc_trace(zm, 'p')


tmp <- as.mcmc(jags.fit)

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
  save(results.df, results, g.diag, sum.res, bugs.data, tmp,
       file = paste0('RData/augmented_Jags_v2_parallel_out_',
                     Sys.Date(), '.RData'))

# update(jm, n.iter = n.iter)

toc <- Sys.time()

