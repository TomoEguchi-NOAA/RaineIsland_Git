#N_estimation


# try to estimate the total abundance of nesting females from 
# Petersen estimates, Tally counts, and Nesting success estimations

# Tally counts and nesting success are used as a combined dataset
# to model the total abundance, whereas Petersen estimates are used
# independently but the total abundance (N) is shared between the
# two models.

rm(list=ls())

run.date <- Sys.Date()
sys.info <- Sys.info()

source('RaineIsl_functions.R')
library(rjags)
#library(runjags)
library(bayesplot)
library(R2jags)

save.results <- T

model.v.I <- 'v15'
model.v.T <- 'v7'
model.v.P <- 'v2'

petersen.file <- 'data/Petersen_data.csv'
nest.file = 'data/nesting_success.csv'
tally.file = 'data/Tally.csv'

modelfile.I <- paste0("models/Model_Petersen_Tally_AllYrs_", 
                      model.v.I, ".txt")
modelfile.P <- paste0("models/Model_Petersen_Binom_", 
                      model.v.P, ".txt" )
modelfile.T <- paste0("models/Model_Tally_", 
                      model.v.T, ".txt")

tally.list.3 <- get.tally.data.v3(nestfile = nest.file,
                                  tallyfile = tally.file)

all.yrs <- data.frame(season = tally.list.3$uniq.yrs)

TURT.df <- as.data.frame(tally.list.3$TURT) %>%
  mutate(., season = V1) %>%
  select(., -contains('V1')) %>%
  right_join(., all.yrs, by = "season")%>%
  select(., -contains('season')) 

n_sec.df <- as.data.frame(tally.list.3$n_sectors) %>%
  mutate(., season = V1) %>%
  select(., -contains('V1')) %>%
  right_join(., all.yrs, by = "season")%>%
  select(., -contains('season')) 

clutches.df <- as.data.frame(tally.list.3$clutches)%>%
  mutate(., season = V1) %>%
  select(., -contains('V1')) %>%
  right_join(., all.yrs, by = "season")%>%
  select(., -contains('season')) 

n_dates.df <- as.data.frame(tally.list.3$n_dates) %>%
  transmute(., season = V1, n_dates = V2) %>%
  right_join(., all.yrs, by = "season")%>%
  select(., -contains('season')) 

n_dates.df[is.na(n_dates.df$n_dates), 'n_dates'] <- 0

TC.df <- data.frame(season = tally.list.3$uniq.yrs,
                    tally.list.3$TC) %>%
  right_join(., all.yrs, by = 'season')%>%
  select(., -contains('season')) 

TC_n.df <- data.frame(season = tally.list.3$uniq.yrs,
                      n = tally.list.3$tally.n)%>%
  right_join(., all.yrs, by = "season")%>%
  select(., -contains('season')) 

tally.n.df <- left_join(all.yrs, data.frame(n = tally.list.3$tally.n,
                                            season = tally.list.3$uniq.yrs),
                        by = 'season') 

# Lincoln-Petersen part
LP.list <- get.LP.data(petersen.file)

####################################################################
# Jags starts here:
tic <- Sys.time()
bugs.data.P <- list(n = LP.list$n.paint.mat,
                    m = LP.list$m.paint.mat,
                    M = LP.list$M.paint.mat,
                    nt = LP.list$LP_paint_n$paint,
                    Nmin = apply(LP.list$n.paint.mat, 
                                 FUN=max, 
                                 MARGIN = 1, na.rm=T) + 1,
                    T = length(LP.list$uniq.yrs))
parameters.P <- c("N", "deviance")

jags.fit <- jags.parallel(data = bugs.data.P,
                          #inits = initsFunction(),
                          parameters.to.save = parameters.P,
                          model.file = modelfile.P,
                          n.burnin = 100000,
                          n.chains = 5,
                          n.iter = 200000,
                          n.thin = 50,
                          jags.module = c('dic', 'glm', 'lecuyer'))

tmp.P <- as.mcmc(jags.fit)

toc <- Sys.time()
dif.time <- toc - tic
save(bugs.data.P, parameters.P, tmp.P, dif.time, sys.info,
     file = paste0("RData/jags_out_Petersen_", model.v.P, "_",
                   run.date, "_P.RData"))

print('Done with Petersen')

tic <- Sys.time()
bugs.data.T <- list(clutches = as.matrix(clutches.df),
                    TURT = as.matrix(TURT.df),
                    Nmin = apply(tally.list.3$TC, 
                                 FUN=max, 
                                 MARGIN = 1, na.rm=T) + 1,
                    T = length(LP.list$uniq.yrs),
                    TC_n = TC_n.df$n,
                    n_dates = n_dates.df$n_dates,
                    TC = as.matrix(TC.df))

parameters.T <- c( 's', 'N', "TC_se",
                   "s_alpha", "s_beta",
                   "q", "q_se",
                   "q_beta0", "q_beta1",
                   "F_island", "deviance")

jags.fit <- jags.parallel(data = bugs.data.T,
                          parameters.to.save = parameters.T,
                          model.file = modelfile.T,
                          n.burnin = 550000,
                          n.chains = 5,
                          n.iter = 750000,
                          n.thin = 75,
                          jags.module = c('dic', 'glm', 'lecuyer'))

tmp.T <- as.mcmc(jags.fit)

toc <- Sys.time()
dif.time <- toc - tic
save(bugs.data.T, parameters.T, tmp.T, dif.time, sys.info,
     file = paste0("RData/jags_out_Tally_", model.v.T, "_", run.date, "_P.RData"))
print('Done with Tally')

tic <- Sys.time()
bugs.data.I <- list(n = LP.list$n.paint.mat,
                    m = LP.list$m.paint.mat,
                    M = LP.list$M.paint.mat,
                    nt = LP.list$LP_paint_n$paint,
                    clutches = as.matrix(clutches.df),
                    TURT = as.matrix(TURT.df),
                    Nmin = apply(tally.list.3$TC, 
                                 FUN=max, 
                                 MARGIN = 1, na.rm=T) + 1,
                    T = length(LP.list$uniq.yrs),
                    TC_n = TC_n.df$n,
                    n_dates = n_dates.df$n_dates,
                    TC = as.matrix(TC.df))

parameters.I <- c( 's', 'N', "TC_se",
                   "s_alpha", "s_beta", 
                   "q", "q_se", 
                   "q_beta0", "q_beta1", 
                   "F_island", "deviance")

# n.burnin < n.iter. Otherwise, n.iter has to be a positive integer error occurs
jags.fit <- jags.parallel(data = bugs.data.I,
                          parameters.to.save = parameters.I,
                          model.file = modelfile.I,
                          n.burnin = 700000,
                          n.chains = 5,
                          n.iter = 1000000,
                          n.thin = 100,
                          jags.module = c('dic', 'glm', 'lecuyer'))

tmp.I <- as.mcmc(jags.fit)

toc <- Sys.time()
dif.time <- toc - tic
save(bugs.data.I, parameters.I, tmp.I, dif.time, sys.info,
     file = paste0("RData/jags_out_Integrated_", model.v.I, "_", 
                   run.date, "_P.RData"))

