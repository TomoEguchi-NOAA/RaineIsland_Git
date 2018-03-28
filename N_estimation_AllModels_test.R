#N_estimation
# tests if lack of CMR data will end up in the same resultls as
# Tally data only when the integrated model contains all data
# sources. 

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

model.v.I <- 'v13'
model.v.T <- 'v7'
model.v.P <- 'v2'

petersen.file <- 'data/Petersen_data.csv'
nest.file = 'data/nesting_success.csv'
tally.file = 'data/Tally.csv'

modelfile.I <- paste0("models/Model_Petersen_Tally_AllYrs_", 
                      model.v.I, ".txt")
modelfile.P <- paste0("models/Model_Petersen_Binom_", model.v.P, ".txt" )
modelfile.T <- paste0("models/Model_Tally_", model.v.T, ".txt")

tally.list.3 <- get.tally.data.v3(nestfile = nest.file,
                                  tallyfile = tally.file)


TURT.df <- as.data.frame(tally.list.3$TURT)
colnames(TURT.df) <- c('SEASON', 1:6)

n_sec.df <- as.data.frame(tally.list.3$n_sectors)
colnames(n_sec.df) <- c('SEASON', 1:6)

TC.df <- as.data.frame(tally.list.3$TC)
TC.df$SEASON <- tally.list.3$uniq.yrs

clutches.df <- as.data.frame(tally.list.3$clutches)
colnames(clutches.df)  <- c('SEASON', 1:6)

# merge TC and TURT to create empty rows for years without
# TURT counts
TC.TURT <- left_join(TC.df, TURT.df,
                     by = 'SEASON')

TURT <- select(TC.TURT, -starts_with('V'))
TURT <- select(TURT, -contains('SEASON'))

# Do the same with n_section:
TC.n_sec <- left_join(TC.df, n_sec.df,
                     by = 'SEASON')

# also with clutches
TC.clutches <- left_join(TC.df, clutches.df,
                      by = 'SEASON')
clutches <- select(TC.clutches, -starts_with('V'))
clutches <- select(clutches, -contains('SEASON'))

n_sectors <- select(TC.n_sec, -starts_with('V'))
n_sectors <- select(n_sectors, -contains('SEASON'))

n_dates.df <- as.data.frame(tally.list.3$n_dates)
colnames(n_dates.df) <- c('SEASON', 'n_dates')
TC.n_dates <- left_join(TC.df, n_dates.df,
                        by = 'SEASON')
n_dates <- select(TC.n_dates, -starts_with('V'))
n_dates <- select(n_dates, -contains('SEASON'))

n_dates[is.na(n_dates)] <- 0

TC <- tally.list.3$TC
TC_n <- tally.list.3$tally.n

tally.n.df <- data.frame(n = tally.list.3$tally.n,
                         season = tally.list.3$uniq.yrs)

LP.list <- get.LP.data(petersen.file)
M.mat.df <- as.data.frame(LP.list$M.mat)
M.mat.df$season <- LP.list$uniq.yrs

m.mat.df <- as.data.frame(LP.list$m.mat)
m.mat.df$season <- LP.list$uniq.yrs

n.mat.df <- as.data.frame(LP.list$n.mat)
n.mat.df$season <- LP.list$uniq.yrs

# merge with Tally data and remove n and season columns:
n.mat.df <- left_join(tally.n.df, n.mat.df,
                        by = 'season')[, 3:(2+max(LP.list$LP.nt.vec))]

m.mat.df <- left_join(tally.n.df, m.mat.df,
                      by = 'season')[, 3:(2+max(LP.list$LP.nt.vec))]
M.mat.df <- left_join(tally.n.df, M.mat.df,
                      by = 'season')[, 3:(2+max(LP.list$LP.nt.vec))]
LP.n <- rowSums(!is.na(n.mat.df))

uniq.yrs <- tally.list.3$uniq.yrs

Nmin <- apply(tally.list.3$TC, 
              FUN=max, 
              MARGIN = 1, na.rm=T) + 1


# convert all non-NAs into NAs
n.mat.df[!is.na(n.mat.df)] <- NA
m.mat.df[!is.na(n.mat.df)] <- NA
M.mat.df[!is.na(n.mat.df)] <- NA
LP.n[LP.n>0]<-0

tic <- Sys.time()
bugs.data.T <- list(clutches = clutches,
                    TURT = TURT,
                    n_sectors = n_sectors,
                    Nmin = Nmin,
                    T = length(uniq.yrs),
                    TC_n = TC_n,
                    n_dates = n_dates$n_dates,
                    TC = TC)

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
     file = paste0("RData/jags_out_Tally_", model.v.T, "_", 
                   run.date, "_P_test.RData"))

# tic <- Sys.time()
bugs.data.I <- list(n = n.mat.df,
                    m = m.mat.df,
                    M = M.mat.df,
                    nt = LP.n,
                    clutches = clutches,
                    TURT = TURT,
                    Nmin = Nmin,
                    T = length(uniq.yrs),
                    TC_n = TC_n,
                    n_dates = n_dates,
                    TC = TC)

parameters.I <- c( 's', 'N', "TC_se",
                   "s_alpha", "s_beta",
                   "q", "q_se",
                   "q_beta0", "q_beta1",
                   "F_island", "deviance")
# 
# # just in case the parallel version doesn't run as before...
# load.module("dic")
# load.module("glm")
# load.module("lecuyer")
# jm <- jags.model(modelfile.I, 
#                  data = bugs.data.I, 
#                  #inits = initsFunction(),
#                  n.chains = 5,
#                  n.adapt = 600000)
# 
# tmp.I <- coda.samples(jm, 
#                       variable.names = parameters.I,
#                       n.iter = 800000,
#                       thin = 80)
# toc <- Sys.time()
# dif.time <- toc - tic
# save(bugs.data.I, parameters.I, tmp.I, dif.time, sys.info,
#      file = paste0("RData/jags_out_Integrated_", model.v.I, "_", 
#                    run.date, "_NP.RData"))

tic <- Sys.time()
# n.burnin < n.iter. Otherwise, n.iter has to be a positive integer error occurs
jags.fit <- jags.parallel(data = bugs.data.I,
                          parameters.to.save = parameters.I,
                          model.file = modelfile.I,
                          n.burnin = 600000,
                          n.chains = 5,
                          n.iter = 800000,
                          n.thin = 80,
                          jags.module = c('dic', 'glm', 'lecuyer'))

tmp.I <- as.mcmc(jags.fit)

toc <- Sys.time()
dif.time <- toc - tic
save(bugs.data.I, parameters.I, tmp.I, dif.time, sys.info,
     file = paste0("RData/jags_out_Integrated_", model.v.I, "_", 
                   run.date, "_P_test.RData"))

