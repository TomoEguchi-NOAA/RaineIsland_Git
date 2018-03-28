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

model.v.P <- 'v2AR'

petersen.file <- 'data/Petersen_data.csv'

modelfile.P <- paste0("models/Model_Petersen_Binom_", 
                      model.v.P, ".txt" )

LP.list <- get.LP.data(petersen.file)

all.yrs <- data.frame(season = seq(from = min(LP.list$uniq.yrs), 
                                   to = max(LP.list$uniq.yrs), 
                                   by = 1))

M.mat.df <- data.frame(season = LP.list$uniq.yrs,
                       LP.list$M.mat) %>%
  right_join(., all.yrs, by = 'season')%>%
  #filter(., season > 1990)
  select(., -contains('season')) 

m.mat.df <- data.frame(season = LP.list$uniq.yrs,
                       LP.list$m.mat) %>%
  right_join(., all.yrs, by = 'season')%>%
  select(., -contains('season')) 

n.mat.df <- data.frame(season = LP.list$uniq.yrs,
                       LP.list$n.mat) %>%
  right_join(., all.yrs, by = 'season')%>%
  select(., -contains('season')) 

# # merge with Tally data and remove n and season columns:
# n.mat.df <- left_join(tally.n.df, n.mat.df,
#                         by = 'season') %>%
#   select(., -contains('season')) 
# 
# m.mat.df <- left_join(tally.n.df, m.mat.df,
#                       by = 'season')[, 3:(2+max(LP.list$LP.nt.vec))]
# M.mat.df <- left_join(tally.n.df, M.mat.df,
#                       by = 'season')[, 3:(2+max(LP.list$LP.nt.vec))]
LP.n <- rowSums(!is.na(n.mat.df))
LP.n[LP.n == 0] <- 1
uniq.yrs <- as.vector(all.yrs$season)
n.mat.df.0 <- n.mat.df
n.mat.df.0[is.na(n.mat.df.0)] <- 0
Nmin.vec <- apply(n.mat.df.0, FUN=max,
                  MARGIN = 1, na.rm=T) + 1
n.mat.df.100 <- n.mat.df
n.mat.df.100[is.na(n.mat.df)] <- 1000

M.mat.df.0 <- M.mat.df
M.mat.df.0[is.na(M.mat.df.0)] <- 0
M.min.vec <- apply(M.mat.df.0, 
                   FUN = max, 
                   MARGIN = 1, 
                   na.rm=T)

M.mat.df.1 <- M.mat.df
M.mat.df.1[is.na(M.mat.df.1)] <- 1

m.mat.df.0 <- m.mat.df
m.mat.df.0[is.na(m.mat.df.0)] <- 0

zero.mat <- matrix(data = 0, nrow = nrow(M.mat.df), ncol = ncol(M.mat.df))
####################################################################
# Jags starts here:
tic <- Sys.time()
bugs.data.P <- list(n = as.matrix(n.mat.df.100),
                    m = as.matrix(m.mat.df.0),
                    M = as.matrix(M.mat.df.1),
                    nt = LP.n,
                    T = length(uniq.yrs),
                    Nmin = Nmin.vec,
                    zero = zero.mat)

parameters.P <- c("N", "deviance", "alpha_1", "N_se")

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

