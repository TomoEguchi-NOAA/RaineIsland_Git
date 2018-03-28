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

model.v.I <- 'v15AR'
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

all.yrs <- data.frame(season = seq(from = min(tally.list.3$uniq.yrs), 
                                   to = max(tally.list.3$uniq.yrs), 
                                   by = 1))

TURT.df <- as.data.frame(tally.list.3$TURT) %>%
  mutate(., season = V1) %>%
  select(., -contains('V1')) %>%
  right_join(., all.yrs, by = "season")%>%
  select(., -contains('season')) 

#colnames(TURT.df) <- c('SEASON', 1:6)

n_sec.df <- as.data.frame(tally.list.3$n_sectors) %>%
  mutate(., season = V1) %>%
  select(., -contains('V1')) %>%
  right_join(., all.yrs, by = "season")%>%
  select(., -contains('season')) 

#colnames(n_sec.df) <- c('SEASON', 1:6)

# TC.df <- as.data.frame(tally.list.3$TC)
# TC.df$SEASON <- tally.list.3$uniq.yrs

clutches.df <- as.data.frame(tally.list.3$clutches)%>%
  mutate(., season = V1) %>%
  select(., -contains('V1')) %>%
  right_join(., all.yrs, by = "season")%>%
  select(., -contains('season')) 

#colnames(clutches.df)  <- c('SEASON', 1:6)

# merge TC and TURT to create empty rows for years without
# TURT counts
# TC.TURT <- left_join(TC.df, TURT.df,
#                      by = 'SEASON')

# TURT <- select(TC.TURT, -starts_with('V'))
# TURT <- select(TURT, -contains('SEASON'))

# Do the same with n_section:
# TC.n_sec <- left_join(TC.df, n_sec.df,
#                      by = 'SEASON')

# also with clutches
# TC.clutches <- left_join(TC.df, clutches.df,
#                       by = 'SEASON')
# clutches <- select(TC.clutches, -starts_with('V'))
# clutches <- select(clutches, -contains('SEASON'))
# 
# n_sectors <- select(TC.n_sec, -starts_with('V'))
# n_sectors <- select(n_sectors, -contains('SEASON'))

#n_dates are for clutches
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

TC_n.df[is.na(TC_n.df), 'n'] <- 0

tally.n.df <- left_join(all.yrs, data.frame(n = tally.list.3$tally.n,
                                            season = tally.list.3$uniq.yrs),
                        by = 'season') 

tally.n.df[is.na(tally.n.df$n), 'n'] <- 0

# Lincoln-Petersen part
LP.list <- get.LP.data(petersen.file)

M.mat.df <- data.frame(season = LP.list$uniq.yrs,
                       LP.list$M.mat) %>%
  right_join(., all.yrs, by = 'season')%>%
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
uniq.yrs <- as.vector(all.yrs$season)
n.mat.df.0 <- n.mat.df
n.mat.df.0[is.na(n.mat.df.0)] <- 0
Nmin.vec <- apply(n.mat.df.0, FUN=max,
                  MARGIN = 1, na.rm=T) + 1

M.mat.df.0 <- M.mat.df
M.mat.df.0[is.na(M.mat.df.0)] <- 0
M.min.vec <- apply(M.mat.df.0, FUN = max, MARGIN = 1, na.rm=T)
####################################################################
# Jags starts here:

tic <- Sys.time()
bugs.data.I <- list(n = as.matrix(n.mat.df),
                    m = as.matrix(m.mat.df),
                    M = as.matrix(M.mat.df),
                    nt = LP.n,
                    #M.max = M.min.vec,
                    clutches = as.matrix(clutches.df),
                    TURT = as.matrix(TURT.df),
                    Nmin = Nmin.vec,
                    T = length(uniq.yrs),
                    TC_n = as.vector(TC_n.df$n),
                    n_dates = as.vector(n_dates.df$n_dates),
                    TC = as.matrix(TC.df))

parameters.I <- c( 's', 'N', "TC_se",
                   "s_alpha", "s_beta", 
                   "q", "q_se", 
                   "q_beta0", "q_beta1", 
                   "F_island", "N_se",
                   "alpha_1", "deviance")

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

