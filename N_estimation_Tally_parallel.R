#N_estimation


# try to estimate the total abundance of nesting females from 
# Petersen estimates, Tally counts, and Nesting success estimations

# Tally counts and nesting success are used as a combined dataset
# to model the total abundance, whereas Petersen estimates are used
# independently but the total abundance (N) is shared between the
# two models.


rm(list=ls())
tic <- Sys.time()
source('RaineIsl_functions.R')
library(rjags)
library(runjags)
library(mcmcplots)
library(R2jags)

model.v <- 'v7'
nest.file <- 'data/nesting_success.csv'
tally.file <- 'data/Tally.csv'
modelName <- paste0("models/Model_Tally_", model.v, ".txt")

run.date <- Sys.Date()

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

uniq.yrs <- tally.list.3$uniq.yrs

# save(tally.list, LP.list,
#      file = paste0('RData/N_estimaion_AllYrs_data_v7_',
#                    Sys.Date(), '.RData'))

# jags starts here:
n.chains <- 5
n.adapt <- 900000
n.iter <- 100000
n.thin <- 50

load.module("dic")
load.module("glm")
load.module("lecuyer")

Nmin <- apply(tally.list.3$TC, 
              FUN=max, 
              MARGIN = 1, na.rm=T) + 1

bugs.data <- list(clutches = clutches,
                  TURT = TURT,
                  n_sectors = n_sectors,
                  Nmin = Nmin,
                  T = length(uniq.yrs),
                  TC_n = TC_n,
                  n_dates = as.vector(n_dates$n_dates),
                  TC = TC)  

parameters <- c( 's', 'N', "TC_se",
                 "s_alpha", "s_beta", 
                 "q", "q_se", 
                 "q_beta0", "q_beta1", 
                 "F_island", "deviance")

# jm <- jags.model(modelName, 
#                  data = bugs.data, 
#                  inits = initsFunction(),
#                  n.chains = n.chains,
#                  n.adapt = n.adapt)

# tmp <- coda.samples(jm, 
#                     variable.names = parameters,
#                     n.iter = n.iter,
#                     thin = n.thin)

jags.fit <- jags.parallel(data = bugs.data,
                          parameters.to.save = parameters,
                          model.file = modelName,
                          n.burnin = 550000,
                          n.chains = 5,
                          n.iter = 750000,
                          n.thin = 75,
                          jags.module = c('dic', 'glm', 'lecuyer'))

tmp <- as.mcmc(jags.fit)

# check for convergence
g.diag <- gelman.diag(tmp)

sum.tmp <- summary(tmp)
N.df <- data.frame(sum.tmp$quantiles[grep(pattern = 'N[/[]',
                                          row.names(sum.tmp$quantiles)), 
                                     c('2.5%', '50%', '97.5%')])
colnames(N.df) <- c('lowN', 'modeN', 'highN')

# R.df <- data.frame(sum.tmp$quantiles[grep(pattern = 'R',
#                                           row.names(sum.tmp$quantiles)), 
#                                      c('2.5%', '50%', '97.5%')])
# colnames(R.df) <- c('lowR', 'modeR', 'highR')

q.df <- data.frame(sum.tmp$quantiles[grep(pattern = 'q',
                                          row.names(sum.tmp$quantiles)), 
                                     c('2.5%', '50%', '97.5%')])

s_alpha.df <- data.frame(sum.tmp$quantiles[grep(pattern = 's_alpha',
                                          row.names(sum.tmp$quantiles)), 
                                     c('2.5%', '50%', '97.5%')])

s_beta.df <- data.frame(sum.tmp$quantiles[grep(pattern = 's_beta',
                                                row.names(sum.tmp$quantiles)), 
                                           c('2.5%', '50%', '97.5%')])
# results.s.df<- data.frame(sum.tmp$quantiles[grep(pattern = 's',
#                                                  row.names(sum.tmp$quantiles)), 
#                                             c('2.5%', '50%', '97.5%')])

N.df$season <- uniq.yrs

dt <- 6
mv.sums <- vector(mode = 'numeric', 
                  length = ceiling(nrow(N.df)))
for (k in 1:(nrow(N.df) - dt)){
  mv.sums[k] <- sum(N.df[k:(k+dt-1), 'modeN'])
}

toc <- Sys.time()

save(tally.list.3,
     file = paste0('RData/N_estimaion_Tally_data_', 
                   model.v, '_',
                   run.date, '.RData'))

save(N.df, q.df, tmp, mv.sums,
     bugs.data, sum.tmp, tic, toc,
     file = paste0('RData/Tally_', model.v, '_Jags_',
                   run.date, '.RData'))
