#N_BaileysBinomial
# Estimates N from Bailey's binomial model and 
# Hypergeometric model from raw data

# 15 December 2017


# rm(list=ls())
# source('RaineIsl_functions.R')

BaileyHyperGeo <- function(data.file = 'data/Petersen_data.csv'){
  
  LP.list <- get.LP.data(data.file)
  
  M <- LP.list$M.mat
  n <- LP.list$n.mat
  m <- LP.list$m.mat
  
  # These formula can be found in Seber (1982)
  # the weighted averages and their SE formula can be found here:
  # https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
  
  # Bailey's binomial model:  (I think this is what Andy Dunstan used)
  Baileys.N <- M * (n + 1) / (m + 1)
  log.Baileys.N.var <- 2 * log(M) + log(n + 1) + log(n - m) - 2 * log(m + 1) - log(m + 2)
  Baileys.N.var <- exp(log.Baileys.N.var)
  
  # Hypergeometric model:
  HG.N <- ((M + 1) * (n + 1) / (m + 1)) - 1
  HG.N.vars <- (M + 1) * (n + 1) * (M - m) * (n - m) / (((m + 1)^2) * (m + 2))
  
  Baileys.N.Wavg <- Baileys.N.Wvar <- Baileys.N.avg <- Baileys.SE <- vector(mode = 'numeric', 
                                                              length = length(LP.list$uniq.yrs))
  HG.N.Wavg <- HG.N.Wvar <- HG.N.avg <- HG.SE <- vector(mode = 'numeric', 
                                               length = length(LP.list$uniq.yrs))
  
  for (k in 1:length(LP.list$uniq.yrs)){
    Baileys.Nhats <- Baileys.N[k, 1:LP.list$LP.nt.vec[k]]
    Baileys.vars <-  Baileys.N.var[k, 1:LP.list$LP.nt.vec[k]]
    Baileys.N.Wavg[k] <- weighted.mean(Baileys.Nhats, Baileys.vars)
    Baileys.N.Wvar[k] <- 1/sum(1/Baileys.vars)
    Baileys.N.avg[k] <- mean(Baileys.Nhats, na.rm = T)
    Baileys.SE[k] <- ifelse(length(Baileys.Nhats) > 1, 
                            SE(Baileys.Nhats),
                            sqrt(Baileys.vars))
    
    HG.Nhats <- HG.N[k, 1:LP.list$LP.nt.vec[k]]
    HG.vars <-  HG.N.vars[k, 1:LP.list$LP.nt.vec[k]]
    HG.N.Wavg[k] <- weighted.mean(HG.Nhats, HG.vars)
    HG.N.Wvar[k] <- 1/sum(1/HG.vars)
    HG.N.avg[k] <- mean(HG.Nhats, na.rm = T)
    HG.SE[k] <- ifelse(length(HG.Nhats) > 1, 
                       SE(HG.Nhats),
                       sqrt(HG.vars))
  }
  
  N.avg.df <- data.frame(Nhat.Baileys.w = Baileys.N.Wavg,
                         N.SE.Baileys.w = sqrt(Baileys.N.Wvar),
                         Nhat.Baileys = Baileys.N.avg,
                         N.SE.Baileys = Baileys.SE,
                         Nhat.HG.w = HG.N.Wavg,
                         N.SE.HG.w = sqrt(HG.N.Wvar),
                         Nhat.HG = HG.N.avg,
                         N.SE.HG = HG.SE,
                         season = LP.list$uniq.yrs,
                         n = LP.list$LP.nt.vec)
  return(N.avg.df)
}
