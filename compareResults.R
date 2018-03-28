# Results comparison


# Compare results from different approaches.


# load results from different approaches:
rm(list=ls())
source('RaineIsl_functions.R')
library(bayesplot)
library(tidyverse)
library(ggridges)

model.v.I <- 'v15'
model.v.T <- 'v7'
model.v.P <- 'v2'

# 1/31/2018 analyses used data pooled by paint colors rather than
# treating them all different for Petersen estimates. All previous analyses were based on
# per sampling. 
run.date.int <- '2018-01-31'
run.date.P <- '2018-01-31'
run.date.Tally <- '2018-01-31'

saveFig <- F

# input data files 
petersen.file <- 'data/Petersen_data.csv'
nest.file <- 'data/nesting_success.csv'
tally.file <- 'data/Tally.csv'

# jags output files are here: _P designated parallel runs:
petersen.outfile <- paste0('RData/jags_out_Petersen_', 
                           model.v.P,'_', run.date.P, '_P.RData')
integrated.outfile <- paste0('RData/jags_out_Integrated_', model.v.I, '_', 
                             run.date.int,  '_P.RData')
tally.outfile <- paste0('RData/jags_out_Tally_', model.v.T, '_',
                        run.date.Tally, '_P.RData')

tally.list.3 <- get.tally.data.v3(nestfile = nest.file,
                                  tallyfile = tally.file)
LP.list <- get.LP.data(petersen.file)

add.season <- function(df.in, conv.table){
  df.in$season <- NA
  for (k in 1:nrow(df.in)){
    df.in$season[k] <- filter(conv.table, 
                              parameter == df.in[k, 'parameter'])$season
  }
  return(df.in)
}

params2seasons.N <- data.frame(parameter = c("N[1]", "N[2]", "N[3]", "N[4]", 
                                           "N[5]", "N[6]", "N[7]", "N[8]", "N[9]", 
                                           "N[10]", "N[11]", "N[12]", "N[13]", 
                                           "N[14]", "N[15]", "N[16]", "N[17]", 
                                           "N[18]", "N[19]", "N[20]", "N[21]", 
                                           "N[22]", "N[23]", "N[24]", "N[25]", 
                                           "N[26]", "N[27]", "N[28]", "N[29]", 
                                           "N[30]", "N[31]", "N[32]", "N[33]", 
                                           "N[34]", "N[35]", "N[36]", "N[37]", 
                                           "N[38]", "N[39]", "N[40]"),
                             season = tally.list.3$uniq.yrs)

params2seasons.s <- data.frame(parameter = c("s[1]", "s[2]", "s[3]", "s[4]", 
                                             "s[5]", "s[6]", "s[7]", "s[8]", "s[9]", 
                                             "s[10]", "s[11]", "s[12]", "s[13]", 
                                             "s[14]", "s[15]", "s[16]", "s[17]", 
                                             "s[18]", "s[19]", "s[20]", "s[21]", 
                                             "s[22]", "s[23]", "s[24]", "s[25]", 
                                             "s[26]", "s[27]", "s[28]", "s[29]", 
                                             "s[30]", "s[31]", "s[32]", "s[33]", 
                                             "s[34]", "s[35]", "s[36]", "s[37]", 
                                             "s[38]", "s[39]", "s[40]"),
                               season = tally.list.3$uniq.yrs)

params2seasons.F_i <- data.frame(parameter = c("F_island[1]", "F_island[2]", 
                                               "F_island[3]", "F_island[4]", 
                                             "F_island[5]", "F_island[6]", 
                                             "F_island[7]", "F_island[8]", 
                                             "F_island[9]", 
                                             "F_island[10]", "F_island[11]", 
                                             "F_island[12]", "F_island[13]", 
                                             "F_island[14]", "F_island[15]", 
                                             "F_island[16]", "F_island[17]", 
                                             "F_island[18]", "F_island[19]", 
                                             "F_island[20]", "F_island[21]", 
                                             "F_island[22]", "F_island[23]", 
                                             "F_island[24]", "F_island[25]", 
                                             "F_island[26]", "F_island[27]", 
                                             "F_island[28]", "F_island[29]", 
                                             "F_island[30]", "F_island[31]", 
                                             "F_island[32]", "F_island[33]", 
                                             "F_island[34]", "F_island[35]", 
                                             "F_island[36]", "F_island[37]", 
                                             "F_island[38]", "F_island[39]", 
                                             "F_island[40]"),
                               season = tally.list.3$uniq.yrs)

params2seasons.q <- data.frame(parameter = c("q[1]", "q[2]", "q[3]", "q[4]", 
                                               "q[5]", "q[6]", "q[7]", "q[8]", "q[9]", 
                                               "q[10]", "q[11]", "q[12]", "q[13]", 
                                               "q[14]", "q[15]", "q[16]", "q[17]", 
                                               "q[18]", "q[19]", "q[20]", "q[21]", 
                                               "q[22]", "q[23]", "q[24]", "q[25]", 
                                               "q[26]", "q[27]", "q[28]", "q[29]", 
                                               "q[30]", "q[31]", "q[32]", "q[33]", 
                                               "q[34]", "q[35]", "q[36]", "q[37]", 
                                               "q[38]", "q[39]", "q[40]"),
                                 season = tally.list.3$uniq.yrs)
load(integrated.outfile)

integrated.jags.out <- tmp.I

integrated.quant.df <- as.data.frame(summary(tmp.I)$quantiles) %>%
  mutate(., parameter = row.names(summary(tmp.I)$quantiles))
integrated.stats.df <- as.data.frame(summary(tmp.I)$statistics) %>%
  mutate(., parameter = row.names(summary(tmp.I)$statistics))

integrated.N <- select(integrated.stats.df[grep(integrated.stats.df$parameter, 
                                                pattern = 'N[/[]'),], 
                       Mean, SD) %>%
  cbind(integrated.quant.df[grep(integrated.quant.df$parameter, 
                                 pattern = 'N[/[]'),]) %>%
  add.season(., params2seasons.N)

integrated.s <- select(integrated.stats.df[grep(integrated.stats.df$parameter,
                                                pattern = 's[/[]'), ],
                       Mean, SD) %>% 
  cbind(integrated.quant.df[grep(integrated.quant.df$parameter,
                                 pattern = 's[/[]'), ])%>%
  add.season(., params2seasons.s)

integrated.F_island <- select(integrated.stats.df[grep(integrated.stats.df$parameter,
                                                       pattern = 'F_island[/[]'), ],
                              Mean, SD) %>% 
  cbind(integrated.quant.df[grep(integrated.quant.df$parameter,
                                 pattern = 'F_island[/[]'), ])%>%
  add.season(., params2seasons.F_i)

integrated.s_alpha <- select(integrated.stats.df[grep(integrated.stats.df$parameter,
                                                      pattern = 's_alpha'), ],
                             Mean, SD) %>% 
  cbind(integrated.quant.df[grep(integrated.quant.df$parameter,
                                 pattern = 's_alpha'), ])

integrated.s_beta <- select(integrated.stats.df[grep(integrated.stats.df$parameter,
                                                      pattern = 's_beta'), ],
                             Mean, SD) %>% 
  cbind(integrated.quant.df[grep(integrated.quant.df$parameter,
                                 pattern = 's_beta'), ])

if (model.v.I != 'v14'){
  integrated.q <- select(integrated.stats.df[grep(integrated.stats.df$parameter,
                                                  pattern = 'q[/[]'), ],
                         Mean, SD) %>% 
    cbind(integrated.quant.df[grep(integrated.quant.df$parameter,
                                   pattern = 'q[/[]'), ])%>%
    add.season(., params2seasons.q)
  
} else {
  integrated.q <- select(integrated.stats.df[grep(integrated.stats.df$parameter,
                                                  pattern = 'q'),],
                         Mean, SD) %>% 
    cbind(integrated.quant.df[grep(integrated.quant.df$parameter,
                                   pattern = 'q'),])
}

# get just Petersen estimates:
load(petersen.outfile)
params2seasons.CMR <- data.frame(parameter = c("N[1]", "N[2]", "N[3]", "N[4]", "N[5]", 
                                               "N[6]", "N[7]", "N[8]", "N[9]", "N[10]", 
                                               "N[11]", "N[12]", "N[13]", "N[14]", "N[15]", 
                                               "N[16]", "N[17]", "N[18]", "N[19]", "N[20]", 
                                               "N[21]", "N[22]", "N[23]", "N[24]"),
                                 season = LP.list$uniq.yrs)

Petersen.jags.out <- tmp.P
summary.N <- as.data.frame(summary(tmp.P)$statistics[grep(row.names(summary(tmp.P)$statistics),
                                                   pattern = 'N[/[]'),]) 
summary.N$parameter <- row.names(summary.N)  
summary.N <- add.season(summary.N, params2seasons.CMR)
  
Bayes.Petersen.quant <- as.data.frame(summary(tmp.P)$quantiles) %>% 
  mutate(parameter = row.names(summary(tmp.P)$quantiles))

Bayes.Petersen.stats <- as.data.frame(summary(tmp.P)$statistics) %>% 
  mutate(parameter = row.names(summary(tmp.P)$statistics))

Bayes.Petersen.N <- select(Bayes.Petersen.stats[grep(Bayes.Petersen.stats$parameter, 
                                              pattern = 'N[/[]'),], 
                     Mean, SD) %>%
  cbind(Bayes.Petersen.quant[grep(Bayes.Petersen.quant$parameter, 
                                 pattern = 'N[/[]'),]) %>%
  add.season(., params2seasons.CMR)

# just Tally data:
TC.df <- data.frame(mean = rowMeans(tally.list.3$TC, na.rm = TRUE),
                    SE = apply(tally.list.3$TC, MARGIN = 1, FUN = SE),
                    season = tally.list.3$uniq.yrs,
                    TC.n = tally.list.3$tally.n)

TC.df[TC.df$TC.n == 1, 'SE'] <- 0

# Tally estimates:
load(tally.outfile)

Tally.jags.out <- tmp.T

Tally.quant <- as.data.frame(summary(tmp.T)$quantiles) %>%
  mutate(., parameter = row.names(summary(tmp.T)$quantiles))

Tally.stats <- as.data.frame(summary(tmp.T)$statistics) %>%
  mutate(., parameter = row.names(summary(tmp.T)$statistics))

Tally.N <- select(Tally.stats[grep(Tally.stats$parameter, 
                                   pattern = 'N[/[]'),], 
                  Mean, SD) %>%
  cbind(Tally.quant[grep(Tally.quant$parameter, 
                                 pattern = 'N[/[]'),]) %>%
  mutate(season = tally.list.3$uniq.yrs)

Tally.s <- select(Tally.stats[grep(Tally.stats$parameter, 
                                   pattern = 's[/[]'),], 
                  Mean, SD) %>%
  cbind(Tally.quant[grep(Tally.quant$parameter, 
                         pattern = 's[/[]'),]) %>%
  mutate(season = tally.list.3$uniq.yrs)

Tally.F_island <- select(Tally.stats[grep(Tally.stats$parameter, 
                                   pattern = 'F_island[/[]'),], 
                  Mean, SD) %>%
  cbind(Tally.quant[grep(Tally.quant$parameter, 
                         pattern = 'F_island[/[]'),]) %>%
  mutate(season = tally.list.3$uniq.yrs)

Tally.s_alpha <- select(Tally.stats[grep(Tally.stats$parameter, 
                                   pattern = 's_alpha'),], 
                  Mean, SD) %>%
  cbind(Tally.quant[grep(Tally.quant$parameter, 
                         pattern = 's_alpha'),])

Tally.s_beta <- select(Tally.stats[grep(Tally.stats$parameter, 
                                         pattern = 's_beta'),], 
                        Mean, SD) %>%
  cbind(Tally.quant[grep(Tally.quant$parameter, 
                         pattern = 's_beta'),])

if (model.v.T != 'v7'){
  Tally.q <- select(Tally.stats[grep(Tally.stats$parameter, 
                                     pattern = 'q[/[]'),], 
                    Mean, SD) %>%
    cbind(Tally.quant[grep(Tally.quant$parameter, 
                           pattern = 'q[/[]'),]) %>%
    mutate(season = N.df$season)
  
} else {
  Tally.q <- select(Tally.stats[grep(Tally.stats$parameter, 
                                     pattern = 'q'),], 
                    Mean, SD) %>%
    cbind(Tally.quant[grep(Tally.quant$parameter, 
                           pattern = 'q'),])
  
}

# get Bailey and Hypergeo results here:
Bailey_HG <- BaileyHyperGeo('data/Petersen_data.csv')

# comparison between nesting success estimates between Tally only
# and Integrated analysis:
dat.int.Tally.s <- data.frame(Tally.med = Tally.s$`50%`,
                              Tally.low = Tally.s$`2.5%`,
                              Tally.high = Tally.s$`97.5%`,
                              Tally.mean = Tally.s$Mean,
                              Tally.SE = Tally.s$SD,
                              int.med = integrated.s$`50%`,
                              int.low = integrated.s$`2.5%`,
                              int.high = integrated.s$`97.5%`,
                              int.mean = integrated.s$Mean,
                              int.SE = integrated.s$SD,
                              season = integrated.s$season)

plot.int.Tally.s <- ggplot(data = dat.int.Tally.s) + 
  geom_point(aes(x = Tally.med, y = int.med))

plot(Tally.s$`50%`, integrated.s$`50%`)
plot(Tally.s$`97.5%`, integrated.s$`97.5%`)

# comparison between the number of females on the island between
# Tally and integrated analysis
plot(Tally.F_island$Mean, integrated.F_island$Mean)

# comparison of Nhats among three approaches:
plot(Bailey_HG$season, Bailey_HG$Nhat.Baileys)
points(integrated.N$season, integrated.N$`50%`, col = 'red')

plot(Bailey_HG$Nhat.Baileys, Bayes.Petersen.N$`50%`)
plot(integrated.N$season, integrated.N$`50%`)
points(Bayes.Petersen.N$season, Bayes.Petersen.N$`50%`, col = 'red')
