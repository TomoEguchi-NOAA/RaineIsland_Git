#plot_Petersen_TAlly_results


rm(list=ls())
source('RaineIsl_functions.R')
library(bayesplot)

model.v <- 'v14'
model.v.Tally <- 'v7'

saveFig <- F
run.date <- '2018-01-05'
load(file = paste0('RData/Petersen_Tally_AllYr_', model.v, 
                   '_parallel_out_', 
                   run.date, '.RData'))

summary.params.df <- as.data.frame(sum.tmp$quantiles)
summary.params.df$parameter <- row.names(sum.tmp$quantiles)
summary.params.df$SE <- as.data.frame(sum.tmp$statistics)$SD

summary.N <- summary.params.df[grep(summary.params.df$parameter, 
                                    pattern = 'N[/[]'),]
# 
# summary.r <- summary.params.df[grep(summary.params.df$parameter, 
#                                     pattern = 'r[/[]'), ]
summary.s <- summary.params.df[grep(summary.params.df$parameter,
                                    pattern = 's[/[]'), ]

# summary.alpha_s <- sum.tmp$quantiles['alpha_s',]
# summary.beta_s <- sum.tmp$quantiles['beta_s',]
if (model.v != 'v.14'){
  summary.q <- summary.params.df[grep(summary.params.df$parameter,
                                      pattern = 'q[/[]'), ]
  
} else {
  summary.q <- summary.params.df[grep(summary.params.df$parameter,
                                      pattern = 'q'),]
}

summary.stats.df <- data.frame(season = N.df$season,
                               N_low_all = N.df[,'lowN'],
                               N_med_all = N.df[,'modeN'],
                               N_high_all = N.df[,'highN'],
                               N_SE = summary.N[, 'SE'],
                               q_low_all = summary.q['q','2.5%'],
                               q_med_all = summary.q['q','50%'],
                               q_high_all = summary.q['q','97.5%']) %>%
  mutate(., N_width_all = N_high_all - N_low_all)

# data used in all-year analysis
load(paste0('RData/N_estimaion_AllYrs_data_', model.v, '_parallel_',
            run.date, '.RData'))
nest.suc.df <- na.omit(data.frame(season = tally.list.3$nest.raw$SEASON,
                                  NSUC1 = tally.list.3$nest.raw$NSUC1))

# get just Petersen estimates:
load('RData/Petersen_Jags_v2_out_2018-01-02.RData')
summary.N <- as.data.frame(sum.res$statistics[grep(row.names(sum.res$statistics),
                                     pattern = 'N[/[]'),]) %>%
  mutate(season = results.df$season)

results.Petersen.df <- mutate(results.df,
                              LP.n = bugs.data$nt) %>%
  right_join(., summary.N, by = 'season') %>%
  right_join(., summary.stats.df,
             by = "season") %>%
  select(., -starts_with('q_')) %>%
  mutate(LP_width = highN - lowN) %>%
  na.omit() %>%
  mutate(SE = SD)
  
# just Tally data:
TC.df <- data.frame(mean = rowMeans(tally.list.3$TC, na.rm = TRUE),
                    SE = apply(tally.list.3$TC, MARGIN = 1, FUN = SE),
                    season = tally.list.3$uniq.yrs,
                    TC.n = tally.list.3$tally.n)

TC.df[TC.df$TC.n == 1, 'SE'] <- 0

results.all <- left_join(results.Petersen.df,
                         TC.df,
                         by = 'season')
results.all[is.na(results.all$LP.n), 'LP.n'] = 0

# Tally estimates:
load(paste0('RData/Tally_', model.v.Tally, 
            '_Jags_2018-01-08.RData'))
summary.Tally <- as.data.frame(sum.tmp$quantiles) %>%
  mutate(., parameter = row.names(sum.tmp$quantiles))
summary.Tally.s <- summary.Tally[grep(summary.Tally$parameter,
                                    pattern = 's[/[]'), ]
N.Tally.df <- N.df
if (model.v.Tally != 'v7'){
  summary.Tally.q <- summary.Tally[grep(summary.Tally$parameter,
                                        pattern = 'q[/[]'), ]
  
} else {
  summary.Tally.q <- summary.Tally[grep(summary.Tally$parameter,
                                        pattern = 'q'), ]
  
}

stats.Tally <- as.data.frame(sum.tmp$statistics) %>%
  mutate(., parameter = row.names(sum.tmp$statistics))

stats.Tally.N <- stats.Tally[grep(stats.Tally$parameter,
                                  pattern = 'N[/[]'),]

summary.Tally.stats.df <- data.frame(season = N.Tally.df$season,
                                     N_low = N.Tally.df[,'lowN'],
                                     N_med = N.Tally.df[,'modeN'],
                                     N_high = N.Tally.df[,'highN'],
                                     s_low = summary.Tally.s[,'2.5%'],
                                     s_med = summary.Tally.s[,'50%'],
                                     s_high = summary.Tally.s[,'97.5%'],
                                     SE = stats.Tally.N$SD) %>%
  mutate(., N_width = N_high - N_low)


# get Bailey and Hypergeo results here:
Bailey_HG <- BaileyHyperGeo('data/Petersen_data.csv')

# s.df <- as.data.frame(summary.s)
# s.df <- select(s.df, -starts_with('parameter'))
# s.df$season <- N.df$season
# colnames(s.df) <- c('Slow', 'S25', 'Smode',
#                     'S75', 'Shigh', 'season')
#tally.df <- cbind(TC.df, s.df)

p.Bailey <- ggplot(data = Bailey_HG) + 
  labs(y = 'Estimated abundance',
       title = 'Bailey\'s estimates') +
  geom_errorbar(aes(x = season, 
                    ymax = Nhat.Baileys + N.SE.Baileys, 
                    ymin = Nhat.Baileys - N.SE.Baileys),
                size = 1,
                color = 'black') + 
  geom_point(aes(x = season, y = Nhat.Baileys),
             size = 2.5, color = 'black',
             shape = 15) +
  
  geom_errorbar(aes(x = season+0.5, 
                    ymax = Nhat.Baileys.w + N.SE.Baileys.w, 
                    ymin = Nhat.Baileys.w - N.SE.Baileys.w),
                size = 1,
                color = 'red') + 
  geom_point(aes(x = season+0.5, y = Nhat.Baileys.w),
             size = 2.5, color = 'red',
             shape = 15) +
  
  geom_text(aes(x = season, 
                y = Nhat.Baileys.w + N.SE.Baileys.w + 2000,
                label = n,
                vjust = 'bottom'),
            color = 'black',
            fontface = 'bold',
            size = 4) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

p.Petersen <- ggplot(data = results.Petersen.df) + 
  labs(y = 'Estimated abundance',
       title = 'Petersen estimates') +
  geom_errorbar(aes(x = season,
                    ymax = modeN + SE, 
                    ymin = modeN - SE),
                color = 'blue',
                size = 1) +
  geom_point(aes(x = season, y = modeN),
             size = 2.5, color = 'blue',
             shape = 15) +
  geom_text(aes(x = season, 
                y = modeN + SE + 2000,
                label = LP.n,
                vjust = 'bottom'),
            color = 'blue',
            fontface = 'bold',
            size = 4) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

p.Tally <- ggplot(data = summary.Tally.stats.df) + 
  labs(y = 'Estimated abundance',
       title = 'Tally counts estimates') +
  geom_errorbar(aes(x = season,
                    ymax = N_med + SE, 
                    ymin = N_med - SE),
                color = 'red',
                size = 1) +
  geom_point(aes(x = season, y = N_med),
             size = 2.5, color = 'red',
             shape = 15) +
  # geom_text(aes(x = season, 
  #               y = N_med + SE + 2000,
  #               label = LP.n,
  #               vjust = 'bottom'),
  #           color = 'blue',
  #           fontface = 'bold',
  #           size = 4) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

p.all <- ggplot() + 
  labs(y = 'Estimated abundance',
       x = "Season") +
  geom_errorbar(data = Bailey_HG, 
                aes(x = season, 
                    ymax = Nhat.Baileys + N.SE.Baileys, 
                    ymin = Nhat.Baileys - N.SE.Baileys),
                size = 1,
                color = 'black') + 
  geom_point(data = Bailey_HG,
             aes(x = season, y = Nhat.Baileys),
             size = 2.5, color = 'black',
             shape = 15) +
  geom_errorbar(data = results.Petersen.df,
                aes(x = season+0.2,
                    ymax = modeN + SE, 
                    ymin = modeN - SE),
                color = 'blue',
                size = 1) +
  geom_point(data = results.Petersen.df,
             aes(x = season+0.2, y = modeN),
             size = 2.5, color = 'blue',
             shape = 15) +
  geom_errorbar(data = summary.Tally.stats.df,
                aes(x = season+0.4,
                    ymax = N_med + SE, 
                    ymin = N_med - SE),
                color = 'red',
                size = 1) +
  geom_point(data = summary.Tally.stats.df,
             aes(x = season+0.4, y = N_med),
             size = 2.5, color = 'red',
             shape = 15)+
  geom_errorbar(data = summary.stats.df,
                aes(x = season+0.6,
                    ymax = N_med_all + N_SE, 
                    ymin = N_med_all - N_SE),
                color = 'gold',
                size = 1) +
  geom_point(data = summary.stats.df,
             aes(x = season+0.6, y = N_med_all),
             size = 2.5, color = 'gold',
             shape = 15)
  

  # theme(plot.title = element_text(hjust = 0.5),
  #       axis.title = element_text(size = 12),
  #       axis.text = element_text(size = 12),
  #       legend.position = c(0.1, 0.6),
  #       legend.title.align = 0.5)

if (saveFig){
  ggsave(filename = paste0('figures/estimated_abundance_',
                           model.v, '.png'),
         dpi = 600,
         plot = p1)
         
}

dat.tally.df <- data.frame(season = N.df$season,
                           N_med = N.df$modeN,
                           tally = apply(tally.list.3$TC, 
                                         MARGIN = 1,
                                         FUN = max,
                                         na.rm = T))

p3 <- ggplot(data = dat.tally.df) + 
  labs(y = 'log(Estimated abundance)', 
       x = 'log(Maximum tally counts)') + 
  geom_point(aes(y = log(N_med),
                 x = log(tally))) + 
  # geom_errorbar(aes(x = N_med,
  #                   ymax = N_high,
  #                   ymin = N_low)) + 
  theme(legend.position = c(0.9, 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

dat.5 <- na.omit(full_join(dat.tally.df, 
                           nest.suc.df, 
                           by = 'season'))


p6 <- ggplot(data = results.all) + 
  geom_point(aes(x = season, 
                 y = N_width, size = LP.n))

# this is an interesting plot:
p7 <- ggplot(data = results.all) + 
  geom_point(aes(x = TC.n, 
                 y = LP.n, size = N_width)) + 
  scale_size_continuous(name = '95% CI width',
                        breaks = c(1000, 5000, 10000, 50000, 80000)) + 
  xlab('Tally count sample size') + 
  ylab('Lincoln-Petersen sample size') + 
  theme(legend.position = c(0.8, 0.8),
        axis.text = element_text(size =12),
        axis.title = element_text(size = 12))



# don't have alpha and beta for nesting success any more
# dat.S.pdf <- data.frame(x = seq(0, 1, by=0.01),
#                         y = dbeta(seq(0, 1, by=0.01), 
#                                   summary.alpha_s[3],
#                                   summary.beta_s[3]))
# 
# p6 <- ggplot(data = dat.S.pdf,
#              aes(x = x, y = y)) +
#   geom_line(size = 2) + 
#   labs(title = 'Nesting success',
#        x = 'Nesting success',
#        y = 'Probability density') + 
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12))
