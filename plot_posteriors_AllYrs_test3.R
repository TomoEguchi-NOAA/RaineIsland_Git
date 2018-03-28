#plot_Petersen_TAlly_results


rm(list=ls())
source('RaineIsl_functions.R')
library(bayesplot)
library(ggridges)
library(tidyverse)
library(rjags)
#library(reshape2)
saveFig <- F

model.v.I <- 'v15'
model.v.T <- 'v7'
model.v.P <- 'v2'

run.date.int <- '2018-01-25'
run.date.P <- '2018-01-12'
run.date.Tally <- '2018-01-25'

petersen.file <- 'data/Petersen_data.csv'
nest.file = 'data/nesting_success.csv'
tally.file = 'data/Tally.csv'

# jags output files are here: _P designated parallel runs:
petersen.outfile <- paste0('RData/jags_out_Petersen_', 
                           model.v.P,'_', run.date.P, '_P.RData')
integrated.outfile <- paste0('RData/jags_out_Integrated_', model.v.I, '_', 
                             run.date.int,  '_P_test3.RData')
tally.outfile <- paste0('RData/jags_out_Tally_', model.v.T, '_',
                        run.date.Tally, '_P_test3.RData')

tally.list.3 <- get.tally.data.v3(nestfile = nest.file,
                                  tallyfile = tally.file)
LP.list <- get.LP.data(petersen.file)

nest.suc.df <- na.omit(data.frame(season = tally.list.3$nest.raw$SEASON,
                                  NSUC1 = tally.list.3$nest.raw$NSUC1))
load(integrated.outfile)

N.names <- c("N[1]", "N[2]", "N[3]", "N[4]", "N[5]", "N[6]", "N[7]", "N[8]", "N[9]", 
             "N[10]", "N[11]", "N[12]", "N[13]", "N[14]", "N[15]", "N[16]", "N[17]", 
             "N[18]", "N[19]", "N[20]", "N[21]", "N[22]", "N[23]", "N[24]")

s.names <- c("s[1]", "s[2]", "s[3]", "s[4]", "s[5]", "s[6]", "s[7]", "s[8]", "s[9]", 
             "s[10]", "s[11]", "s[12]", "s[13]", "s[14]", "s[15]", "s[16]", "s[17]", 
             "s[18]", "s[19]", "s[20]", "s[21]", "s[22]", "s[23]", "s[24]")

q.names <- c("q[1]", "q[2]", "q[3]", "q[4]", 
             "q[5]", "q[6]", "q[7]", "q[8]", "q[9]", 
             "q[10]", "q[11]", "q[12]", "q[13]", 
             "q[14]", "q[15]", "q[16]", "q[17]", 
             "q[18]", "q[19]", "q[20]", "q[21]", 
             "q[22]", "q[23]", "q[24]")

uniq.yrs <- LP.list$uniq.yrs
valuecol <- "sample"
keycol <- 'index'

df.integrated <- do.call(rbind.data.frame, tmp.I)
sum.tmp.I <- summary(tmp.I)
N.df.I <- data.frame(sum.tmp.I$quantiles[grep(pattern = 'N[/[]',
                                          row.names(sum.tmp.I$quantiles)), 
                                     c('2.5%', '50%', '97.5%')])
colnames(N.df.I) <- c('lowN', 'modeN', 'highN')

df.N.integrated <- select(df.integrated, starts_with('N')) %>%
  gather_(key_col = keycol, value = valuecol, N.names) %>%
  mutate(season = rep(uniq.yrs, each = nrow(df.integrated))) %>%
  mutate(f.season = as.factor(season))

df.s.integrated <- select(df.integrated, starts_with('s[')) %>%
  gather_(key_col = keycol, value = valuecol, s.names) %>%
  mutate(season = rep(uniq.yrs, each = nrow(df.integrated))) %>%
  mutate(f.season = as.factor(season))

df.q.integrated <- select(df.integrated, starts_with('q[')) 
df.s_alpha.integrated <- select(df.integrated, 's_alpha')
df.s_beta.integrated <- select(df.integrated, 's_beta')

p.N.integrated <- ggplot(data = df.N.integrated, 
                       aes(x = sample, y = f.season)) + 
  geom_density_ridges2(fill = 'red') + 
  scale_y_discrete(limits = rev(levels(df.N.integrated$f.season))) + 
  #scale_x_continuous(limits = c(0, 100000)) + 
  labs(x = 'Abundance (Integrated)', y = '') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

yr <- 2001
p.N.yr <- ggplot(data = filter(df.N.integrated, 
                                            season == yr)) + 
  geom_density(aes(x = sample),
               fill = 'red', color = 'black') +
  labs(x = 'Abundance (Integrated)', y = '')  + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)
  
if (saveFig){
  ggsave(filename = paste0('figures/estimated_abundance_integrated.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.N.integrated)
  
  ggsave(filename = paste0('figures/estimated_abundance_integrated_2001.png'),
         dpi = 600,
         height = 3,
         width = 4.8,
         units = 'in',
         plot = p.N.integrated.2001)

  ggsave(filename = paste0('figures/estimated_abundance_integrated_2007.png'),
         dpi = 600,
         height = 3,
         width = 4.8,
         units = 'in',
         plot = p.N.integrated.2007)
  
}

p.s.integrated <- ggplot(data = df.s.integrated, 
                         aes(x = sample, y = f.season)) + 
  geom_density_ridges2(fill = 'red') + 
  scale_y_discrete(limits = rev(levels(df.s.integrated$f.season))) + 
  #scale_x_continuous(limits = c(0, 100000)) + 
  labs(x = 'Nesting success probability (Integrated)', y = '') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

# get just Petersen estimates:
LP.list <- get.LP.data(petersen.file)

load(petersen.outfile)

N.names.Petersen <- c("N[1]", "N[2]", "N[3]", "N[4]", "N[5]", 
                      "N[6]", "N[7]", "N[8]", "N[9]", "N[10]", 
                      "N[11]", "N[12]", "N[13]", "N[14]", "N[15]", 
                      "N[16]", "N[17]", "N[18]", "N[19]", "N[20]", 
                      "N[21]", "N[22]", "N[23]", "N[24]")
df.Petersen <- do.call(rbind.data.frame, tmp.P)
sum.tmp.P <- summary(tmp.P)
df.N.Petersen <- select(df.Petersen, starts_with('N')) %>%
  gather_(key_col = keycol, value = valuecol, N.names.Petersen) %>%
  mutate(season = rep(LP.list$uniq.yrs, each = nrow(df.Petersen))) %>%
  mutate(f.season = as.factor(season))

p.Petersen <- ggplot(data = df.N.Petersen, 
                     aes(x = sample, y = f.season)) + 
  geom_density_ridges2(fill = 'blue') + 
  scale_y_discrete(limits = rev(levels(df.N.Petersen$f.season))) + 
  #scale_x_continuous(limits = c(0, 200000)) + 
  labs(x = 'Abundance (Bayesian binomial)', y = '') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

yr <- 2001
p.N.Petersen <- ggplot(data = filter(df.N.Petersen,
                                     season == 2001)) + 
  geom_density(aes(x = sample), 
               fill = 'blue', color = 'black') +  
  labs (x = 'Abundance (Petersen)', y = '')  + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

if (saveFig){
  ggsave(filename = paste0('figures/estimated_abundance_BayesPetersen.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.Petersen)
  
}
# Tally estimates:
load(tally.outfile)

sum.tmp.T <- summary(tmp.T)
df.Tally <- do.call(rbind.data.frame, tmp.T)
df.N.Tally <- select(df.Tally, starts_with('N')) %>%
  gather_(key_col = keycol, value = valuecol, N.names) %>%
  mutate(season = rep(uniq.yrs, each = nrow(df.Tally))) %>%
  mutate(f.season = as.factor(season))

df.q.Tally <- select(df.Tally, 'q')
df.s_alpha.Tally <- select(df.Tally, 's_alpha')
df.s_beta.Tally <- select(df.Tally, 's_beta')

df.s.Tally <- select(df.Tally, starts_with('s[')) %>%
  gather_(key_col = keycol, value = valuecol, s.names) %>%
  mutate(season = rep(uniq.yrs, each = nrow(df.Tally))) %>%
  mutate(f.season = as.factor(season))

p.N.Tally <- ggplot(data = df.N.Tally, 
                  aes(x = sample, y = f.season)) + 
  geom_density_ridges2(fill = 'gold') +  
  scale_y_discrete(limits = rev(levels(df.N.Tally$f.season))) + 
  #scale_x_continuous(limits = c(0, 100000)) + 
  labs(x = 'Abundance (Counts)', y = '') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

yr <- 2001
p.N.yr <- ggplot(data = filter(df.N.Tally,
                                  season == yr)) + 
  labs (x = 'Abundance (Counts)', y = '',
        title = yr)  + 
  geom_density(aes(x = sample), 
               fill = 'gold', color = 'black') + 
  scale_x_continuous(labels=function(n){format(n, scientific = FALSE)}) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

yr <- 2000
p.N <- ggplot() + 
  geom_density(data = filter(df.N.integrated, 
                             season == yr), 
               aes(x = sample), 
               fill = 'red', color = 'black', alpha = 0.5) + 
  geom_density(data = filter(df.N.Petersen,
                              season == yr),
               aes(x = sample), 
               fill = 'blue', color = 'black', alpha = 0.5) + 
  geom_density(data = filter(df.N.Tally,
                              season == yr), 
               aes(x = sample), 
               fill = 'gold', color = 'black', alpha = 0.5) + 
  
  labs(x = 'Abundance', y = '',
       title = yr)  + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)
yr <- 2000
p.s <- ggplot() + 
  geom_density(data = filter(df.s.integrated, 
                             season == yr), 
               aes(x = sample), 
               fill = 'red', color = 'black', alpha = 0.5) + 
  geom_density(data = filter(df.s.Tally,
                             season == yr), 
               aes(x = sample), 
               fill = 'gold', color = 'black', alpha = 0.5) + 
  
  labs(x = 's', y = '',
       title = yr)  + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

plot.var1 <- function(varname){
  p1 <- ggplot() + 
    geom_density(data = select(df.integrated, varname), 
                 aes(x = get(varname, df.integrated)), 
                 fill = 'red', color = 'black', alpha = 0.5) + 
    geom_density(data = select(df.Tally, varname), 
                 aes(x = get(varname, df.Tally)), 
                 fill = 'gold', color = 'black', alpha = 0.5) + 
    geom_density(data = filter(df.N.Petersen,
                               season == yr),
                 aes(x = sample), 
                 fill = 'blue', color = 'black', alpha = 0.5) + 
    labs(x = varname, y = '')  + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.y = element_blank(),
          legend.position = c(0.1, 0.6),
          legend.title.align = 0.5)
  return(p1)
}

if (saveFig){
  ggsave(filename = paste0('figures/estimated_abundance_Tally_',
                           model.v.Tally, '.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.Tally)
  
}

# get Bailey and Hypergeo results here:
Bailey_HG <- BaileyHyperGeo('data/Petersen_data.csv')

# s.df <- as.data.frame(summary.s)
# s.df <- select(s.df, -starts_with('parameter'))
# s.df$season <- N.df$season
# colnames(s.df) <- c('Slow', 'S25', 'Smode',
#                     'S75', 'Shigh', 'season')
#tally.df <- cbind(TC.df, s.df)
Bailey_HG$f.season <- as.factor(Bailey_HG$season)

# the function inside scale_x_continuous was found here:
# https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot

p.Bailey <- ggplot(data = Bailey_HG) + 
  labs(x = 'Abundance (Bailey\'s)',
       y = '') +
  geom_errorbarh(aes(y = f.season, x = Nhat.Baileys,
                    xmax = Nhat.Baileys + N.SE.Baileys, 
                    xmin = Nhat.Baileys - N.SE.Baileys),
                size = 1,
                color = 'black') + 
  geom_point(aes(y = f.season, x = Nhat.Baileys),
             size = 2.5, color = 'black',
             shape = 15) + 

  # geom_errorbarh(aes(x = Nhat.Baileys.w, y = f.season, 
  #                   xmax = Nhat.Baileys.w + N.SE.Baileys.w, 
  #                   xmin = Nhat.Baileys.w - N.SE.Baileys.w),
  #               size = 1,
  #               color = 'red', alpha = 0.5) + 
  # geom_point(aes(y = f.season, x = Nhat.Baileys.w),
  #            size = 2.5, color = 'red',
  #            shape = 15) +

  geom_text(aes(y = f.season, 
                x = Nhat.Baileys + N.SE.Baileys + 2000,
                label = n,
                vjust = 'bottom'),
            color = 'black',
            fontface = 'bold',
            size = 4) +
  scale_y_discrete(limits = rev(levels(Bailey_HG$f.season))) + 
  scale_x_continuous(labels=function(n){format(n, scientific = FALSE)}) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

if (saveFig){
  ggsave(filename = paste0('figures/estimated_abundance_Bailey.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.Bailey)
  
}

  # theme(plot.title = element_text(hjust = 0.5),
  #       axis.title = element_text(size = 12),
  #       axis.text = element_text(size = 12),
  #       legend.position = c(0.1, 0.6),
  #       legend.title.align = 0.5)

