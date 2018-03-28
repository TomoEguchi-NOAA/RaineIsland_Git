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
model.v.T <- 'v7'  # v8 and v9 did not converge
model.v.P <- 'v2'
run.date <- '2018-01-25'  # 1/31/2018 analyses used data pooled by paint colors rather than
# treating them all different for Petersen estimates. All previous analyses were based on
# per sampling. 

# input data files 
petersen.file <- 'data/Petersen_data.csv'
nest.file <- 'data/nesting_success.csv'
tally.file <- 'data/Tally.csv'

counts.color <- c(192, 80, 77)

# jags output files are here:
petersen.outfile <- paste0('RData/jags_out_Petersen_', 
                           model.v.P,'_', run.date, '_P.RData')
integrated.outfile <- paste0('RData/jags_out_Integrated_', model.v.I, '_', 
                             run.date,  '_P.RData')
tally.outfile <- paste0('RData/jags_out_Tally_', model.v.T, '_',
                        run.date, '_P.RData')

change.parameter.name <- function(df.in, conv.table){
  col.names <- vector(mode = 'character', length = ncol(df.in))
  for (k in 1:ncol(df.in)){
    colname1 <- colnames(df.in)[k]
    col.names[k] <-filter(conv.table, 
                          parameter == colname1)$season
  }
  colnames(df.in) <- col.names
  return(df.in)
}

tally.list.3 <- get.tally.data.v3(nestfile = nest.file,
                                  tallyfile = tally.file)

nest.suc.df <- na.omit(data.frame(season = tally.list.3$nest.raw$SEASON,
                                  NSUC1 = tally.list.3$nest.raw$NSUC1))
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

seq.1 <- seq(from = 1, to = 40)
N.names <- paste0('N[', seq.1, ']')
s.names <- paste0('s[', seq.1, ']')

uniq.yrs <- tally.list.3$uniq.yrs
valuecol <- "sample"
keycol <- 'index'

df.integrated <- do.call(rbind.data.frame, tmp.I)
sum.tmp.I <- summary(tmp.I)
N.df.I <- data.frame(sum.tmp.I$quantiles[grep(pattern = 'N[/[]',
                                          row.names(sum.tmp.I$quantiles)), 
                                     c('2.5%', '50%', '97.5%')])
colnames(N.df.I) <- c('lowN', 'modeN', 'highN')

df.N.integrated <- select(df.integrated, starts_with('N')) %>%
  change.parameter.name(., params2seasons.N) %>%
  gather(1974:1984, factor_key = TRUE) 

colnames(df.N.integrated) <- c("season", "sample") 

df.N.integrated <- mutate(df.N.integrated, 
                          f.season = as.factor(season))

df.N.integrated.yr <- group_by(df.N.integrated, 
                               season) %>%
  summarise(., median = median(sample),
            mean = mean(sample),
            SE = SE(sample))

df.s.integrated <- select(df.integrated, starts_with('s[')) %>%
  gather_(key_col = keycol, value = valuecol, s.names) %>%
  mutate(season = rep(uniq.yrs, each = nrow(df.integrated))) %>%
  mutate(f.season = as.factor(season))

if (model.v.I != 'v15' & model.v.I != 'v16' ){
  df.q.integrated <- select(df.integrated, 'q') 
  
} else {
  q.names <- paste0('q[', seq.1, ']')
  
  df.q.integrated <- select(df.integrated, starts_with('q[')) %>%
    gather_(key_col = keycol, value = valuecol, q.names) %>%
    mutate(season = rep(uniq.yrs, each = nrow(df.integrated))) %>%
    mutate(f.season = as.factor(season))
}
df.s_alpha.integrated <- select(df.integrated, 's_alpha')
df.s_beta.integrated <- select(df.integrated, 's_beta')

p.N.integrated <- ggplot(data = df.N.integrated, 
                       aes(x = sample, y = f.season)) + 
  geom_density_ridges2(fill = rgb(red = 196,
                                  green = 189,
                                  blue = 151,
                                  maxColorValue = 255)) + 
  scale_y_discrete(limits = rev(levels(df.N.integrated$f.season))) + 
  #coord_flip() + 
  #scale_x_continuous(limits = c(0, 100000)) + 
  labs(x = 'Abundance (Integrated)', y = '') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

if (saveFig)
  ggsave(filename = paste0('figures/estimated_abundance_integrated_',
                           model.v.I, '.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.N.integrated)

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

if (model.v.I == 'v15' | model.v.I == 'v16' ){
  
  p.q.integrated <- ggplot(data = df.q.integrated, 
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
  
}

#yr <- 2001
plot.N.integrated.1yr <- function(yr){
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
  return(p.N.yr)
}

yr <- 1989
plot.s.integrated.1yr <- function(yr){
  p.s.1yr <- ggplot() + 
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
  return(p.s.1yr)
}

# get just Petersen estimates:
datafile <- 'data/Petersen_data.csv'
LP.list <- get.LP.data(datafile)

load(petersen.outfile)
N.names.2 <- seq(from = 1, to = 24)
N.names.Petersen <- paste0("N[", N.names.2, "]")

df.Petersen <- do.call(rbind.data.frame, tmp.P)
sum.tmp.P <- summary(tmp.P)
df.N.Petersen <- select(df.Petersen, starts_with('N')) %>%
  gather_(key_col = keycol, value = valuecol, N.names.Petersen) %>%
  mutate(season = rep(LP.list$uniq.yrs, each = nrow(df.Petersen))) %>%
  mutate(f.season = as.factor(season))

p.Petersen <- ggplot() + 
  geom_density_ridges2(data = df.N.Petersen, 
                       aes(x = sample, y = f.season),
                       fill = rgb(red = 79, 
                                  green = 129, 
                                  blue = 189, 
                                  maxColorValue = 255)) + 
  scale_y_discrete(limits = rev(levels(df.N.Petersen$f.season))) + 
  geom_text(data = data.frame(n = LP.list$LP.nt.vec,
                              f.season = unique(df.N.Petersen$f.season)),
            aes(y = f.season, x = 200000,
                label = n, vjust = 'bottom'),
            color = 'black', fontface = 'bold',
            size = 6) +
  #scale_x_continuous(limits = c(0, 200000)) + 
  labs(x = 'Abundance (capture-recapture)', y = '') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

if (saveFig)
  ggsave(filename = paste0('figures/estimated_abundance_BayesPetersen.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.Petersen)


yr <- 2001
p.N.Petersen <- function(yr){
  p.N <- ggplot(data = filter(df.N.Petersen,
                              season == yr)) + 
  geom_density(aes(x = sample), 
               fill = 'blue', color = 'black') +  
  labs (x = 'Abundance (Petersen)', y = '')  + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)
  return(p.N)
}

# Tally estimates:
load(tally.outfile)
sum.tmp.T <- summary(tmp.T)
df.Tally <- do.call(rbind.data.frame, tmp.T)
df.N.Tally <- select(df.Tally, starts_with('N')) %>%
  gather_(key_col = keycol, value = valuecol, N.names) %>%
  mutate(season = rep(uniq.yrs, each = nrow(df.Tally))) %>%
  mutate(f.season = as.factor(season))

if (model.v.T != 'v8' & model.v.T != 'v9' ){
  df.q.Tally <- select(df.Tally, 'q')
  
} else {
  df.q.Tally <- select(df.Tally, starts_with('q[')) %>%
    gather_(key_col = keycol, value = valuecol, q.names) %>%
    mutate(season = rep(uniq.yrs, each = nrow(df.Tally))) %>%
    mutate(f.season = as.factor(season))
}

df.s_alpha.Tally <- select(df.Tally, 's_alpha')
df.s_beta.Tally <- select(df.Tally, 's_beta')

df.s.Tally <- select(df.Tally, starts_with('s[')) %>%
  gather_(key_col = keycol, value = valuecol, s.names) %>%
  mutate(season = rep(uniq.yrs, each = nrow(df.Tally))) %>%
  mutate(f.season = as.factor(season))

p.N.Tally <- ggplot(data = df.N.Tally, 
                  aes(x = sample, y = f.season)) + 
  geom_density_ridges2(fill = rgb(red = 149,
                                  green = 55,
                                  blue = 53,
                                  maxColorValue = 255)) +  
  scale_y_discrete(limits = rev(levels(df.N.Tally$f.season))) + 
  #scale_x_continuous(limits = c(0, 100000)) + 
  labs(x = 'Abundance (counts)', y = '') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

if (saveFig)
  ggsave(filename = paste0('figures/estimated_abundance_Tally_',
                           model.v.T, '.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.N.Tally)

p.s.Tally <- ggplot(data = df.s.Tally, 
                         aes(x = sample, y = f.season)) + 
  geom_density_ridges2(fill = 'gold') + 
  scale_y_discrete(limits = rev(levels(df.s.integrated$f.season))) + 
  #scale_x_continuous(limits = c(0, 100000)) + 
  labs(x = 'Nesting success probability (Counts)', y = '') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

p.N.Tally.1yr <- function(yr){
  p.N.Tally <- ggplot(data = filter(df.N.Tally,
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
  return(p.N.Tally)
}

if (model.v.T == 'v9' | model.v.T == 'v8'){
  
  p.q.Tally <- ggplot(data = df.q.Tally, 
                      aes(x = sample, y = f.season)) + 
    geom_density_ridges2(fill = 'gold') + 
    scale_y_discrete(limits = rev(levels(df.q.Tally$f.season))) + 
    #scale_x_continuous(limits = c(0, 100000)) + 
    labs(x = 'Nesting success probability (Integrated)', y = '') + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.position = c(0.1, 0.6),
          legend.title.align = 0.5)
  
}

plot.var1 <- function(varname){
  p1 <- ggplot() + 
    geom_density(data = select(df.integrated, varname), 
                 aes(x = get(varname, df.integrated)), 
                 fill = 'red', color = 'black', alpha = 0.5) + 
    geom_density(data = select(df.Tally, varname), 
                 aes(x = get(varname, df.Tally)), 
                 fill = 'gold', color = 'black', alpha = 0.5) + 
    labs(x = varname, y = '')  + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.y = element_blank(),
          legend.position = c(0.1, 0.6),
          legend.title.align = 0.5)
  return(p1)
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
                    xmax = Nhat.Baileys + N.SE.Baileys * 2, 
                    xmin = Nhat.Baileys - N.SE.Baileys * 2),
                size = 1,
                color = 'black') + 
  geom_point(aes(y = f.season, x = Nhat.Baileys),
             size = 2.5, color = 'black',
             shape = 15) + 
  geom_text(aes(y = f.season, 
                x = Nhat.Baileys + N.SE.Baileys * 2 + 2000,
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

# Create random samples using the estimate and SE - normal distributions:
Bailey.stats <- select(Bailey_HG, Nhat.Baileys, N.SE.Baileys)
Bailey.samples <- as.data.frame(apply(Bailey.stats, MARGIN = 1, 
                        FUN = function(x) rnorm(n = 1500,
                                                mean = x[1],
                                                sd = x[2]))) 
colnames(Bailey.samples) <- Bailey_HG$f.season
Bailey.samples <- reshape2::melt(Bailey.samples)
colnames(Bailey.samples) <- c('season', 'sample')
Bailey.samples$f.season <- as.factor(Bailey.samples$season)

p.N.Bailey <- ggplot(data = Bailey.samples, 
                     aes(x = sample, y = f.season)) + 
  geom_density_ridges2(fill = rgb(red = 119,
                                  green = 147,
                                  blue = 60,
                                  maxColorValue = 255)) + 
  scale_y_discrete(limits = rev(levels(Bailey.samples$f.season))) + 
  #scale_x_continuous(limits = c(0, 200000)) + 
  labs(x = 'Abundance (Bailey\'s)', y = '') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.1, 0.6),
        legend.title.align = 0.5)

if (saveFig)
  ggsave(filename = paste0('figures/estimated_abundance_Bailey.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.N.Bailey)

plot.N.all.1yr <- function(yr, saveFig = F, x.limits = NA){
  # generate normal random numbers for Bailey's estimate:
  Bailey_data <- filter(Bailey_HG, season == yr)
  if (nrow(Bailey_data) > 0){
    sample.Bailey <- data.frame(sample = rnorm(n = 1500, 
                                               mean = Bailey_data$Nhat.Baileys,
                                               sd = Bailey_data$N.SE.Baileys))
  } else {
    sample.Bailey <- data.frame(sample = numeric(length = 0))
  }
  
  p.N <- ggplot() + 
    # geom_density(data = sample.Bailey,
    #              aes(x = sample),
    #              fill = 'green',
    #              color = 'black', alpha = 0.5) +
    geom_density(data = filter(df.N.Petersen,
                               season == yr),
                 aes(x = sample), 
                 fill = rgb(red = 79,
                            green = 129,
                            blue = 189,
                            maxColorValue = 255), 
                 color = 'black', alpha = 0.5) + 
    # geom_density(data = filter(df.N.Tally,
    #                            season == yr), 
    #              aes(x = sample), 
    #              fill = rgb(red = 149,
    #                         green = 55,
    #                         blue = 53,
    #                         maxColorValue = 255), 
    #              color = 'black', alpha = 0.5) + 
    geom_density(data = filter(df.N.integrated, 
                               season == yr), 
                 aes(x = sample), 
                 fill = rgb(red = 196,
                            green = 189,
                            blue = 151,
                            maxColorValue = 255), 
                 color = 'black', alpha = 0.7) + 
    labs(x = 'Abundance', y = '',
         title = yr)  + 
    scale_x_continuous(limits = x.limits) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.y = element_blank(),
          legend.position = c(0.1, 0.6),
          legend.title.align = 0.5)
  if (saveFig)
    ggsave(filename = paste0('figures/estimated_abundance_PetersenVsIntegrated_',
                             yr, '.png'),
           dpi = 600,
           height = 3.5, width = 4.8,
           units = 'in', plot = p.N)
  return(p.N)
}


if (saveFig){
  ggsave(filename = paste0('figures/estimated_abundance_BayesPetersen.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.Petersen)
  
  ggsave(filename = paste0('figures/estimated_abundance_Tally_',
                           model.v.T, '.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.N.Tally)
  
  ggsave(filename = paste0('figures/estimated_abundance_Bailey.png'),
         dpi = 600,
         height = 7,
         units = 'in',
         plot = p.Bailey)
  
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

  # theme(plot.title = element_text(hjust = 0.5),
  #       axis.title = element_text(size = 12),
  #       axis.text = element_text(size = 12),
  #       legend.position = c(0.1, 0.6),
  #       legend.title.align = 0.5)

# geom_errorbarh(aes(x = Nhat.Baileys.w, y = f.season, 
#                   xmax = Nhat.Baileys.w + N.SE.Baileys.w, 
#                   xmin = Nhat.Baileys.w - N.SE.Baileys.w),
#               size = 1,
#               color = 'red', alpha = 0.5) + 
# geom_point(aes(y = f.season, x = Nhat.Baileys.w),
#            size = 2.5, color = 'red',
#            shape = 15) +
