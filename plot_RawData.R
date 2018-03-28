#plot_RawData


rm(list=ls())
source('RaineIsl_functions.R')

library(tidyverse)
#library(ggridges)

saveFig <- F

# input data files 
petersen.file <- 'data/Petersen_data.csv'
nest.file <- 'data/nesting_success.csv'
tally.file <- 'data/Tally.csv'

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
n_dates <- select(n_dates, 'n_dates')
n_dates[is.na(n_dates)] <- 0

TC <- tally.list.3$TC
TC_n <- tally.list.3$tally.n

tally.n.df <- data.frame(n = tally.list.3$tally.n,
                         season = tally.list.3$uniq.yrs)

TC.df.long <- na.omit(reshape2::melt(TC.df, id.vars = 'SEASON',
                                     value.name = 'Counts'))

p.TC <- ggplot(data = TC.df.long) + 
  geom_point(aes(x = SEASON, y = Counts),
             size = 2.5,
             shape = 19,
             color = rgb(red = 149, 
                         green = 55, 
                         blue = 53, 
                         maxColorValue = 255)) + 
  labs(x = '', y = 'Counts') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

if (saveFig)
  ggsave(filename = paste0('figures/TallyCounts.png'),
         dpi = 600,
         height = 3.5,
         width = 4.8,
         units = 'in',
         plot = p.TC)

# get Bailey and Hypergeo results here:
Bailey_HG <- BaileyHyperGeo(petersen.file)

p.N.Bailey <- ggplot(data = Bailey_HG) + 
  geom_errorbar(aes(x = season,
                    ymin = Nhat.HG - N.SE.HG,
                    ymax = Nhat.HG + N.SE.HG),
                color = rgb(red = 79,
                            green = 129,
                            blue = 189,
                            maxColorValue = 255),
                size = 1.0) +
  geom_point(aes(x = season, y = Nhat.HG),
             shape = 19,
             size = 2.5,
             color = rgb(red = 79,
                         green = 129,
                         blue = 189,
                         maxColorValue = 255)) + 
  scale_y_continuous(labels = function(n) format(n, scientific = FALSE)) +
  labs(y = 'Mean Bailey\'s estimates (+/- SE)', x = '') + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

if (saveFig)
  ggsave(filename = paste0('figures/BaileysEstimates.png'),
         dpi = 600,
         height = 3.5,
         width = 4.8,
         units = 'in',
         plot = p.N.Bailey)

TC.df.mean <- group_by(TC.df.long, SEASON) %>%
  summarize(., TC.mean = mean(Counts), TC.SE = SE(Counts)) %>%
  transmute(season = SEASON, TC.mean = TC.mean, TC.SE = TC.SE)

Bailey.df <- select(Bailey_HG, Nhat.HG, N.SE.HG, season) %>%
  transmute(., season = season, Nhat = Nhat.HG, Nhat_SE = N.SE.HG)

TC.Bailey <- na.omit(data.frame(left_join(TC.df.mean, 
                                          Bailey.df, by = 'season')))

p.Bailey.TC <- ggplot(TC.Bailey) + 
  geom_point(aes(x = TC.mean, y = Nhat),
             shape = 19, size = 2.5) + 
  geom_errorbar(aes(x = TC.mean, 
                    ymin = Nhat - Nhat_SE,
                    ymax = Nhat + Nhat_SE),
                size = 1.0,
                color = rgb(red = 79,
                            green = 129,
                            blue = 189,
                            maxColorValue = 255)) +
  geom_errorbarh(aes(y = Nhat, x = TC.mean,
                     xmin = TC.mean - TC.SE,
                     xmax = TC.mean + TC.SE),
                 size = 1.0,
                 color = rgb(red = 149, 
                             green = 55, 
                             blue = 53, 
                             maxColorValue = 255)) +
  scale_x_continuous(labels = function(n) format(n, scientific = FALSE)) +
  scale_y_continuous(labels = function(n) format(n, scientific = FALSE)) +
  labs(y = 'Mean Bailey\'s estimates (+/- SE)', x = 'Mean counts (+/- SE)') + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

if (saveFig)
  ggsave(filename = paste0('figures/BaileysVsTally.png'),
         dpi = 600,
         height = 6.8,
         width = 4.8,
         units = 'in',
         plot = p.Bailey.TC)


# nest success part
clutches.df <- data.frame(tally.list.3$clutches) 
colnames(clutches.df) <- c('season', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6')
clutches.df.long <- na.omit(reshape2::melt(clutches.df, id.vars = 'season',
                                   value.name = 'Clutches')) %>%
  group_by(., season)

turtles.df <- data.frame(tally.list.3$TURT)
colnames(turtles.df) <- c('season', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6')
turtles.df.long <- na.omit(reshape2::melt(turtles.df, id.vars = 'season',
                                  value.name = 'Turtles')) %>%
  group_by(., season)

nest.success.mean.df <- left_join(clutches.df.long, 
                             turtles.df.long, by = 'season') %>%
  mutate(., nest_success = Clutches/Turtles) %>%
  select(., season, Clutches, Turtles, nest_success) %>%
  group_by(., season) %>%
  summarise(., mean.nest.suc = mean(nest_success),
            SE.nest.suc = SE(nest_success)) %>%
  data.frame()

nest.success.TC.df <- na.omit(left_join(nest.success.mean.df, TC.Bailey, 
                                by = 'season'))

p.nest.success.TC <- ggplot(data = nest.success.TC.df) + 
  geom_point(aes(x = TC.mean, 
                 y = mean.nest.suc),
             size = 2.5,shape = 19) + 
  geom_errorbarh(aes(y = mean.nest.suc, 
                     x = TC.mean,
                     xmin = TC.mean - TC.SE,
                     xmax = TC.mean + TC.SE),
                 size = 1.0,
                 color = rgb(red = 149, 
                             green = 55, 
                             blue = 53, 
                             maxColorValue = 255)) + 
  geom_errorbar(aes(x = TC.mean, 
                    ymin = mean.nest.suc - SE.nest.suc,
                    ymax = mean.nest.suc + SE.nest.suc),
                size = 1.0,
                color = rgb(red = 119,
                            green = 147,
                            blue = 60,
                            maxColorValue = 255)) +
    scale_x_continuous(labels = function(n) format(n, scientific = FALSE)) +
    scale_y_continuous(labels = function(n) format(n, scientific = FALSE)) +
    labs(y = 'Mean nest success (+/- SE)', x = 'Mean counts (+/- SE)') + 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
  
if (saveFig)
  ggsave(filename = paste0('figures/NestSuccessVsTC.png'),
         dpi = 600,
         height = 3.5,
         width = 4.8,
         units = 'in',
         plot = p.nest.success.TC)

p.nest.success.Bailey <- ggplot(data = nest.success.TC.df) + 
  geom_point(aes(x = Nhat, 
                 y = mean.nest.suc),
             size = 2.5,shape = 19) + 
  geom_errorbarh(aes(y = mean.nest.suc, 
                     x = Nhat,
                     xmin = Nhat - Nhat_SE,
                     xmax = Nhat + Nhat_SE),
                 size = 1.0,
                 color = rgb(red = 79, 
                             green = 129, 
                             blue = 189, 
                             maxColorValue = 255)) + 
  geom_errorbar(aes(x = Nhat, 
                    ymin = mean.nest.suc - SE.nest.suc,
                    ymax = mean.nest.suc + SE.nest.suc),
                size = 1.0,
                color = rgb(red = 119,
                            green = 147,
                            blue = 60,
                            maxColorValue = 255)) +
  scale_x_continuous(labels = function(n) format(n, scientific = FALSE)) +
  scale_y_continuous(labels = function(n) format(n, scientific = FALSE)) +
  labs(y = 'Mean nest success (+/- SE)', x = 'Bailey\'s estimate (+/- SE)') + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

if (saveFig)
  ggsave(filename = paste0('figures/NestSuccessVsBailey.png'),
         dpi = 600,
         height = 3.5,
         width = 4.8,
         units = 'in',
         plot = p.nest.success.Bailey)
