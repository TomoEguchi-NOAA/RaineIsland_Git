#get_data_tally_fcn.R

# a function to get data for tally count and nesting success data
# to be used in jags for estimating abundance.

get.tally.data.v4 <- function(nestfile = 'data/nesting_success.csv',
                              tallyfile = 'data/Tally.csv'){
  library(tidyverse)
  library(reshape2)
  # R is the average internesting period in days
  #r <- 1; dat <- dat.nest.tally
  extract.dates <- function(dat){
    dat$INCLUDE <- 1
    for (r in 1:nrow(dat)){
      date1 <- dat[r, 'date']
      yr <- year(date1)
      if (date1 < as.Date(paste0(yr, '-11-29')) | date1 > as.Date(paste0(yr, '-12-14'))){
        dat[r, 'INCLUDE'] <- 0
      }
    }	
    return(dat)
  }
  
  # Get data first:
  dat.nest.success <- na.omit(read.csv(file = nestfile,
                                       header = T))
  dat.nest.success$date <- as.Date(dat.nest.success$CDATE, 
                                    format = '%m/%d/%Y')
  dat.nest.success <- extract.dates(dat.nest.success)
  dat.nest.success.raw <- dat.nest.success
  dat.nest.success <- dat.nest.success[dat.nest.success$INCLUDE > 0,] %>% 
    mutate(., NSUC = CLUTCHES/TURTWATER) %>%
    select(., SEASON, TURTWATER, CLUTCHES, NSUC)
  
  # TURTWATER matrix:
  # get sample sizes first:
  dat.nest.success %>% group_by(SEASON) %>% count(., SEASON) %>%
   as.data.frame()-> n.dat.NSUC
  
  TURTWATER.mat <- CLUTCHES.mat <- NSUC.mat <- matrix(nrow = nrow(n.dat.NSUC),
                                                      ncol = max(n.dat.NSUC$n))
  k<-1
  for (k in 1:nrow(n.dat.NSUC)){
    tmp <- filter(dat.nest.success, SEASON == n.dat.NSUC[k, 'SEASON'])
    TURTWATER.mat[k, 1:n.dat.NSUC[k, 'n']] <- tmp$TURTWATER
    CLUTCHES.mat[k, 1:n.dat.NSUC[k, 'n']] <- tmp$CLUTCHES
    NSUC.mat[k, 1:n.dat.NSUC[k, 'n']] <- tmp$NSUC
  }
  
  NSUC.n.vec <- n.dat.NSUC$n
  uniq.season.NSUC <- n.dat.NSUC$SEASON
  
  # Tally counts:
  dat.tally <- na.omit(read.csv(file = tallyfile))
  dat.tally$date <- as.Date(dat.tally$CDATE, 
                            format = '%m/%d/%Y')
  
  # then extract dates - Nov 29 to Dec 14:
  dat.tally <- extract.dates(dat.tally)
  dat.tally <- dat.tally[dat.tally$INCLUDE > 0, 
                         c('SEASON', 'TALLY', 'date')]
  # find sampling dates
  unique.dates.TC <- unique(dat.tally$date)
  
  # there was a duplicate of 2011-12-04 (removed on 10/13/2017)
  
  dat.4 <- data.frame(DATE = unique.dates.TC,
                      TALLY = NA,
                      SEASON = NA)
  
  k <- 1
  options(warn = 2)  # this stops when a warning happens
  for (k in 1:length(unique.dates.TC)){
    tmp <- dat.tally[dat.tally$date == unique.dates.TC[k],]
    dat.4[k, 'TALLY'] <- sum(tmp$TALLY, na.rm=T)
    dat.4[k, 'SEASON'] <- year(unique.dates.TC[k])
    
  }
  
  uniq.yrs <- unique(dat.tally$SEASON)
  
  tally.n.vec <- aggregate(dat.tally,
                           by = list(dat.tally$SEASON),
                           FUN = length)[, c('TALLY')]
  
  TC.mat <- matrix(nrow = length(uniq.yrs),
                   ncol = max(tally.n.vec))
  t <- 1  
  for (t in 1:length(uniq.yrs)){
    dat1 <- dat.4[dat.4$SEASON == uniq.yrs[t],]
    #dt <- as.integer(dat1$DATE - dat1$DATE[1])
    #s_dt <- 1 + sum(cumprod((1-s)))
    TC <- dat1$TALLY
    TC[TC == 0] <- NA
    
    TC.mat[uniq.yrs == uniq.yrs[t], 
           1:tally.n.vec[t]] <- TC
    #s_dt.vec[t] <- s_dt
  }

  TC.df <- as.data.frame(TC.mat)  
  # TC.df$SEASON <- uniq.yrs
  # TC.nest.data <- left_join(TC.df, 
  #                           dat.nest.success, 
  #                           by='SEASON')
  # dat.nest <- TC.nest.data[, c('SEASON', 'TURTWATER', 
  #                              'CLUTCHES', 'NSUC1')]
  
  unique.seasons.TURT <- unique(dat.nest.success.raw$SEASON)
  
  TURT.1 <- arrange(aggregate(TURTWATER ~ SEASON + date, 
                      data = dat.nest.success.raw,
                      FUN = sum),
                    by = SEASON) %>% group_by(date)
  
  n_sectors.1 <- arrange(aggregate(TURTWATER ~ SEASON + date,
                         data = dat.nest.success.raw,
                         FUN = length),
                       by = SEASON) %>% group_by(date)
  
  n_TURT <- aggregate(date ~ SEASON, data = TURT.1, 
                      FUN = length)
  
  n_sectors <- matrix(nrow = length(unique.seasons.TURT),
                      ncol = (max(n_TURT$date) + 1))
  
  n_dates <- matrix(nrow = length(unique.seasons.TURT),
                    ncol = 2)
  k <- 1
  for (k in 1:length(unique.seasons.TURT)){
    tmp.TURT <- filter(TURT.1, SEASON == unique.seasons.TURT[k])
    tmp.sect <- filter(n_sectors.1, SEASON == unique.seasons.TURT[k] )
    #tmp.clutches <- filter(dat.nest.success.raw, SEASON == unique.seasons.TURT[k])
    
    unique.dates.TURT <- unique(tmp.TURT$date)
    n_sectors[k, 1] <- n_dates[k, 1] <-unique.seasons.TURT[k]
    n_dates[k, 2] <- length(unique.dates.TURT)
    #k1 <- 1
    for (k1 in 1:length(unique.dates.TURT)){
      #tmp1.TURT <- filter(tmp.TURT, date == unique.dates.TURT[k1])
      tmp1.sect <- filter(tmp.sect, date == unique.dates.TURT[k1])
      #tmp1.clutches <- filter(tmp.clutches, date == unique.dates.TURT[k1])
      n_sectors[k, k1+1] <- tmp1.sect$TURTWATER
      #TURT.mat[k, k1+1] <- sum(tmp1.TURT$TURTWATER)
      #clutches.mat[k, k1+1] <- sum(tmp1.clutches$CLUTCHES)
    }
  }
  
  return.list.tally <- list(TC = TC.mat,
                            tally.n = tally.n.vec,
                            nest.raw = dat.nest.success.raw,
                            turtwater.mat = TURTWATER.mat,
                            nsuc.mat = NSUC.mat,
                            clutches.mat = CLUTCHES.mat,
                            uniq.yrs = uniq.yrs,
                            dat.4 = dat.4,
                            n_dates = n_dates,
                            n_sectors = n_sectors)
  
  return(return.list.tally)
}


