#get_data_tally_fcn.R

# a function to get data for tally count and nesting success data
# to be used in jags for estimating abundance.

get.tally.data.v2 <- function(R=12, 
                              nestfile = 'data/nesting_success.csv',
                              tallyfile = 'data/Tally.csv'){
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
  dat.nest.success <- dat.nest.success[dat.nest.success$INCLUDE > 0,]
  dat.nest.success.raw <- dat.nest.success
  
  # compute the mean nest success rate per season
  dat.nest.success <- aggregate(dat.nest.success[,c('TURTWATER',
                                                    'CLUTCHES')],
                                by = list(dat.nest.success$SEASON),
                                FUN = sum, na.rm = T)
  dat.nest.success$NSUC <- dat.nest.success$CLUTCHES/dat.nest.success$TURTWATER
  colnames(dat.nest.success) <- c('season', 'TURTWATER', 'CLUTCHES', 'NSUC1')
  
  dat.tally <- na.omit(read.csv(file = tallyfile))
  dat.tally$date <- as.Date(dat.tally$CDATE, 
                            format = '%m/%d/%Y')
  
  # then extract dates - Nov 29 to Dec 14:
  dat.tally <- extract.dates(dat.tally)
  dat.tally <- dat.tally[dat.tally$INCLUDE > 0, 
                         c('SEASON', 'TALLY', 'date')]
  # find sampling dates
  unique.dates <- unique(dat.tally$date)
  
  dat.4 <- data.frame(DATE = unique.dates,
                      TALLY = NA,
                      SEASON = NA,
                      NSUC1 = NA)
  k <- 1
  options(warn = 2)  # this stops when a warning happens
  for (k in 1:length(unique.dates)){
    tmp <- dat.tally[dat.tally$date == unique.dates[k],]
    dat.4[k, 'TALLY'] <- sum(tmp$TALLY, na.rm=T)
    dat.4[k, 'SEASON'] <- year(unique.dates[k])

  }
  
  uniq.yrs <- unique(dat.tally$SEASON)
  tally.n.vec <- aggregate(dat.tally,
                           by = list(dat.tally$SEASON),
                           FUN = length)[, c('TALLY')]
  
  TURT.mat <- TC.mat <- matrix(nrow = length(uniq.yrs),
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
  TC.df$season <- uniq.yrs
  TC.nest.data <- left_join(TC.df, 
                            dat.nest.success, 
                            by='season')
  dat.nest <- TC.nest.data[, c('season', 'TURTWATER', 
                               'CLUTCHES', 'NSUC1')]
  return.list.tally <- list(TC = TC.mat,
                            tally.n = tally.n.vec,
                            nest = dat.nest,
                            nest.raw = dat.nest.success.raw,
                            uniq.yrs = uniq.yrs,
                            dat.4 = dat.4)
  
  return(return.list.tally)
}


