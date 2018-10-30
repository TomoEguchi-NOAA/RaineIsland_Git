#get_data_tally_fcn.R

# a function to get data for tally count and nesting success data
# to be used in jags for estimating abundance.

get.tally.data <- function(R=12, 
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
  dat.nest.success$dates <- as.Date(dat.nest.success$CDATE, 
                                    format = '%m/%d/%Y')
  
  dat.tally <- na.omit(read.csv(file = tallyfile))
  #dat.tally$CDATE <- as.factor(dat.tally$CDATE)
  dat.tally$dates <- as.Date(dat.tally$CDATE, 
                             format = '%m/%d/%Y')
  
  # combine the two datasets:
  dat.nest.tally <- full_join(dat.nest.success, 
                              dat.tally, 
                              by = 'dates')
  
  colnames(dat.nest.tally) <- c('CDATE', 'SEASON', 
                                'TURTWATER',
                                'CLUTCHES', 'NSUC1', 
                                'SECTOR', 
                                'date',
                                'CDATE.1', 'SEASON.1', 
                                'TALLY')
  
  dat.nest.tally$season <- year(dat.nest.tally$date)
  
  dat.nest.tally <- arrange(select(.data = dat.nest.tally, 
                                   CDATE, SEASON, 
                                   TURTWATER, CLUTCHES, NSUC1,
                                   SECTOR, date, CDATE.1, 
                                   SEASON.1, TALLY, season), 
                            -desc(date))
  
  # look at mean turtwater and mean tally counts
  # for using a regression analysis - not used in the end...
  dat.nest.tally$NSUC2 <- dat.nest.tally$CLUTCHES/dat.nest.tally$TURTWATER
  dat.tmp <- aggregate(dat.nest.tally[, c('TALLY', 
                                          'TURTWATER')],
                       by = list(dat.nest.tally$season),
                       FUN = mean, na.rm = T)
  
  #plot(log(dat.tmp$TALLY), log(dat.tmp$TURTWATER))
  lm.fit.1 <- lm(log(TALLY) ~ log(TURTWATER), 
                 data = dat.tmp)
  
  # then extract dates - Nov 29 to Dec 14:
  dat.2 <- extract.dates(dat.nest.tally)
  dat.3 <- dat.2[dat.2$INCLUDE > 0,]
  
  # find sampling dates
  unique.dates <- unique(dat.3$date)
  
  # prepare a data frame 
  dat.4 <- data.frame(DATE = unique.dates,
                      TURTWATER = NA,
                      CLUTCHES = NA,
                      SUC1 = NA,
                      TALLY = NA,
                      season = NA)
  
  k <- 126
  options(warn = 2)  # this stops when a warning happens
  for (k in 1:length(unique.dates)){
    tmp <- dat.3[dat.3$date == unique.dates[k],]
    dat.4[k, 'TURTWATER']  <- sum(tmp$TURTWATER, na.rm = T)
    dat.4[k, 'CLUTCHES'] <- sum(tmp$CLUTCHES, na.rm = T)
    dat.4[k, 'NSUC1'] <- sum(tmp$CLUTCHES, na.rm = T)/sum(tmp$TURTWATER, na.rm = T)
    dat.4[k, 'TALLY'] <- sum(tmp$TALLY, na.rm=T)
    dat.4[k, 'season'] <- year(unique.dates[k])
  }
  
  mean.NSUC1 <- mean(dat.4$NSUC1, na.rm = T)
  
  # sort by DATE in the ascending order
  dat.4 <- arrange(select(.data = dat.4, 
                          DATE, TURTWATER,CLUTCHES,NSUC1,
                          TALLY, season), 
                   -desc(DATE))
  
  uniq.yrs <- unique(dat.4$season)
  
  tally.n <- aggregate(dat.4,
                       by = list(dat.4$season),
                       FUN = length)[, c('Group.1', 
                                         'season')]
  
  tally.n.vec <- vector(mode = 'numeric',
                        length = length(uniq.yrs))
  TURT.mat <- nest.mat <- TC.mat <- matrix(nrow = length(uniq.yrs),
                                           ncol = max(tally.n$season))
  
  count.tally.mean <- count.tally.pred <-vector(mode = 'numeric', 
                                                length = length(uniq.yrs))
  t <- 20
  for (t in 1:length(uniq.yrs)){
    dat1 <- dat.4[dat.4$season == uniq.yrs[t],]
    # fill in the missing values:
    if (any(dat1[, 'TALLY'] == 0)){
      # if there are more than one Tally counts, then 
      # take the mean of those to fill in the blanks
      if (length(dat1[,'TALLY']>0) > 1){
        dat1[dat1$TALLY == 0, 'TALLY'] <- mean(dat1$TALLY[dat1$TALLY > 0])
        count.tally.mean[t] <- count.tally.mean[t] + 1
      } else {
        # if just one or less tally count, then
        # use the regression to predict
        pred <- predict(lm.fit.1,
                        newdata = data.frame(TURTWATER = dat1[dat1$TALLY == 0, 'TURTWATER']))
        dat1[dat1$TALLY == 0, 'TALLY'] <- exp(pred)
        count.tally.pred[t] <- count.tally.pred[t] + 1
      }
    }
    
    tally.n.vec[uniq.yrs == uniq.yrs[t]] <- tally.n$season[tally.n$Group.1 == uniq.yrs[t]]
    
    #dt <- as.integer(dat1$DATE - dat1$DATE[1])
    #s_dt <- 1 + sum(cumprod((1-s)))
    TC <- dat1$TALLY
    TC[TC == 0] <- NA
    
    TURT.mat[uniq.yrs == uniq.yrs[t], 
             1:tally.n.vec[t]] <- dat1$TURTWATER
    
    nest.mat[uniq.yrs == uniq.yrs[t], 
             1:tally.n.vec[t]] <- dat1$CLUTCHES
    TC.mat[uniq.yrs == uniq.yrs[t], 
           1:tally.n.vec[t]] <- TC
    #s_dt.vec[t] <- s_dt
  }
  tally.impute <- data.frame(mean = count.tally.mean,
                             pred = count.tally.pred,
                             season = uniq.yrs)
  
  return.list.tally <- list(TC = TC,
                            TURT.mat = TURT.mat,
                            nest.mat = nest.mat,
                            TC.mat = TC.mat,
                            tally.impute = tally.impute,
                            tally.n.vec = tally.n.vec,
                            dat.4 = dat.4,
                            lm.fit = lm.fit.1,
                            uniq.yrs = uniq.yrs)
  
  return(return.list.tally)
}


