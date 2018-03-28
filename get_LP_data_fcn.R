#Petersen_jags

get.LP.data <- function(datafile){
  
  datafile <- 'data/Petersen_data.csv'
  
  #library(dplyr)
  # data
  # get the raw data first
  Petersen_cols <- cols(
    CDATE = col_datetime("%m/%d/%Y"),
    SEASON = col_integer(),
    PINDEX = col_double(),
    PINDEXSE = col_double(),
    M = col_integer(),
    NSMALL = col_integer(),
    MSMALL = col_integer(),
    PAINT = col_character())
  
  raw.data <- readr::read_csv(datafile,
                              col_types = Petersen_cols)
  
  # summarize data by season
  # M = the number of marked individuals with particular paint color
  # n = the total number of turtles sighted per sampling
  # m = the number of resighted individuals with particular color 
  raw.data %>% select(., SEASON, M, NSMALL, MSMALL, PAINT) %>%
    mutate(., season = SEASON, M = M, n = NSMALL, 
           m = MSMALL, paint = PAINT) %>%
    group_by(season, paint) %>%
    #summarise(., M = first(M), n = sum(n), m = sum(m)) %>%
    #mutate(., fSeason = as.factor(season)) 
    group_by(., season) %>%
    arrange(., .by_group = TRUE) -> LP_data.1
  
  # mutate returns with warning messages about converting 
  # into characters... but the line below works. 
  LP_data.1$fSeason <- as.factor(LP_data.1$season)
  
  # pool by colors:
  raw.data %>% select(., SEASON, M, NSMALL, MSMALL, PAINT) %>%
    transmute(., season = SEASON, 
              M = M, 
              n = NSMALL, 
              m = MSMALL, paint = PAINT) %>%
    group_by(season, paint) %>%
    summarise(., M = first(M), 
              n = sum(n), m = sum(m)) -> LP_data_paint
  
  # the number of different colors per year
  LP.sample.n <- aggregate(LP_data.1,
                           by = list(LP_data.1$fSeason),
                           FUN = length)%>%
    select(., Group.1, fSeason) %>%
    transmute(., season = Group.1, fSeason = fSeason)
  
  LP.paint.n <- aggregate(LP_data_paint,
                          by = list(as.factor(LP_data_paint$season)),
                          FUN = length) %>%
    select(., Group.1, paint) %>%
    transmute(., season = Group.1, paint = paint)
  
  uniq.yrs <- unique(LP_data.1$season)
  
  LP.nt.vec <- vector(mode = 'numeric', 
                      length = length(uniq.yrs))
  
  n.mat <- M.mat <- m.mat <- matrix(nrow = length(uniq.yrs),
                                    ncol = max(LP.sample.n$fSeason))
  n.paint.mat <- M.paint.mat <- m.paint.mat <- matrix(nrow = length(uniq.yrs),
                                    ncol = max(LP.paint.n$paint))
  t <- 1
  for (t in 1:length(uniq.yrs)){
    if (sum(LP.sample.n$season == uniq.yrs[t]) > 0){
      dat0 <- LP_data.1[LP_data.1$season == uniq.yrs[t],]
      
      LP.nt.vec[t] <- LP.sample.n$fSeason[LP.sample.n$season == uniq.yrs[t]]
      M.mat[t, 1:LP.nt.vec[t]] <- dat0$M
      m.mat[t, 1:LP.nt.vec[t]] <- dat0$m 
      n.mat[t, 1:LP.nt.vec[t]] <- dat0$n
      
      dat1 <- filter(LP_data_paint, season == uniq.yrs[t])
      tmp1 <- filter(LP.paint.n, season == uniq.yrs[t])
      M.paint.mat[t, 1:tmp1$paint] <- dat1$M
      m.paint.mat[t, 1:tmp1$paint] <- dat1$m
      n.paint.mat[t, 1:tmp1$paint] <- dat1$n
    }
  }

  return.list <- list(n.mat = n.mat,
                      m.mat = m.mat,
                      M.mat = M.mat,
                      LP.nt.vec = LP.nt.vec,
                      uniq.yrs = uniq.yrs,
                      LP_data_1 = LP_data.1,
                      LP_sample_n = LP.sample.n,
                      LP_paint = LP_data_paint,
                      LP_paint_n = LP.paint.n,
                      n.paint.mat = n.paint.mat,
                      m.paint.mat = m.paint.mat,
                      M.paint.mat = M.paint.mat)
  return(return.list)
  
}


