#RaineIsl_functions



sysInfo <- Sys.info()
ifelse(sysInfo[1] == 'Linux',
       source('~/Documents/R/tools/TomosFunctions.R'),
       source('~/R/tools/TomosFunctions.R'))

library(tidyverse)
library(lubridate)
# library(readr)
# library(dplyr)
# library(ggplot2)

source('get_LP_data_fcn.R')
source('get_tally_data_fcn.R')
source('get_tally_data_fcn_v2.R')
source('get_tally_data_fcn_v3.R')
source('get_tally_data_fcn_v4.R')
source('N_estimation_BaileyHyperGeo.R')

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
