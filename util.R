utilFormat <- function(prefix, data.vec, width, format, flag, sep){
  return(
    paste(prefix, 
          formatC(data.vec, 
                  width = width, 
                  format = format, 
                  flag = flag
          ), 
          sep = sep
    )
  )
}



# Classify a list of timeseries using TWDTW
#
# @param weight.fun       A weight function
# @param timeseries.list  A list of zoo objects
# @param                  A twdtwTimeSeries object
# @return                 A list of twdtwMatches objects
classifyTwdtw <- function(weight.fun, timeseries.list, patterns_ts){
  match.list <- lapply(timeseries.list, 
                       function(ts.zoo, patterns_ts, fun_weight){
                         ts <- twdtwTimeSeries(ts.zoo, labels="Time series")
                         matches <- twdtwApply(x = ts, y = patterns_ts, 
                                               weight.fun = fun_weight, keep=TRUE
                         )
                         return(matches)
                       }, 
                       patterns_ts = patterns_ts, 
                       fun_weight = weight.fun
  )
  return(match.list)
}



# extends a TWDTW pattern into nyears
#
# @param ts.zoo A zoo timeseries object
# @param nyears An integer. The number of years to extend
# @return       A zoo timeseries
extendTS <- function(ts.zoo, nyears){
  p <- as.numeric(difftime(index(ts.zoo)[2], index(ts.zoo)[1], units = "days")) # time inbetween observations
  yvec <- rep(NA, times = 365/p)                                                # one year of NA observations
  yvec[1:nrow(coredata(ts.zoo))] <- coredata(ts.zoo)[, 1]                       # replace the actual observations in the vector
  tsvec <- rep(yvec, times = (nyears + 1))                                       # repeat the observations nyears
  seqDate <- seq.Date(from = index(ts.zoo)[1], by = p, length.out = length(yvec)) # sequence of one year of dates
  seqDate <- lapply(as.character(seqDate), function(d, nyears){                    # build sequence of dates
    md <- substr(d, 5, nchar(d))
    y <- as.numeric(substr(d, 1, 4))
    as.Date(paste(y + 0:nyears, md, sep = ""))
  }, 
  nyears = nyears
  )
  seqDate <- as.Date(as.vector(do.call("rbind", seqDate)))
  tsmat <- as.matrix(tsvec)
  colnames(tsmat) <- colnames(coredata(ts.zoo))
  return(zoo(x = tsmat, order.by = seqDate))
}



# Filter the inner attributes of zoo time series
#
# @param zoo.list A list of zoo objects
# @param att      A character. The name of the attribute (column) to keep
# @param return   A list of zoo objects
filterTS <- function(zoo.list, att){
  res <- lapply(zoo.list, function(z, att){
    if(is.zoo(z) == FALSE){cat("Not a zoo object"); return(NA)}
    data.fil <- as.matrix(coredata(z)[, att])
    colnames(data.fil) <- att
    return(zoo(x = data.fil, order.by = index(z)))
  }, att = att)
  names(res) <- names(zoo.list)
  return(res)
}


# Count the number of matches found in a TWDTW classification
#
# @param match.list A list of dtwSat objects resulting from twdtwApply.
# @return           A matrix. The number of matches found
matchmatrix <- function(match.list){
  match.mat <- matrix(NA, ncol = length(match.list), nrow = length(match.list))
  colnames(match.mat) <- do.call("rbind", strsplit(names(match.list), split = "-"))[,1]
  rownames(match.mat) <- paste("TS", names(match.list), sep = " ")
  for(i in 1:length(match.list)){
    ts.name <- names(match.list)[i]         # pattern repeated to build the time series
    m <- match.list[[i]]                    # twdtwMatches object. Result of applying TWDTW
    for(j in 1:length(m@alignments[[1]])){
      n <- m@alignments[[1]][[j]]
      ts.match <- as.vector(n$label)        # name of the pattern found
      ts.match.times <- length(n$matching)  # number of times the pattern was found. Not necessarily the best match
      match.mat[i, j] <- ts.match.times
    }
  }
  return(match.mat)
}



# measurement of quality of TWDTW patterns.
#
# @param match.list A list of dtwSat objects resulting from twdtwApply.
# @return           A data.frame. 
qualitymatrix <- function(match.list){
  # TODO: This functions fails with too much noise. Perhaps the lack of matches make the function fail
  ts.name_match_metric <- data.frame()
  for(i in 1:length(match.list)){
    ts.name <- names(match.list)[i]                                             # name of pattern repeated to build the time series
    m <- match.list[[i]]                                                        # twdtwMatches object. Result of applying TWDTW
    ts.match_metric <- data.frame()
    twdtw.dist <- NA
    for(j in 1:length(m@alignments[[1]])){
      n <- m@alignments[[1]][[j]]                                               # alignment to one of the timeseries
      twdtw.dist <- n$distance                                                  # TWDTW distance
      ts.match <- as.vector(n$label)                                            # name of the pattern found
      len.metric <- rep(NA, times = length(n$matching))
      len.pattern <- NA
      len.alignment <- NA
      len.metric <- NA
      if(length(n$matching) > 0){
        for(k in 1:length(n$matching)){
          len.pattern <- length(n$matching[[k]][[1]])
          len.alignment <- n$matching[[k]][[2]][length(n$matching[[k]][[2]])] - n$matching[[k]][[2]][1] + 1
          len.metric[k] <- len.alignment / len.pattern         # length metric
        }
      }
      ts.match_metric <- rbind(ts.match_metric, as.data.frame(cbind(ts.match, len.metric, twdtw.dist)))
    }
    ts.name_match_metric <- rbind(ts.name_match_metric, cbind(ts.name, ts.match_metric))
  }
  #ts.name_match_metric <- as.numeric(ts.name_match_metric$len.metric)
  return(ts.name_match_metric)
}



# compute the duration metric
#
# @param match.list A list of dtwSat objects resulting from twdtwApply.
# @return           A list of 2 data.frames. The first with the duration metric and the second with the TWDTW distance
computeMetrics <- function(match.list){
  # get the metrics
  qlen.df <- qualitymatrix(match.list) # get the duration metric in time steps
  qlen.df$len.metric <- as.numeric(levels(qlen.df$len.metric))[qlen.df$len.metric] # R sucks!
  qlen.df$twdtw.dist <- as.numeric(levels(qlen.df$twdtw.dist))[qlen.df$twdtw.dist] # R sucks!
  names(qlen.df) <- c("time_series", "match", "duration_metric", "twdtw_dist")
  # compute statistics
  qlen.df.durmean <- aggregate(duration_metric ~ time_series * match, data = qlen.df, FUN = mean)
  qlen.df.dursd  <- aggregate(duration_metric ~ time_series * match, data = qlen.df, FUN = sd)
  qlen.df.durmin <- aggregate(duration_metric ~ time_series * match, data = qlen.df, FUN = min)
  qlen.df.durmax <- aggregate(duration_metric ~ time_series * match, data = qlen.df, FUN = max)
  qlen.df.dursum <- aggregate(duration_metric ~ time_series * match, data = qlen.df, FUN = sum)
  qlen.df.count   <- aggregate(duration_metric ~ time_series * match, data = qlen.df, FUN = length)
  # -
  qlen.df.dismean <- aggregate(twdtw_dist ~ time_series * match, data = qlen.df, FUN = mean)
  qlen.df.dissd   <- aggregate(twdtw_dist ~ time_series * match, data = qlen.df, FUN = sd)
  qlen.df.dismin  <- aggregate(twdtw_dist ~ time_series * match, data = qlen.df, FUN = min)
  qlen.df.dismax  <- aggregate(twdtw_dist ~ time_series * match, data = qlen.df, FUN = max)
  qlen.df.dissum  <- aggregate(twdtw_dist ~ time_series * match, data = qlen.df, FUN = sum)
  # change column names
  colnames(qlen.df.durmean) [colnames(qlen.df.durmean)  == "duration_metric"] <- "dm.mean"
  colnames(qlen.df.dursd)   [colnames(qlen.df.dursd)    == "duration_metric"] <- "dm.sd"
  colnames(qlen.df.durmin)  [colnames(qlen.df.durmin)   == "duration_metric"] <- "dm.min"
  colnames(qlen.df.durmax)  [colnames(qlen.df.durmax)   == "duration_metric"] <- "dm.max"
  colnames(qlen.df.dursum)  [colnames(qlen.df.dursum)   == "duration_metric"] <- "dm.sum"
  colnames(qlen.df.count)   [colnames(qlen.df.count)    == "duration_metric"] <- "alignments"
  # -
  colnames(qlen.df.dismean) [colnames(qlen.df.dismean)  == "twdtw_dist"] <- "td.mean"
  colnames(qlen.df.dissd)   [colnames(qlen.df.dissd)    == "twdtw_dist"] <- "td.sd"
  colnames(qlen.df.dismin)  [colnames(qlen.df.dismin)   == "twdtw_dist"] <- "td.min"
  colnames(qlen.df.dismax)  [colnames(qlen.df.dismax)   == "twdtw_dist"] <- "td.max"
  colnames(qlen.df.dissum)  [colnames(qlen.df.dissum)   == "twdtw_dist"] <- "td.sum"
  # join
  durmetric.df <- qlen.df.count
  durmetric.df["dm.mean"]  <- qlen.df.durmean$dm.mean
  durmetric.df["dm.sd"]    <- qlen.df.dursd$dm.sd
  durmetric.df["dm.min"]   <- qlen.df.durmin$dm.min
  durmetric.df["dm.max"]   <- qlen.df.durmax$dm.max
  durmetric.df["dm.sum"]   <- qlen.df.dursum$dm.sum
  durmetric.df <- durmetric.df[order(durmetric.df$time_series, durmetric.df$match), ] # sort
  # -
  dismetric.df <- qlen.df.count
  dismetric.df["td.mean"]  <- qlen.df.dismean$td.mean
  dismetric.df["td.sd"]    <- qlen.df.dissd$td.sd
  dismetric.df["td.min"]   <- qlen.df.dismin$td.min
  dismetric.df["td.max"]   <- qlen.df.dismax$td.max
  dismetric.df["td.sum"]   <- qlen.df.dissum$td.sum
  dismetric.df <- dismetric.df[order(dismetric.df$time_series, dismetric.df$match), ] # sort
  return(list(durmetric.df, dismetric.df))
}