## ----chunk_preliminar, include=FALSE-------------------------------------
#-------------------------------------------------------------------------------
# preliminar
#-------------------------------------------------------------------------------
require(dtwSat)
require(knitr)
source("util.R")
att <- "evi"                                                                    # attribute used to filter the given time-series
knit_hooks$set(purl = hook_purl) # Save chunks to Rscript file  

## ----chunk_plotPatterns, echo=FALSE--------------------------------------
#-------------------------------------------------------------------------------
# plot sample patterns
#-------------------------------------------------------------------------------
#patterns.list <- yearly_patterns_mt                                            # switch to the other example set
patterns_ts <-  twdtwTimeSeries(filterTS(patterns.list, att))                   # get timeseries of a single attribute
plot(patterns_ts, type = "patterns")

## ----chunk_plotSampleTS, echo=FALSE--------------------------------------
#-------------------------------------------------------------------------------
# plot sample target timeseries
#-------------------------------------------------------------------------------
ts <- twdtwTimeSeries(filterTS(list(example_ts), att), labels="Time series")
plot(ts, type = "timeseries") + annotate(
  geom = "text", x = example_ts_labels$from+90, 
  y = 0.98, label = example_ts_labels$label, size = 2
)

## ----chunk_sampleClassification, echo=FALSE, message=FALSE, warning=FALSE----
#-------------------------------------------------------------------------------
# classify the target TS using the sample patterns
#-------------------------------------------------------------------------------
matches <- twdtwApply(x = ts, y = patterns_ts, 
                      weight.fun = logisticWeight(alpha = -0.1, beta = 100), 
                      keep=TRUE
)
slotNames(matches)
show(matches)
print(plot(matches, type="alignments", attr = att, threshold = 3.0))

## ----chunk_buildTSpatterns, echo=FALSE-----------------------------------
#-------------------------------------------------------------------------------
# build timeseries by repeating the patterns
#-------------------------------------------------------------------------------
nyears <- 10
patt.list <- filterTS(zoo.list = patterns.list, att = att)                      # attribute-filtered patterns
timeseries.list <- lapply(patt.list, extendTS, nyears = nyears)                 # build repeated time series
timeseries.list <- lapply(timeseries.list, na.approx)                           # replace NA with linear interpolation
pname.list <- as.list(names(timeseries.list))
names(timeseries.list) <- paste(names(timeseries.list), "-ts", sep = "")

## ----chunk_plotSampleTSpatterns, echo=FALSE------------------------------
#-------------------------------------------------------------------------------
# plot a sample pattern and its timeseries
#-------------------------------------------------------------------------------
plot(patt.list[[3]], main = paste(pname.list[[3]], "pattern", sep = " "), ylab = "value")
plot(timeseries.list[[3]], main = paste("Pattern timeseries", "-", names(timeseries.list)[3], sep = " "), ylab = "value")

## ----cunk_plotWeightFunctions, echo=FALSE--------------------------------
#-------------------------------------------------------------------------------
# TWDTW weight functions
#-------------------------------------------------------------------------------
lf.alfa <- -0.1                                                                 # alpha parameter of the logistic weight function
lf.beta <- 100                                                                  # beta parameter of the logistic weight function. It is half the time to get th highest weight using the logistic function
log_weight <- logisticWeight(alpha = lf.alfa, beta = lf.beta)                   # logistic function
lin_weight <- linearWeight(a = 1/(lf.beta * 2), b = 0)                          # linear function
#-------------------------------------------------------------------------------
# plot TWDTW weight functions
#-------------------------------------------------------------------------------
ndays <- 0:(lf.beta * 2)                                                        # number of days to plot
log_res <- log_weight(psi = ndays)
lin_res <- lin_weight(psi = ndays)
plot(x = ndays, y = log_res, type = "l", 
     xlab = "Time difference (days)", ylab = "Weight",
     main = "TWDTW weight functions"
)
lines(x = ndays, y = lin_res, type = "l", lty = 2)
legend('topleft', c("Logistic function", "Linear function"), 
       lty = c(1, 2), bty='n'
)

## ----chunk_twdtwPatternTSclassification, echo=FALSE----------------------
#-------------------------------------------------------------------------------
# pattern-timeseries classification
#-------------------------------------------------------------------------------
match.list <- lapply(timeseries.list, 
                     function(ts.zoo, patterns_ts, fun_weight){
                       ts <- twdtwTimeSeries(ts.zoo, labels="Time series")
                       matches <- twdtwApply(x = ts, y = patterns_ts, 
                                             weight.fun = fun_weight, keep=TRUE
                       )
                       return(matches)
                     }, 
                     patterns_ts = patterns_ts, 
                     fun_weight = log_weight
)
#-------------------------------------------------------------------------------
# plot the results of the classification
#-------------------------------------------------------------------------------
for(i in 1:length(names(timeseries.list))){
  print(
    plot(match.list[[i]], type="alignments", threshold = 20.0, attr = "evi") + 
      ggtitle(paste("Baseline for", names(timeseries.list)[i], sep = " "))
  )
}
#-------------------------------------------------------------------------------
# show the summary of the classification
#-------------------------------------------------------------------------------
match.mat <- matchmatrix(match.list = match.list)
kable(match.mat)

## ----chunk_computeMetrics, echo=FALSE------------------------------------
#-------------------------------------------------------------------------------
# compute metrics
#-------------------------------------------------------------------------------
metrics.list <- computeMetrics(match.list)
durmetric.df <- metrics.list[[1]]
dismetric.df <- metrics.list[[2]]

## ----cunk_displayDuration, echo=FALSE------------------------------------
#-------------------------------------------------------------------------------
# show duration metric
#-------------------------------------------------------------------------------
kable(durmetric.df, digits = 2)

## ----cunk_displayTWDTWdistance, echo=FALSE-------------------------------
#-------------------------------------------------------------------------------
# show TWDTW distance metric 
#-------------------------------------------------------------------------------
kable(dismetric.df, digits = 2)

## ----cunk_indexComputationDTW, echo=FALSE--------------------------------
#-------------------------------------------------------------------------------
# compute index
#-------------------------------------------------------------------------------
durmetric.df["dm.index"] <- durmetric.df$dm.sum/nyears
index.df <- durmetric.df[, c("time_series", "match", "dm.index")]
index.df["td.mean"] <- dismetric.df$td.mean
kable(index.df, digits = 2)

## ----chunk_twdtwPatternTSclassificationNOISE, echo=FALSE-----------------
#-------------------------------------------------------------------------------
# add noise to time-series
# NOTE: too much noise breaks the function qualitymatrix
#-------------------------------------------------------------------------------
sd.vec <- c(0.01, 0.05, 0.1, 0.15) # NOTE: keep it belo 0.2
chosenTS <- 3                                                                   # chosen time series to test
ts.z <- timeseries.list[[chosenTS]]
match.list <- list()
index.list <- list()
for(j in 1:length(sd.vec)){
  # plot
  ts.zoo <- zoo(x = coredata(ts.z) + rnorm(
    n = nrow(coredata(ts.z)), 
    mean = 0, sd = sd.vec[j]
  ), 
  order.by = index(ts.z)
  )
  ts <- twdtwTimeSeries(ts.zoo, labels="Time series")
  matches <- twdtwApply(x = ts, y = patterns_ts, weight.fun = log_weight, keep=TRUE)
  print(
    plot(matches, type="alignments", attr = "evi", threshold = 3.0) + 
      ggtitle(paste(names(timeseries.list)[chosenTS], "Noise SD = ", sd.vec[j], sep = " "))
  )
  # compute metric
  match.list[[names(timeseries.list)[chosenTS]]] <- matches
  metrics.list <- computeMetrics(match.list)
  durmetric.df <- metrics.list[[1]]
  dismetric.df <- metrics.list[[2]]
  durmetric.df["dm.index"] <- durmetric.df$dm.sum/nyears
  index.df <- durmetric.df[, c("time_series", "match", "alignments", "dm.index")]
  index.df["td.mean"] <- dismetric.df$td.mean
  index.list[[j]] <- index.df
  print(kable(index.df, digits = 2, caption = paste(names(timeseries.list)[chosenTS], "Noise SD = ", sd.vec[j], sep = " ")))
}


