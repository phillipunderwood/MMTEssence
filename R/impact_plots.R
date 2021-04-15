GetPeriodIndices <- function(period, times) {
  # Computes indices belonging to a period in data.
  #
  # Args:
  #   period:  two-element vector specifying start and end of a period, having
  #            the same data type as `times. The range from `period[1]` to
  #            `period[2]` must have an intersect with `times`.
  #   times:   vector of time points; can be of integer or of POSIXct type.
  #
  # Returns:
  #   A two-element vector with the indices of the period start and end within
  #   `times`.

  # Check input
  assert_that(length(period) == 2)
  assert_that(!anyNA(times))
  assert_that(identical(class(period), class(times)) ||
                (is.numeric(period) && is.numeric(times)))
  # Check if period boundaries are in the right order, and if `period` has an
  # overlap with `times`.
  assert_that(period[1] <= period[2])
  assert_that(period[1] <= tail(times, 1), period[2] >= times[1])

  # Look up values of start and end of period in `times`; also works if the
  # period start and end time are not exactly present in the time series.
  indices <- seq_along(times)
  is.period <- (period[1] <= times) & (times <= period[2])
  # Make sure the period does match any time points.
  assert_that(any(is.period),
              msg = "The period must cover at least one data point")
  period.indices <- range(indices[is.period])
  return(period.indices)
}

CreateDataFrameForPlot <- function(impact) {
  # Creates a long-format data frame for CreateImpactPlot().
  #
  # Args:
  #   impact: \code{CausalImpact} results object
  #
  # Returns:
  #   data frame of: time, response, mean, lower, upper, metric

  # Check input
  assert_that((class(impact) == "CausalImpact"))
  assert_that(!isTRUE(all(is.na(impact$series[, -c(1, 2)]))),
              msg = "inference was aborted; cannot create plot")

  # Create data frame from zoo series
  data <- as.data.frame(impact$series)
  data <- cbind(time = time(impact$series), data)

  # Reshape data frame
  tmp1 <- data[, c("time", "response", "point.pred", "point.pred.lower",
                   "point.pred.upper"), drop = FALSE]
  names(tmp1) <- c("time", "response", "mean", "lower", "upper")
  tmp1$baseline <- NA
  tmp1$metric <- "original"
  tmp2 <- data[, c("time", "response", "point.effect", "point.effect.lower",
                   "point.effect.upper"), drop = FALSE]
  names(tmp2) <- c("time", "response", "mean", "lower", "upper")
  tmp2$baseline <- 0
  tmp2$metric <- "pointwise"
  tmp2$response <- NA
  tmp3 <- data[, c("time", "response", "cum.effect", "cum.effect.lower",
                   "cum.effect.upper"), drop = FALSE]
  names(tmp3) <- c("time", "response", "mean", "lower", "upper")
  tmp3$metric <- "cumulative"
  tmp3$baseline <- 0
  tmp3$response <- NA
  data <- rbind(tmp1, tmp2, tmp3)
  data$metric <- factor(data$metric, c("original", "pointwise", "cumulative"))
  rownames(data) <- NULL
  return(data)
}

CreatePeriodMarkers <- function(pre.period, post.period, times) {
  # Creates a vector of period markers to display.
  #
  # Args:
  #   pre.period:  vector of 2 time points that define the pre-period.
  #   post.period: vector of 2 time points that define the post-period.
  #   times:       vector of time points.
  #
  # Returns:
  #   Vector of period markers that should be displayed, generally depicting the
  #   first and last time points of pre- and post-period. The start of the pre-
  #   period is not shown if it coincides with the first time point of the time
  #   series; similarly, the last time point of the post-period is not shown if
  #   it coincides with the last time point of the series. If there is no gap
  #   between pre- and post-period, the start marker of the post-period is
  #   omitted.

  pre.period.indices <- GetPeriodIndices(pre.period, times)
  post.period.indices <- GetPeriodIndices(post.period, times)
  markers <- NULL
  if (pre.period.indices[1] > 1) {
    markers <- c(markers, times[pre.period.indices[1]])
  }
  markers <- c(markers, times[pre.period.indices[2]])
  if (pre.period.indices[2] < post.period.indices[1] - 1) {
    markers <- c(markers, times[post.period.indices[1]])
  }
  if (post.period.indices[2] < length(times)) {
    markers <- c(markers, times[post.period.indices[2]])
  }
  markers <- as.numeric(markers)
  return(markers)
}

# Tell R CMD check to treat columns of data frames used in `ggplot` functions
# as global variables; this avoids false positives of "no visible binding for
# global variable ..." during the check.
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("baseline", "lower", "response", "upper"))
}

CreateImpactPlot <- function(impact, metrics = c("original", "pointwise",
                                                 "cumulative"), period_include = c("pre", "post")) {
  # Creates a plot of observed data and counterfactual predictions.
  #
  # Args:
  #   impact:  \code{CausalImpact} results object returned by
  #            \code{CausalImpact()}.
  #   metrics: Which metrics to include in the plot. Can be any combination of
  #            "original", "pointwise", and "cumulative".
  #
  # Returns:
  #   A ggplot2 object that can be plotted using plot().

  # Create data frame of: time, response, mean, lower, upper, metric
  data <- CreateDataFrameForPlot(impact)

  # Select metrics to display (and their order)
  assert_that(is.vector(metrics))
  metrics <- match.arg(metrics, several.ok = TRUE)
  data <- data[data$metric %in% metrics, , drop = FALSE]
  data$metric <- factor(data$metric, metrics)

  if(period_include == "pre") {
    data <- data[which(data[,1] <= impact$model$pre.period[2]),,drop = F]
  } else if(period_include == "post"){
    data <- data[which(data[,1] >= impact$model$post.period[1]),,drop = F]
  } else if (period_include == "all") {
    data <- data[which(data[,1] < impact$model$post.period[2]),,drop = F]
  }

  # Initialize plot
  q <- ggplot(data, aes(x = time, group = 1)) + theme_bw(base_size = 30)
  q <- q + xlab("") + ylab("")
  if (length(metrics) > 1) {
    q <- q + facet_grid(metric ~ ., scales = "free_y")
  }

  # Add prediction intervals
  q <- q + geom_ribbon(aes(ymin = lower, ymax = upper),
                       data, fill = "slategray2")

  # Add pre-period markers
  xintercept <- CreatePeriodMarkers(impact$model$pre.period,
                                    impact$model$post.period,
                                    time(impact$series))
  q <- q + geom_vline(xintercept = xintercept,
                      colour = "darkgrey", size = 4, linetype = "dashed")

  # Add zero line to pointwise and cumulative plot
  q <- q + geom_line(aes(y = baseline),
                     colour = "darkgrey", size = 3, linetype = "solid",
                     na.rm = TRUE)

  # Add point predictions
  q <- q + geom_line(aes(y = mean), data,
                     size = 2, colour = "darkblue", linetype = "dashed",
                     na.rm = TRUE)

  # Add observed data
  q <- q + geom_line(aes(y = response), size = 2,  na.rm = TRUE)
  return(q)
}

plot.CausalImpact <- function(x, ...) {
  # Creates a plot of observed data and counterfactual predictions.
  #
  # Args:
  #   x:   A \code{CausalImpact} results object, as returned by
  #        \code{CausalImpact()}.
  #   ...: Can be used to specify \code{metrics}, which determines which panels
  #        to include in the plot. The argument \code{metrics} can be any
  #        combination of "original", "pointwise", "cumulative". Partial matches
  #        are allowed.
  #
  # Returns:
  #   A ggplot2 object that can be plotted using plot().
  #
  # Examples:
  #   \dontrun{
  #   impact <- CausalImpact(...)
  #
  #   # Default plot:
  #   plot(impact)
  #
  #   # Customized plot:
  #   impact.plot <- plot(impact) + ylab("Sales")
  #   plot(impact.plot)
  #   }

  return(CreateImpactPlot(x, ...))
}
