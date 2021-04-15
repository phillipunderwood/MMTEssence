# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

####import base packages

#' my title
#'
#' @name CausalImpactBackTest
#'
#' @param data a data frame
#' @param metric conversion metric
#' @param pairs a list of market pairs
#' @param pre_period vector of two dates representing start and end of pre period
#' @param post_period vector of two dates representing start and end of post period
#' @param p_value_thresh the minumum p-value acceptable for valid test
#' @param plots Boolean for plot inclusion
#'
#' @return list of valid tests
#' @importFrom CausalImpact CausalImpact
#' @importFrom utils combn tail
#' @importFrom stats time
#' @importFrom zoo zoo
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 aes ggplot geom_ribbon theme_bw xlab ylab facet_grid geom_line geom_vline
#' @export
causal_impact_test <- function(data, metric, pairs, pre_period, post_period, p_value_thresh, plots=F) {

  results <- vector("list", length = length(pairs))

  if(plots==T) {
    dir.create(paste0(getwd(),"/Test Plots"))
    ###with plots
    for(j in 1:length(pairs)) {
      print(paste0("Backtesting | Pair ", j, "/", length(pairs), " | ", pairs[[j]][1], "-", pairs[[j]][2]))
      combs <- combn(pairs[[j]],2)
      valid <- TRUE
      fails <- NA
      error <- 0
      for (i in 1:ncol(combs)) {
        x1 <- as.data.frame(data[which(data[,"DMA"] == combs[1,i]),c("Date", metric)])
        y <- as.data.frame(data[which(data[,"DMA"] == combs[2,i]),c("Date", metric)])
        time.points <- x1[,1]

        dat <- zoo(cbind(x1[,2,drop=F], y[,2,drop=F]), time.points)
        colnames(dat) <- c('y','x')
        pre.period <- pre_period
        post.period <- post_period

        set.seed(500)
        impact <- CausalImpact(dat, pre.period, post.period)
        error <- impact$summary$p[1]


        if(impact$summary$p[1] <  p_value_thresh) {
          valid <- FALSE
          fails <- combs[,i]
        }
        if(valid == FALSE) {
          i <- ncol(combs) + 1
        }
        if(valid == TRUE & i == ncol(combs)){
          for(k in 1:ncol(combs)){
            png(paste0(getwd(),'/Test Plots/',combs[1,k], '_', combs[2,k], '.png'), width=480*3*1.5,height=480*2*1.5)
            plot(impact, c("original", "pointwise"), period_include = 'pre')
            x <- dev.off()
          }
        }
      }
      results[[j]] <- list(Markets=pairs[[j]],valid=valid,error=error)
    }

  } else {
    ###without plots
    for(j in 1:length(pairs)) {
      print(paste0("Backtesting | Pair ", j, "/", length(pairs), " | ", pairs[[j]][1], "-", pairs[[j]][2]))
      combs <- combn(pairs[[j]],2)
      valid <- TRUE
      fails <- NA
      error <- 0
      for (i in 1:ncol(combs)) {
        x1 <- as.data.frame(data[which(data[,"DMA"] == combs[1,i]),c("Date", metric)])
        y <- as.data.frame(data[which(data[,"DMA"] == combs[2,i]),c("Date", metric)])
        time.points <- x1[,1]

        dat <- zoo(cbind(x1[,2,drop=F], y[,2,drop=F]), time.points)
        colnames(dat) <- c('y','x')
        pre.period <- pre_period
        post.period <- post_period

        set.seed(500)
        impact <- CausalImpact(dat, pre.period, post.period)
        error <- impact$summary$p[1]

        if(impact$summary$p[1] <  p_value_thresh) {
          valid <- FALSE
          fails <- combs[,i]
        }
        if(valid == FALSE) {
          i <- ncol(combs) + 1
        }

      }
      results[[j]] <- list(Markets=pairs[[j]],valid=valid,error=error)
    }
  }
  return(results)
}

#' my title
#'
#' @name GetMarketMatches
#'
#' @param data a data frame
#' @param metric conversion metric
#' @param pre_period vector of two dates representing start and end of pre period
#' @param post_period vector of two dates representing start and end of post period
#' @param corr_thresh the minimum correlation value for acceptable matches
#' @param p_value_thresh the minumum p-value acceptable for valid test
#' @param label a label for your MMT e.g. "YouTube TV Q1 2021"
#' @param rank the amount of matches for each DMA to look at, default is all
#' @param plots Boolean for plot inclusion
#'
#' @return list of valid matches
#' @importFrom CausalImpact CausalImpact
#' @importFrom utils combn write.csv tail
#' @importFrom zoo zoo
#' @importFrom MarketMatching best_matches
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @importFrom stats time
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 aes ggplot geom_ribbon theme_bw xlab ylab facet_grid geom_line geom_vline
#' @export
get_market_matches <- function(data, metric, pre_period, post_period, corr_thresh, p_value_thresh, label=NULL, rank=20, plots=F) {

  if(!(metric%in%colnames(data))) stop('Defined metric not in data set')
  if(!("DMA"%in%colnames(data))) stop('DMA column not found')
  if(!("Date"%in%colnames(data))) stop('Date column not found')
  if(!(inherits(data[,"Date"],'Date'))) stop('Date column in wrong format')
  if(corr_thresh>1||corr_thresh<0) stop('Correlation threshold must be between 0 and 1')
  if(p_value_thresh>1||p_value_thresh<0) stop('P-value threshold must be between 0 and 1')

  dir.create(paste0(getwd(), "/", label))
  setwd(paste0(getwd(), "/", label))

  print("Finding Matches")
  #output dtw relative distance and correlation scores for all regions
  mm <- best_matches(data = data,
                     id_variable = "DMA",
                     date_variable = "Date",
                     matching_variable = metric,
                     parallel = TRUE,
                     warping_limit = 1,
                     dtw_emphasis = 1, #setting dtw_emphasis to 1 ensures scores are output- no impact to chosen control
                     matches = rank, # retrieve scores for all markets for each market
                     start_match_period = pre_period[1],
                     end_match_period = pre_period[2])

  matches <- mm$BestMatches[which(mm$BestMatches$Correlation>=corr_thresh),]
  pairs <- vector('list', nrow(matches))
  for(i in 1:length(pairs)) {pairs[[i]] <- c(matches$DMA[i],matches$BestControl[i])}

  back_test <- (causal_impact_test(data=data,metric=metric,pairs=pairs,pre_period=pre_period,
                                   post_period=post_period,p_value_thresh=p_value_thresh,plots=plots))

  results <- as.data.frame(matrix(NA,length(back_test),8))
  for(i in 1:nrow(results)) {
    results[i,1] <- back_test[[i]]$Markets[1]
    results[i,2] <- back_test[[i]]$Markets[2]
    #df[i,3] <- pop$`Gender = All | Ethnicity = All | Age = All`[which(pop$DMA==back_test[[i]]$markets[1])]
    #df[i,4] <- pop$`Gender = All | Ethnicity = All | Age = All`[which(pop$DMA==back_test[[i]]$markets[2])]

    results[i,5] <- back_test[[i]]$valid
    results[i,6] <- back_test[[i]]$error

    results[i,7] <- matches[i,"Correlation"]
    results[i,8] <- matches[i,"RelativeDistance"]

  }

  colnames(results) <- c("Control", "Exposed", "Control Pop", "Exposed Pop",
                         "Valid", "Train P-value", "Correlation", "Relative Distance")

  write.csv(results,paste0(getwd(),'/Market Pairs.csv'))

  return(list(back_test=back_test,results=results))

}

#' my title
#'
#' @name InterventionAnalysis
#'
#' @param data a data frame
#' @param metric conversion metric
#' @param pairs a list of market pairs
#' @param pre_period vector of two dates representing start and end of pre period
#' @param post_period vector of two dates representing start and end of post period
#' @param label label to match file path used for finding pairs
#' @param plots Boolean for plot inclusion
#'
#' @return list of valid tests
#' @importFrom CausalImpact CausalImpact
#' @importFrom utils combn
#' @importFrom zoo zoo
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @export

intervention <- function(data, metric, pairs, pre_period, post_period, label=NULL, plots=F) {
  impacts <- vector("list", length(pairs))
  if(plots == T) {
    setwd(paste0(getwd(), "/", label))
    dir.create(paste0(getwd(),"Intervention Plots"))
    for(j in 1:length(pairs)) {
      combs <- combn(pairs[[j]],2)
      if(ncol(combs)==1) {
        x1 <- as.data.frame(data[which(data[,"DMA"] == combs[1,1]),c("Date", metric)])
        y <- as.data.frame(data[which(data[,"DMA"] == combs[2,1]),c("Date", metric)])
        time.points <- x1[,1]

        data <- zoo(cbind(x1[,2], y[,2]), time.points)
        colnames(data) <- c('y','x')

        pre.period <- pre_period
        post.period <- post_period

        set.seed(500)
        impact <- CausalImpact(data, pre.period, post.period)

        png(paste0(getwd(),'Intervention Plots/',combs[1,1], ' - ', combs[2,1], '.png'), width=480*3*1.5,height=480*2*1.5)
        CreateImpactPlot(impact)
        x <- dev.off()

        impacts[[j]] <- impact

      } else {
        for (i in 1:ncol(combs)) {
          impacts2 <- vector('list', length(ncol(combs)))
          x1 <- as.data.frame(data[which(data[,"DMA"] == combs[1,i]),c("Date", metric)])
          y <- as.data.frame(data[which(data[,"DMA"] == combs[2,i]),c("Date", metric)])
          time.points <- x1[,1]

          data <- zoo(cbind(x1[,2], y[,2]), time.points)
          colnames(data) <- c('y','x')

          pre.period <- pre_period
          post.period <- post_period

          set.seed(500)
          impact <- CausalImpact(data, pre.period, post.period)

          png(paste0(getwd(),'/Intervention Plots/',combs[1,i], ' - ', combs[2,i], '.png'), width=480*3*1.5,height=480*2*1.5)
          plot(impact)
          x <- dev.off()

          impacts2[[i]] <- impact
        }
        impacts[[j]] <- impacts2
      }
    }

  } else {
    for(j in 1:length(pairs)) {
      combs <- combn(pairs[[j]],2)
      if(ncol(combs)==1) {
        x1 <- as.data.frame(data[which(data[,"DMA"] == combs[1,1]),c("Date", metric)])
        y <- as.data.frame(data[which(data[,"DMA"] == combs[2,1]),c("Date", metric)])
        time.points <- x1[,1]

        data <- zoo(cbind(x1[,2], y[,2]), time.points)
        colnames(data) <- c('y','x')

        pre.period <- pre_period
        post.period <- post_period

        set.seed(500)
        impact <- CausalImpact(data, pre.period, post.period)
        impacts[[j]] <- impact

      } else {
        for (i in 1:ncol(combs)) {
          impacts2 <- vector('list', length(ncol(combs)))
          x1 <- as.data.frame(data[which(data[,"DMA"] == combs[1,i]),c("Date", metric)])
          y <- as.data.frame(data[which(data[,"DMA"] == combs[2,i]),c("Date", metric)])
          time.points <- x1[,1]

          data <- zoo(cbind(x1[,2], y[,2]), time.points)
          colnames(data) <- c('y','x')

          pre.period <- pre_period
          post.period <- post_period

          set.seed(500)
          impact <- CausalImpact(data, pre.period, post.period)
          impacts2[[i]] <- impact
        }
        impacts[[j]] <- impacts2
      }
    }
  }
}

