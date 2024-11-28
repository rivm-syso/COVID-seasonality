##################################
### Functions for all analyses ###
##################################

# make confidence intervals for model parameters
simulate_intervals <- function(gamresult, toreturn = c("ah", "temp"), nrep = 100000) {
  grepterm <- paste0(toreturn, collapse = "|")
  intervalterms <- grep(grepterm, names(coef(gamresult)), value = TRUE)
  if(length(intervalterms) == 0) return(c())
  
  coefsample <- MASS::mvrnorm(nrep, 
                              mu = coef(gamresult)[intervalterms], 
                              Sigma = vcov(gamresult)[intervalterms,
                                                      intervalterms])
  return(apply(coefsample, 2, quantile, probs = c(.025, .5, .975)))
}

# create a readable confidence interval
present_interval <- function(intervals, whichparameter, precision = 3) {
  intervalpresentation <- c(0, " (", 0, ";", 0, ")")
  intervalpresentation[c(1,3,5)] <- round(intervals[, whichparameter], precision)[c(2,1,3)]
  intervalpresentation <- paste0(intervalpresentation, collapse = "")
  
  return(intervalpresentation)
}

# Create Rt predictions from a modelfit, with various options to include/exclude certain model terms. Used in plot_Rtpred below
simulate_Rt <- function(gamresult, newdata = NULL, with_season = TRUE, with_period = TRUE, with_offset = TRUE, nrep = 10000, intervals = TRUE) {
  if(with_season & with_period) {
    termstouse <- names(coef(gamresult))
  } else if(with_period) {
    termstouse <- grep("temp|ah", names(coef(gamresult)), value = TRUE, invert = TRUE)
  } else if(with_season) {
    termstouse <- c(grep("period", names(coef(gamresult)), value = TRUE, invert = TRUE))
  } else {
    termstouse <- c(grep("period|temp|ah", names(coef(gamresult)), value = TRUE, invert = TRUE))
  }
  
  if(is.null(newdata)) {
    toinclude <- rep(TRUE, length(gamresult$residuals))
    modmatall <- model.matrix(gamresult)
  } else {
    if(!with_period) {
      gamresult$coefficients["(Intercept)"] <- 0
    }
    toinclude <- paste0("period", newdata$period) %in% names(coef(gamresult)) | newdata$period == "_01"
    modmatall <- predict.gam(gamresult, newdata[toinclude, ], type = "lpmatrix")
  }
  
  modmatterms <- modmatall[, termstouse, drop = FALSE]     
  modmatoffset <- attr(modmatall, "model.offset")
  
    
  effectsample <- MASS::mvrnorm(nrep, 
                                mu = coef(gamresult)[termstouse], 
                                Sigma = vcov(gamresult)[termstouse,
                                                        termstouse])
  resultsample <- matrix(NA_real_, nrow = nrep, ncol = length(toinclude))
  
  resultsample[, toinclude] <- apply(modmatterms, 1, function(x) colSums(t(effectsample) * x)) + 
    with_offset * rep(modmatoffset, each = nrep)

  if(intervals) resultsample <- apply(resultsample, 2, quantile, probs = c(.025, .5, .975), na.rm = TRUE) %>%
    t() 
  
  return(resultsample)
}

# Create relative Rt predicted from a modelfit with only the weatherterms present in the model
simulate_season <- function(gamresult, newdata = NULL, nrep = 10000, intervals = TRUE) {
  seasonterms <- grep("temp|ah", names(coef(gamresult)), value = TRUE)
 
  if(length(seasonterms) == 0) {
    if(is.null(newdata)) {
      resultsample <- matrix(0, nrow = nrep, ncol = length(gamresult$residuals))
    } else {
      resultsample <- matrix(0, nrow = nrep, ncol = length(newdata$period))
    }
  } else {
    if(is.null(newdata)) {
      modmatahtemp <- model.matrix(gamresult)[, seasonterms, drop = FALSE]   
    } else {
      modmatahtemp <- predict.gam(gamresult, newdata = newdata %>%
                                    mutate(period = "_01"), type = "lpmatrix")[, seasonterms, drop = FALSE] 
    }
    
    effectsample <- matrix(MASS::mvrnorm(nrep, 
                                  mu = coef(gamresult)[seasonterms], 
                                  Sigma = vcov(gamresult)[seasonterms,
                                                          seasonterms]), ncol = length(seasonterms))
    resultsample <- matrix(apply(modmatahtemp, 1, function(x) colSums(t(effectsample) * x)), nrow = nrep)
  }
  
  if(intervals) resultsample <- apply(resultsample, 2, quantile, probs = c(.025, .5, .975)) %>% t()
  
  return(resultsample)
}

# Plot simulated Rt-prediction from modelfit, with options to include/exclude model terms
plot_Rtpred <- function(gamresult, observeddata,  
                              plottitle = "Rt and prediction confidence intervals",
                              plotobserved = FALSE, pointshape = 19,
                              with_period = TRUE, with_season = TRUE, with_offset = TRUE, coloroutliers = FALSE, yrange = c(0,5),
                        returndata = FALSE) {
  plotdata <- observeddata %>%
    dplyr::select(date, logRt) %>%
    bind_cols(matrix(simulate_Rt(gamresult, newdata = observeddata,
                                 with_period = with_period, with_season = with_season, with_offset = with_offset), ncol = 3),
              .name_repair = ~ c("date", "actual", "min", "est", "max")) %>%
    mutate(diff = est - actual,
           diffmin = min - actual,
           diffmax = max - actual,
           outside = diffmin > 0 | diffmax < 0) %>%
    mutate(across(c(actual, min:max), exp)) %>%
    mutate(logRttoplot = if(plotobserved) actual else est) 
  
  if(returndata) return(plotdata)
  
  if(coloroutliers) {
    toplot <- ggplot(aes(x = date, y = actual, color = outside), data = plotdata) +
      scale_color_manual(values = c("black", "red"), limits = c(FALSE, TRUE))
  } else (
    toplot <- ggplot(aes(x = date, y = logRttoplot), data = plotdata)
  )
  toplot +
    geom_point(shape = pointshape) +
    geom_linerange(aes(ymin = min, ymax = max))  +
    theme_light() +
    coord_cartesian(ylim = yrange) +
    labs(title = plottitle,
         x = "Date", y = "Reproduction number") 
  
}

# Plot the simulated relative Rt resulting from changing weather variables
plot_relRt_season <- function(gamresult, dataset, newdata = TRUE,
                              plottitle = "Seasonality of Rt (relative to mean temperature and humidity)",
                              xtitle = "Date") {
  if(newdata) {
    newdata <- dataset
  } else {
    newdata <- NULL
  }
  
  dataset %>% 
    dplyr::select(date) %>%
    bind_cols(
      as_tibble(simulate_season(gamresult, newdata))
    ) %>%
    mutate(across(`2.5%`:`97.5%`, exp)) %>%
    ggplot(aes(x = date, y = `50%`, ymin = `2.5%`, ymax = `97.5%`)) +
    geom_ribbon(fill = "grey") +
    geom_line() +
    theme_light() +
    labs(title = plottitle, x = xtitle, y = "Relative Rt ")
}

# Calculate the simulated relative Rt from a single model fit, or from multiple fits, weighted by Akaike Weights
summarise_relRt_season <- function(gamresults, Akaikeweights = NULL) {
  if(is.null(Akaikeweights)) {
    simres <- simulate_season(gamresults, 
                              newdata = alldata_weekly %>%
                                mutate(ah = ah_avg,
                                       temp = temp_avg,
                                       ah_1 = ah_avg_1,
                                       temp_1 = temp_avg_1,
                                       period = "_01") %>%
                                filter(year(date) == "2021") %>%
                                mutate(across(c(retail:residential, logRvariant:index_zorg, retail_1:index_zorg_1), mean)) %>%
                                dplyr::select(date, ah, temp, ah_1, temp_1, retail:index_zorg, retail_1:index_zorg_1) %>%
                                distinct %>%
                                arrange(date),
                              intervals = FALSE
    )
  } else {
    samplemodels <- rmultinom(1, 10000, Akaikeweights)[, 1]
    gamresults <- gamresults[samplemodels > 0]
    samplemodels <- samplemodels[samplemodels > 0]
    simres <- do.call(rbind, lapply(1:length(samplemodels),
                     function(x) simulate_season(
                       gamresults[[x]],
                       newdata = alldata_weekly %>%
                         mutate(ah = ah_avg,# - ah_reference,
                                temp = temp_avg,# - temp_reference,
                                ah_1 = ah_avg_1,# - ah_reference,
                                temp_1 = temp_avg_1,# - temp_reference, 
                                period = "_01") %>%
                         filter(year(date) == "2021") %>%
                         mutate(across(c(retail:residential, logRvariant:index_zorg, retail_1:index_zorg_1), mean)) %>%
                         dplyr::select(date, ah, temp, ah_1, temp_1, retail:index_zorg, retail_1:index_zorg_1) %>%
                         distinct %>%
                         arrange(date),
                       intervals = FALSE, nrep = samplemodels[x]
                     )))
  }
  simres <- exp(simres)
  yearlymins <- apply(simres, 1, min)
  yearlymaxs <- apply(simres, 1, max)
  yearlyamps <- apply(simres, 1, function(x) max(x)/min(x))
  
  toreturn <- rbind(min = quantile(yearlymins, c(.025, .5, .975)),
                   max = quantile(yearlymaxs, c(.025, .5, .975)),
                   amplitude = quantile(yearlyamps, c(.025, .5, .975)))
  return(toreturn)
}

# Calculate R-squared from a model fit
R_squared <- function(gamresult, adjusted = TRUE, excludeoffset = FALSE) {
  ### R-sq of lm subtracts offset from total sum-of-squares -> excludeoffset = TRUE
  ### R-sq of gam does not substract offset from total sos  -> excludeoffset = FALSE
  
  SSe <- sum(gamresult$weights * (gamresult$residuals)^2)
  obs <- gamresult$residuals + gamresult$fitted.values - excludeoffset * gamresult$offset
  SSt <- sum((gamresult$weights * (obs - weighted.mean(obs, gamresult$weights))^2))
  
  if(adjusted) {
    return(1 - (SSe/gamresult$df.residual)/(SSt/gamresult$df.null))
  } else {
    return(1 - SSe/SSt)
  }

}
